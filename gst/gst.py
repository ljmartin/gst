from __future__ import print_function

"""
Generalized serial tempering. This script is about 98% copy of Muneeb
Sultan's and Peter Eastman's simulatedtempering.py, which is open source.

A portion of the original license:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import simtk.unit as unit
from simtk.openmm import CustomIntegrator
import math
import random
from sys import stdout

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

class SimulatedSoluteTempering(object):
    """This script implements a serial version of generalized-REST.
    gREST: J. Chem. Phys. 149, 072304 (2018); https://doi.org/10.1063/1.5016222
    gREST is generalized in that it can apply tempering to any subset
    of the forces. This is simple to achieve in OpenMM using force groups.

    Arrange the forces such that the force to be tempered is the only force within
    it's group. Must be used in combination with a customintegrator that scales a force group.

    Only other differences to simulatedtempering.py are:
    1) the addition of a cutoff for the weight update factor, which is set to 1e-8 by
    default (the same value as in the original Wang-Landau paper).
    2) you must input the force group to be tempered.
    3) the output now includes the PE of the force group to aid in analysis.

    Otherwise see simualtedtempering.py for usage.
    """

    def __init__(self, simulation, forceGroupSet, cutoff, baseTemp, temperatures=None, numTemperatures=None, minTemperature=None, maxTemperature=None, weights=None, tempChangeInterval=25, reportInterval=1000, reportFile=stdout):
        """Create a new SimulatedTempering.

        Parameters
        ----------
        simulation: Simulation
            The Simulation defining the System, Context, and Integrator to use
        forceGroup: int
            Force group that will be tempered. Only tried with dihedrals so far.
        cutoff: float
            when to stop adjusting weights. When _weightUpdateFactor reduces below this,
            the _updateWeights flag turns off and weights no longer update (i.e. equilibration is over).
        temperatures: list
            The list of temperatures to use for tempering, in increasing order
        numTemperatures: int
            The number of temperatures to use for tempering.  If temperatures is not None, this is ignored.
        minTemperature: temperature
            The minimum temperature to use for tempering.  If temperatures is not None, this is ignored.
        maxTemperature: temperature
            The maximum temperature to use for tempering.  If temperatures is not None, this is ignored.
        weights: list
            The weight factor for each temperature.  If none, weights are selected automatically.
        tempChangeInterval: int
            The interval (in time steps) at which to attempt transitions between temperatures
        reportInterval: int
            The interval (in time steps) at which to write information to the report file
        reportFile: string or file
            The file to write reporting information to, specified as a file name or file object
        """
        self.forceGroupSet=forceGroupSet
        self.cutoff = cutoff
        print(self.cutoff)
        self.baseTemp = baseTemp
        self.simulation = simulation
        if temperatures is None:
            if unit.is_quantity(minTemperature):
                minTemperature = minTemperature.value_in_unit(unit.kelvin)
            if unit.is_quantity(maxTemperature):
                maxTemperature = maxTemperature.value_in_unit(unit.kelvin)
            self.temperatures = [minTemperature*((float(maxTemperature)/minTemperature)**(i/float(numTemperatures-1))) for i in range(numTemperatures)]*unit.kelvin
        else:
            numTemperatures = len(temperatures)
            self.temperatures = [(t.value_in_unit(unit.kelvin) if unit.is_quantity(t) else t)*unit.kelvin for t in temperatures]
            if any(self.temperatures[i] >= self.temperatures[i+1] for i in range(numTemperatures-1)):
                raise ValueError('The temperatures must be in strictly increasing order')
        self.tempChangeInterval = tempChangeInterval
        self.reportInterval = reportInterval
        for i in self.temperatures:
            print('Temperature:', i)

        print('Base temperature:', self.baseTemp)

        self.inverseTemperatures = [1.0/(unit.MOLAR_GAS_CONSTANT_R*t) for t in self.temperatures]
        print('inverseTemperatures:', self.inverseTemperatures)
        self.inverseBaseTemperature = 1.0/(unit.MOLAR_GAS_CONSTANT_R*baseTemp)
        print('inverseBaseTemperature:', self.inverseBaseTemperature)
        #setting all temperatures relative to lowest (not necessarily what you want - might want to simulate lower and higher than existing!
        #self.scalingFactors = [self.inverseTemperatures[i] / self.inverseTemperatures[0] for i in range(len(self.inverseTemperatures))]
        self.scalingFactors = [self.inverseTemperatures[i] / self.inverseBaseTemperature for i in range(len(self.inverseTemperatures))]
        print('scalingFactors:', self.scalingFactors)
        # If necessary, open the file we will write reports to.

        self._openedFile = isinstance(reportFile, str)
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if reportFile.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(reportFile, 'wb', 0))
            elif reportFile.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(reportFile, 'w', 0)
            else:
                self._out = open(reportFile, 'w', 1)
        else:
            self._out = reportFile

        # Initialize the weights.

        if weights is None:
            self._weights = [0.0]*numTemperatures
            self._updateWeights = True
            self._weightUpdateFactor = 1.0
            self._histogram = [0]*numTemperatures
            self._hasMadeTransition = False
        else:
            self._weights = weights
            self._updateWeights = False

        # Select the initial scaling factor.

        self.currentTemperature = 0
        #this was the old way, scaling by the integrator:
        #self.simulation.integrator.setScalingFactor(self.scalingFactors[self.currentTemperature])
        #this is the new way, scaling by parameter within the force:
        self.simulation.context.setParameter('scalingFactor', self.scalingFactors[self.currentTemperature])

        # Add a reporter to the simulation which will handle the updates and reports.

        class STReporter(object):
            def __init__(self, st, fg):
                self.st = st
                self.fg = fg

            def describeNextReport(self, simulation):
                st = self.st
                steps1 = st.tempChangeInterval - simulation.currentStep%st.tempChangeInterval
                steps2 = st.reportInterval - simulation.currentStep%st.reportInterval
                steps = min(steps1, steps2)
                isUpdateAttempt = (steps1 == steps)
                return (steps, False, False, False, False)

            def report(self, simulation, state):
                state = simulation.context.getState(getEnergy=True,getParameterDerivatives=True,groups=self.fg)
                st = self.st
                if st._weightUpdateFactor<st.cutoff:
                    st._updateWeights=False
                if simulation.currentStep%st.tempChangeInterval == 0:
                    st._attemptTemperatureChange(state)
                if simulation.currentStep%st.reportInterval == 0:
                    st._writeReport(state)


        simulation.reporters.append(STReporter(self, self.forceGroupSet))

        # Write out the header line.

        headers = ['Steps', 'Temperature (K)', 'PotentialEnergy', 'Scaling']
        for t in self.temperatures:
            headers.append('%gK Weight' % t.value_in_unit(unit.kelvin))
        print('"%s"' % ('"\t"').join(headers), file=self._out)
        #print('"'+('"\t"').join(headers)+'"', file=self._out)


    def __del__(self):
        if self._openedFile:
            self._out.close()


    @property
    def weights(self):
        return [x-self._weights[0] for x in self._weights]



    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self.simulation.step(steps)

    def _attemptTemperatureChange(self, state):
        """Attempt to move to a different temperature."""
        # Compute the probability for each temperature.  This is done in log space to avoid overflow.
        ##potential energy for a whole force group could be calculated like this:
        #nrg = -1*state.getPotentialEnergy()

        ##Instead, we want to account for situations where a subset of the whole
        ##force is calculated, where the subset is the interacitons that are
        ##scaled by scalingFactor. To do this, take the derivative wrt the scalingFactor,
        ##and the 'rise' of the derivative (being rise/run) is the contribution
        ##of just the scaled terms.
        nrg = state.getEnergyParameterDerivatives()['scalingFactor']/self.scalingFactors[self.currentTemperature]*unit.kilojoule/unit.mole

        logProbability = [(self._weights[i]-self.inverseTemperatures[i]*nrg) for i in range(len(self._weights))]
        maxLogProb = max(logProbability)
        offset = maxLogProb + math.log(sum(math.exp(x-maxLogProb) for x in logProbability))
        probability = [math.exp(x-offset) for x in logProbability]
        r = random.random()
        for j in range(len(probability)):
            if r < probability[j]:
                if j != self.currentTemperature:
                    ##### Rescale the velocities.
                    ##Velocities are not re-scaled in gREST -Lewis.

                    self._hasMadeTransition = True
                    self.currentTemperature = j
                    ###Instead, the scaling factor for the integrator is adjusted: -Lewis
                    #self.simulation.integrator.setScalingFactor(self.scalingFactors[j])
                    ###Instead, the global parameter 'scalingFactor' is adjusted:
                    self.simulation.context.setParameter('scalingFactor', self.scalingFactors[j])
                if self._updateWeights:
                    # Update the weight factors.

                    self._weights[j] -= self._weightUpdateFactor
                    self._histogram[j] += 1
                    minCounts = min(self._histogram)
                    if minCounts > 20 and minCounts >= 0.2*sum(self._histogram)/len(self._histogram):
                        # Reduce the weight update factor and reset the histogram.

                        self._weightUpdateFactor *= 0.5
                        self._histogram = [0]*len(self.temperatures)
                        self._weights = [x-self._weights[0] for x in self._weights]
                    elif not self._hasMadeTransition and probability[self.currentTemperature] > 0.99:
                        # Rapidly increase the weight update factor at the start of the simulation to find
                        # a reasonable starting value.

                        self._weightUpdateFactor *= 2.0
                        self._histogram = [0]*len(self.temperatures)
                return
            r -= probability[j]

    def _writeReport(self, state):
        """Write out a line to the report."""
        temperature = self.temperatures[self.currentTemperature].value_in_unit(unit.kelvin)
        #temperature = self.temperatures[self.currentTemperature]/unit.kelvin
        potEnergy = state.getEnergyParameterDerivatives()['scalingFactor']/self.scalingFactors[self.currentTemperature]
        #potEnergy = potEnergy.value_in_unit(unit.kilojoule/unit.mole)
        values = [temperature]+[potEnergy]+[self.scalingFactors[self.currentTemperature]]+self.weights
        print(('%d\t' % self.simulation.currentStep) + '\t'.join('%g' % v for v in values), file=self._out)
