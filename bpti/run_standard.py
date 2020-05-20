import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
sys.path.append('../')
from gst.gst import grestIntegrator,SimulatedSoluteTempering
import numpy as np

####
###Setup
####
#after equilibration, standardMD and GST-MD will run for this many seconds:
number_ns = 200
one_ns = int(5e5)
number_steps = number_ns*one_ns
dcdstride = 50000

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=310
maxTemp=650
numTemps=6
exchangeInterval=500


def set_dihedral_force_group(system, g=2):
  """Sets the dihedral forcegroup to a number other than 0,
  which will be used by serial tempering"""
  print('Scanning forces:')
  for f in system.getForces():
    if isinstance(f, simtk.openmm.openmm.PeriodicTorsionForce):
      print('Found the torsions - setting group to 2')
      f.setForceGroup(2)
    print(f.getForceGroup(), f.__class__)

def setup_simulation(system, pdb, integrator):
  """Creates a simulation object"""
  #platform = Platform.getPlatformByName('CPU')
  platform = Platform.getPlatformByName('OpenCL')
  prop = {'OpenCLPrecision':'single'}
  
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy(100)
  simulation.context.setVelocitiesToTemperature(310*kelvin)
  print('Created simulation')
  return simulation


output_directory = './'

if __name__ == '__main__':

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')

  prmtop = AmberPrmtopFile('./bpti-ff99SBi.prmtop')
  inpcrd = AmberInpcrdFile('./bpti-ff99SBi.inpcrd')  
  system = prmtop.createSystem(nonbondedMethod=PME, constraints=HBonds,)

  system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))
  set_dihedral_force_group(system)
  
  integrator = grestIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)

  simulation = Simulation(prmtop.topology, system, integrator)
  simulation.context.setPositions(inpcrd.positions)
  simulation.minimizeEnergy(100)
  simulation.context.setVelocitiesToTemperature(310*kelvin)

  simulation.reporters.append(DCDReporter(output_directory+'bpti_standardmd.dcd', 25000))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,totalSteps=number_steps,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.step(number_steps)
  



