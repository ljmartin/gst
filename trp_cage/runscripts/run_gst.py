import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
sys.path.append('../../')
from gst import utils
from gst.gst import SimulatedSoluteTempering
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
maxTemp=500
baseTemp = 310*kelvin
numTemps=11
exchangeInterval=500

def setup_system_implicit(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sb.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,
                                   implicitSolvent=HCT,implicitSolventKappa=1.0/nanometer)
  print('Created system')
  return system, pdb

def setup_simulation(system, pdb, integrator):
  """Creates a simulation object"""
  #platform = Platform.getPlatformByName('CPU')
  platform = Platform.getPlatformByName('OpenCL')
  prop = {'OpenCLPrecision':'single'}
  
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(baseTemp)
  print('Created simulation')
  return simulation


filename_implicit = '../trp-min.pdb'
output_directory = '../trajectories/'

if __name__ == '__main__':

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  system, pdb = setup_system_implicit(filename_implicit)
  integrator = LangevinIntegrator(baseTemp, 1/picosecond, 0.002*picoseconds)


  forces = { force.__class__.__name__ : force for force in system.getForces() }
  torsionforce = forces['PeriodicTorsionForce']
  #choose the solute atoms:
  solute = [int(i.index) for i in pdb.topology.atoms()] #kind of redundant with implicit

  custom_torsion_force = utils.replace_torsion_force(torsionforce, solute)
  ###Now, delete the original periodictorsionforce,
  for count, force in enumerate(system.getForces()):
    if isinstance(force, PeriodicTorsionForce):
        system.removeForce(count)
  ##And add the custom one:
  system.addForce(custom_torsion_force)

  for force in system.getForces():
    if isinstance(force, CustomTorsionForce):
      force.setForceGroup(2)
    print(force.__class__.__name__, force.getForceGroup(), force.__class__)

  simulation = setup_simulation(system, pdb, integrator)  
  #simulation.reporters.append(DCDReporter(output_directory+'trp_gst_equilibration.dcd', 25000))

  ###First, equilibrate the weights:
  #The simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroupSet={2},
                              cutoff=1e-8, 
                              baseTemp=baseTemp,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=dcdstride,
                              #reportFile=output_directory+'/trp_gst_temp_equilibration.dat',
                                reportFile=output_directory+'/trp_gst_temp.dat',
                             )

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'trp_gst_traj.dcd', dcdstride))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,totalSteps=number_steps,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  #this first "pre equilibration" might change the relative free energy of the states
  simulation.step(int(one_ns*20))

  #if that is the case, we re-set the equilibration after 20ns or so, allowing any new
  #folded conformation to have it's own equilibrated weights:
  #reset equilibration:
  print('Resetting weight equilibration:')
  st._updateWeights = True
  st._weightUpdateFactor = 1.0
  st._histogram = [0]*len(st.weights)
  st._hasMadeTransition = False
  simulation.step(int(180*one_ns))
  

  
    
#  stepsDone=0
#  st.currentTemperature=numTemps-1
#  while st._weightUpdateFactor>st.cutoff:
#    simulation.step(exchangeInterval)
#    stepsDone+=exchangeInterval
#    print(st._weightUpdateFactor, st.currentTemperature)
#    
#
#  ###Then, run simulated tempering for real:
#  print('Running GST')
#  #remove the st reporter from the simulation object:
#  print(simulation.reporters)
#  simulation.reporters.pop(0)
#  simulation.reporters.pop(0)
#  print(simulation.reporters)
#  weights_store = np.array(st.weights).copy()
#  #the new simulated tempering object:
#  st = SimulatedSoluteTempering(simulation,
#                              forceGroupSet={2},
#                              cutoff=1e-8, #now back to typical size. 
#                              baseTemp=baseTemp,
#                              numTemperatures=numTemps,
#                              tempChangeInterval=exchangeInterval,
#                              minTemperature=minTemp*kelvin,
#                              maxTemperature=maxTemp*kelvin,
#                              reportInterval=dcdstride,
#                              reportFile=output_directory+'trp_gst_temp.dat',
#                             )
#  #new 'st' objects assume they need to equilibrate first. If the starting state
#  #was close to equilibrium, we might now set _updateWeights to 'False'.
#  #Instead, we have weights equilibrated to the (potentially) unfolded state.
#  #So now we will let it have another period of equilibration - hopefully
#  #it's pretty close and gets there quickly. 
#  st._weights = list(weights_store)
#  st._updateWeights = True
#
#  ###Before running, we re-set the positions and velocities for a fair comparison to the standardMD
#  simulation.context.setPositions(pdb.positions)
#  simulation.minimizeEnergy()
#  simulation.context.setVelocitiesToTemperature(baseTemp)
#  
#  #now add the normal reporters and run!
#  simulation.reporters.append(DCDReporter(output_directory+'trp_gst_traj.dcd', dcdstride))
#  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,totalSteps=number_steps+stepsDone,remainingTime=True,
#                                              potentialEnergy=True, density=True, speed=True))
#  simulation.step(number_steps)
  



