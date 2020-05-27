import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
from simtk import unit

from sys import stdout
import numpy as np

sys.path.append('../')
from gst.gst import SimulatedSoluteTempering
from gst.utils import replace_nonbonded_force

####
###Setup
####
#after equilibration, standardMD and GST-MD will run for this many seconds:
number_ns = 1500
one_ns = int(5e5)
number_steps = number_ns*one_ns
dcdstride = 50000

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=310
maxTemp=700
numTemps=7
exchangeInterval=1000



def setup_system(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
  system.addForce(MonteCarloBarostat(1*bar, 323*kelvin))
  print('Created system')
  return system, pdb

def setup_simulation(system, pdb, integrator):
  """Creates a minimized simulation object"""
  #platform = Platform.getPlatformByName('CPU')
  platform = Platform.getPlatformByName('OpenCL')
  prop = {'OpenCLPrecision':'single'}

  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(323*kelvin)
  print('Created simulation')
  return simulation



filename = './input.pdb'
output_directory = './'

if __name__ == '__main__':

  #########################
  ## Run normal MD first. #
  #########################
  #system, pdb = setup_system(filename)
  #integrator = LangevinIntegrator(323*kelvin, 1/picosecond, 0.002*picoseconds)
  #simulation = setup_simulation(system, pdb, integrator)
  #simulation.reporters.append(DCDReporter(output_directory+'normal.dcd', 1000))
  #simulation.reporters.append(StateDataReporter(output_directory+'normal.dat',
  #                                  stride, step=True,totalSteps=number_steps,remainingTime=True,
  #                                  potentialEnergy=True, density=True, speed=True))
  #simulation.step(number_steps)


  ################################
  # Run reconstructed nonbonded  #
  ################################

  system, pdb = setup_system(filename)

  ###fetch the original nonbonded force:
  forces = { force.__class__.__name__ : force for force in system.getForces() }
  nbforce = forces['NonbondedForce']

  ###now get indices of the three sets of atoms.
  #solvent = [int(i.index) for i in pdb.topology.atoms() if i.residue.name in ['HOH', 'Cl']]
  solute = [int(i.index) for i in pdb.topology.atoms() if not (i.residue.name in ['HOH', 'Cl'])]

  ###replicate the nonbondedforce object using CustomNonbondedForces
  custom_nonbonded_force, custom_bond_force = replace_nonbonded_force(system, nbforce, solute)

  ###Finally, delete the original nbforce,
  for count, force in enumerate(system.getForces()):
    if isinstance(force, simtk.openmm.openmm.NonbondedForce):
        system.removeForce(count)

  ###add the aqueous_cnb:
  system.addForce(custom_nonbonded_force)
  system.addForce(custom_bond_force)

  ###as a sanity check, print all the forces out with their
  ###force groups.
  for force in system.getForces():
    if isinstance(force, CustomNonbondedForce):
      force.setForceGroup(2)
    if isinstance(force, CustomBondForce):
      force.setForceGroup(3)
    print(force.__class__.__name__, force.getForceGroup(), force.__class__)
    
  #carry on as usual:
  integrator = LangevinIntegrator(323*kelvin, 1/picosecond, 0.002*picoseconds)
  simulation = setup_simulation(system, pdb, integrator)

  simulation.reporters.append(DCDReporter(output_directory+'villin_gst_equilibration.dcd', 25000))

  ###First, equilibrate the weights:
  #The simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroupSet={2,3},
                              cutoff=1e-8,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=100,
                              reportFile='./villin_gst_temp_equilibration.dat',
                             )

  stepsDone=0
  st.currentTemperature=numTemps-1
  while st._weightUpdateFactor>st.cutoff:
    simulation.step(250)
    stepsDone+=250
    print(st._weightUpdateFactor, st.currentTemperature, st.weights)




  ###Then, run simulated tempering for real:
  print('Running GST')
  #remove the st reporter from the simulation object:
  print(simulation.reporters)
  simulation.reporters.pop(0)
  simulation.reporters.pop(0)
  print(simulation.reporters)
  weights_store = np.array(st.weights).copy()
  #the new simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroup=2,
                              cutoff=1e-8,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=500,
                              reportFile='./villin_gst_temp.dat',
                             )





  #new st objects assume they need to equilibrate first. We don't.
  #so set the weights to what we determined, and turn off updating.
  st._weights = list(weights_store)
  st._updateWeights = False

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'villin_gst_traj.dcd', dcdstride))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,totalSteps=number_steps+stepsDone,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.step(number_steps)

