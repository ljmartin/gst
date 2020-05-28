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
number_ns = 500
one_ns = int(5e5)
number_steps = number_ns*one_ns

#for simulated tempering, you need min temp, max temp, and number of temps total
baseTemp = 298*kelvin
minTemp=298
maxTemp=650
numTemps=3
exchangeInterval=250

#for metadynamics, you need to set the number of bins and the time to run for. 
num_grid_pts = 25
bias_factor = 3 #simulation will behave as if at temp*bias_factor
number_metad_steps = 200000


def setup_system_implicit(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,)

  #"""prmtop and crd are from https://github.com/choderalab/gibbs/blob/master/openmm/python/parallel-tempering.py"""
  #prmtop = AmberPrmtopFile(filename)
  #system = prmtop.createSystem(nonbondedMethod=app.CutoffNonPeriodic, constraints=app.HBonds, implicitSolvent=app.OBC2)
  
  set_dihedral_force_group(system)
  print('Created system')
  return system, pdb

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
  #platform = Platform.getPlatformByName('CUDA') 
  #prop = {'CudaPrecision':'single'}
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(baseTemp)
  print('Created simulation')
  return simulation


filename ='../alanine-dipeptide-implicit.pdb'
output_directory = '../trajectories'

if __name__ == '__main__':
  print('Running standard MD.')
  #########################
  # Run standardMD first. #
  #########################
  system, pdb = setup_system_implicit(filename)
  
  integrator = LangevinIntegrator(baseTemp, 91/picosecond, 0.002*picoseconds)
  simulation = setup_simulation(system, pdb, integrator)
  
  #Instantiate reporters
  simulation.reporters.append(DCDReporter(output_directory+'diala_standardmd_traj.dcd', 25000))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True, remainingTime=True,totalSteps=number_steps))
  simulation.step(number_steps)

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  system, pdb = setup_system_implicit(filename)
  integrator = LangevinIntegrator(baseTemp, 91/picosecond, 0.002*picoseconds)

  forces = { force.__class__.__name__ : force for force in system.getForces() }
  torsionforce = forces['PeriodicTorsionForce']
  #choose the solute atoms:
  solute = [int(i.index) for i in pdb.topology.atoms()] #kind of redundant with implicit
  
  custom_torsion_force = utils.replace_torsion_force(torsionforce, solute)
  ###Finally, delete the original nbforce,
  for count, force in enumerate(system.getForces()):
    if isinstance(force, simtk.openmm.openmm.NonbondedForce):
        system.removeForce(count)
  ##And add the custom one:
  system.addForce(custom_torsion_force)

  for force in system.getForces():
    if isinstance(force, CustomTorsionForce):
      force.setForceGroup(2)
    #if isinstance(force, PeriodicTorsionForce):
    #  force.setForceGroup(2)
    print(force.__class__.__name__, force.getForceGroup(), force.__class__)
  
  simulation = setup_simulation(system, pdb, integrator)

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
                              reportInterval=50,
                              reportFile=output_directory+'diala_gst_temp_equilibration.dat',
                             )


  stepsDone=0
  while st._weightUpdateFactor>st.cutoff:
    simulation.step(exchangeInterval)
    stepsDone+=exchangeInterval
    print(st._weightUpdateFactor, st.currentTemperature)
    

  ###Then, run simulated tempering for real:
  print('Running GST')
  #remove the st reporter from the simulation object:
  simulation.reporters.pop(0)
  weights_store = np.array(st.weights).copy()
  #the new simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroupSet={2},
                              cutoff=1e-8,
                              baseTemp=baseTemp,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=exchangeInterval,
                              reportFile=output_directory+'diala_gst_temp.dat',
                             )
  #new st objects assume they need to equilibrate first. We don't.
  #so set the weights to what we determined, and turn off updating.
  st._weights = list(weights_store)
  st._updateWeights = False
  
  ##With equilibrated weights, we can now reset the simulation to have a fair go at comparison to standardMD
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(baseTemp)

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'diala_gst_traj.dcd', 25000))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True,totalSteps=number_steps+stepsDone, remainingTime=True))
  simulation.step(number_steps)
  


  ######################
  # Run metadynamics.  #
  ######################
  print('Running metadynamics')
  system, pdb = setup_system_implicit(filename)
  integrator = LangevinIntegrator(baseTemp, 1/picosecond, 0.002*picoseconds)

  cv1 = CustomTorsionForce('theta')
  cv1.addTorsion(1, 6, 8, 14)
  phi = BiasVariable(cv1, -np.pi, np.pi, 0.5, True,gridWidth=num_grid_pts,)
  cv2 = CustomTorsionForce('theta')
  cv2.addTorsion(6, 8, 14, 16)
  psi = BiasVariable(cv2, -np.pi, np.pi, 0.5,  True, gridWidth=num_grid_pts)

  meta = Metadynamics(system, [phi, psi], baseTemp, bias_factor, 1.0*kilojoules_per_mole, 50)
  simulation = setup_simulation(system, pdb, integrator)

  simulation.reporters.append(DCDReporter('./diala_metad_traj.dcd', 500))
  simulation.reporters.append(StateDataReporter(stdout, 10000, step=True,
                                                potentialEnergy=True, density=True, speed=True, totalSteps=number_metad_steps, remainingTime=True))
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(baseTemp)

  meta.step(simulation, number_metad_steps)
