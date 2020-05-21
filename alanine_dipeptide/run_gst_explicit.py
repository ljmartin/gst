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

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=300
maxTemp=600
numTemps=6
exchangeInterval=250

#for metadynamics, you need to set the number of bins and the time to run for. 
num_grid_pts = 25
bias_factor = 3 #simulation will behave as if at temp*bias_factor
number_metad_steps = 200000


def setup_system_explicit(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  #forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
  forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds) 
  
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
  #platform = Platform.getPlatformByName('OpenCL')
  platform = Platform.getPlatformByName('CUDA')
  #prop = {'OpenCLPrecision':'single'}
  prop = {'CudaPrecision':'single'}
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.context.setPeriodicBoxVectors(Vec3(x=2.8651899104606384, y=0.0, z=0.0),
                                           Vec3(x=0.0, y=2.865974822316307, z=0.0),
                                           Vec3(x=0.0, y=0.0, z=2.7781519069102902))
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(300*kelvin)
  print('Created simulation')
  return simulation


filename ='./alanine-dipeptide-explicit.pdb'
#filename ='./alanine-dipeptide-implicit.prmtop'
#coords_filename = './alanine-dipeptide-implicit.inpcrd'
output_directory = './explicit/'

if __name__ == '__main__':
  print('Running standard MD.')
  #########################
  # Run standardMD first. #
  #########################
  system, pdb = setup_system_explicit(filename)
  #inpcrd = AmberInpcrdFile(coords_filename)
  
  integrator = grestIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
  simulation = setup_simulation(system, pdb, integrator)
  
  #Instantiate reporters
  simulation.reporters.append(DCDReporter(output_directory+'diala_standardmd_traj.dcd', 2500))
  simulation.reporters.append(StateDataReporter(output_directory+'diala_standardmd.out', 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True, remainingTime=True,totalSteps=number_steps))
  simulation.step(number_steps)

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  system, pdb = setup_system_explicit(filename)
  integrator = grestIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)

  simulation = setup_simulation(system, pdb, integrator)

  ###First, equilibrate the weights:
  #The simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroup=2,
                              cutoff=1e-4,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=exchangeInterval,
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
                              forceGroup=2,
                              cutoff=1e-8,
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
  simulation.context.setVelocitiesToTemperature(300*kelvin)

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'diala_gst_traj.dcd', 2500))
  simulation.reporters.append(StateDataReporter(output_directory+'diala_gst.out', 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True,totalSteps=number_steps+stepsDone, remainingTime=True))
  simulation.step(number_steps)
  


  ######################
  # Run metadynamics.  #
  ######################
  print('Running metadynamics')
  system, pdb = setup_system_explicit(filename)
  integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

  cv1 = CustomTorsionForce('theta')
  cv1.addTorsion(1, 6, 8, 14)
  phi = BiasVariable(cv1, -np.pi, np.pi, 0.5, True,gridWidth=num_grid_pts,)
  cv2 = CustomTorsionForce('theta')
  cv2.addTorsion(6, 8, 14, 16)
  psi = BiasVariable(cv2, -np.pi, np.pi, 0.5,  True, gridWidth=num_grid_pts)

  meta = Metadynamics(system, [phi, psi], 300*kelvin, bias_factor, 1.0*kilojoules_per_mole, 50)
  simulation = setup_simulation(system, pdb, integrator)

  simulation.reporters.append(DCDReporter(output_directory+'diala_metad_traj.dcd', 500))
  simulation.reporters.append(StateDataReporter(output_directory+'diala_metad.out', 10000, step=True,
                                                potentialEnergy=True, density=True, speed=True, totalSteps=number_metad_steps, remainingTime=True))
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(310*kelvin)

  meta.step(simulation, number_metad_steps)
