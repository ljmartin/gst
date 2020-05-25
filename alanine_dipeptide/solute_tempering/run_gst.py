import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk import unit

from sys import stdout
sys.path.append('../../')
from gst.gst import grestIntegrator,SimulatedSoluteTempering
import numpy as np

####
###Setup
####
#after equilibration, standardMD and GST-MD will run for this many seconds:
number_ns = 500
one_ns = int(5e5)
number_steps = number_ns*one_ns

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=300
maxTemp=650
numTemps=6
exchangeInterval=500

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
  system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
  print('Created system')
  return system, pdb

def setup_simulation(system, pdb, integrator):
  """Creates a simulation object"""
  #platform = Platform.getPlatformByName('OpenCL')
  platform = Platform.getPlatformByName('OpenCL')
  prop = {'OpenCLPrecision':'single'}
  #prop = {'CudaPrecision':'single'}
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.context.setPeriodicBoxVectors(Vec3(x=2.8651899104606384, y=0.0, z=0.0),
                                           Vec3(x=0.0, y=2.865974822316307, z=0.0),
                                           Vec3(x=0.0, y=0.0, z=2.7781519069102902))
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(300*kelvin)
  print('Created simulation')
  return simulation


filename ='../alanine-dipeptide-explicit.pdb'
#coords_filename = './alanine-dipeptide-implicit.inpcrd'
output_directory = './'

if __name__ == '__main__':
  print('Running standard MD.')
  #########################
  # Run standardMD first. #
  #########################
  system, pdb = setup_system_explicit(filename)

  ###setting up the custom nonbonded force:
  forces = { force.__class__.__name__ : force for force in system.getForces() }
  nbforce = forces['NonbondedForce']

  solute = [int(i.index) for i in pdb.topology.atoms() if i.residue.name !='HOH']
  print('solute:', solute)
  solvent = [int(i.index) for i in pdb.topology.atoms() if i.residue.name =='HOH']
  print('solvent:', solvent[:10], solvent[-10:])


  # Determine PME parameters from nonbonded_force
  cutoff_distance = nbforce.getCutoffDistance()
  [alpha_ewald, nx, ny, nz] = nbforce.getPMEParameters()
  if (alpha_ewald/alpha_ewald.unit) == 0.0:
    # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance
    tol = nbforce.getEwaldErrorTolerance()
    alpha_ewald = (1.0/cutoff_distance) * np.sqrt(-np.log(2.0*tol))        
  print(alpha_ewald)

  ONE_4PI_EPS0 = 138.935456
  energy_expression  = "(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*r)/r);"
#  energy_expression  = "(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r);"
  energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
  energy_expression += "sigma = 0.5*(sigma1+sigma2);"
  energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)  # already in OpenMM units
  energy_expression += "chargeprod = charge1*charge2;"
  energy_expression += "alpha_ewald = {:f};".format(alpha_ewald.value_in_unit_system(unit.md_unit_system))
  custom_nonbonded_force = CustomNonbondedForce(energy_expression)
  custom_nonbonded_force.addPerParticleParameter('charge')
  custom_nonbonded_force.addPerParticleParameter('sigma')
  custom_nonbonded_force.addPerParticleParameter('epsilon')

  # Configure force
  custom_nonbonded_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
  custom_nonbonded_force.setCutoffDistance(1*nanometer)

  print('adding particles to custom force')
  for index in range(system.getNumParticles()):
    [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
    custom_nonbonded_force.addParticle([charge, sigma, epsilon])
    if index in solute:
      nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)

  print('adding exceptions:')#custom_nonbonded_force.getNumExclusions()
  print('Number of exceptions:', custom_nonbonded_force.getNumExclusions())
  for index in range(nbforce.getNumExceptions()):
    idx, jdx, c, s, eps = nbforce.getExceptionParameters(index)
    custom_nonbonded_force.addExclusion(idx, jdx)
  print('Number of exceptions:', custom_nonbonded_force.getNumExclusions())

  custom_nonbonded_force.addInteractionGroup(solute, solute)
  custom_nonbonded_force.addInteractionGroup(solvent, solute)

  system.addForce(custom_nonbonded_force)





  integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
  #integrator = grestIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
  simulation = setup_simulation(system, pdb, integrator)
  
  #Instantiate reporters
  simulation.reporters.append(DCDReporter(output_directory+'diala_standardmd_traj.dcd', 25000))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True, remainingTime=True,totalSteps=number_steps))


  for i in range(10000):
    simulation.step(5000)
    print(i)
    print(simulation.context.getState(getEnergy=True, groups={2}).getPotentialEnergy())
    print(simulation.context.getState().getPeriodicBoxVectors())










  simulation.step(number_steps)








  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  system, pdb = setup_system_implicit(filename)
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
                              reportInterval=50,
                              reportFile='./diala_gst_temp_equilibration.dat',
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
                              reportInterval=500,
                              reportFile='./diala_gst_temp.dat',
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
  simulation.reporters.append(DCDReporter(output_directory+'diala_gst_traj.dcd', 25000))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,
                                                potentialEnergy=True, density=True, speed=True,totalSteps=number_steps+stepsDone, remainingTime=True))
  simulation.step(number_steps)
  


  ######################
  # Run metadynamics.  #
  ######################
  print('Running metadynamics')
  system, pdb = setup_system_implicit(filename)
  integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

  cv1 = CustomTorsionForce('theta')
  cv1.addTorsion(1, 6, 8, 14)
  phi = BiasVariable(cv1, -np.pi, np.pi, 0.5, True,gridWidth=num_grid_pts,)
  cv2 = CustomTorsionForce('theta')
  cv2.addTorsion(6, 8, 14, 16)
  psi = BiasVariable(cv2, -np.pi, np.pi, 0.5,  True, gridWidth=num_grid_pts)

  meta = Metadynamics(system, [phi, psi], 300*kelvin, bias_factor, 1.0*kilojoules_per_mole, 50)
  simulation = setup_simulation(system, pdb, integrator)

  simulation.reporters.append(DCDReporter('./diala_metad_traj.dcd', 500))
  simulation.reporters.append(StateDataReporter(stdout, 10000, step=True,
                                                potentialEnergy=True, density=True, speed=True, totalSteps=number_metad_steps, remainingTime=True))
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(310*kelvin)

  meta.step(simulation, number_metad_steps)
