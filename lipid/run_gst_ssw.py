import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
from simtk import unit

from sys import stdout
sys.path.append('../')
from gst.gst import grestIntegrator,SimulatedSoluteTempering
import numpy as np

####
###Setup
####
#after equilibration, standardMD and GST-MD will run for this many seconds:
number_ns = 100
##for 2fs timestep:
one_ns = int(5e5)
number_steps = number_ns*one_ns
dcdstride = 50000

##For 1fs timestep
#one_ns = int(1e6)
#number_steps = number_ns*one_ns
#dcdstride = 100000

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=310
maxTemp=650
numTemps=9
#2fs:
exchangeInterval=500
#1fs
#exchangeInterval=1000

def setup_system(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

  #box vectors from charmm-gui files:
  pdb.topology.setPeriodicBoxVectors((Vec3(5.75760367, 0.0, 0.0),
                                           Vec3(0, 5.75760367, 0.0),
                                           Vec3(0.0, 0.0, 6.0)))
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
  barostat = MonteCarloMembraneBarostat(1*bar, 200*bar*nanometer, 300*kelvin,
    MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree)
  system.addForce(barostat)
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
  simulation.context.setVelocitiesToTemperature(300*kelvin)
  print('Created simulation')
  return simulation

def create_cnb(original_nbforce):
    """Creates a CustomNonbondedForce object that recapitulates the function
    of the original nonbonded force"""
    #Determine PME parameters from nonbonded_force
    cutoff_distance = original_nbforce.getCutoffDistance()
    [alpha_ewald, nx, ny, nz] = original_nbforce.getPMEParameters()
    if (alpha_ewald/alpha_ewald.unit) == 0.0:
      # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance
      tol = original_nbforce.getEwaldErrorTolerance()
      alpha_ewald = (1.0/cutoff_distance) * np.sqrt(-np.log(2.0*tol))
    print(alpha_ewald)

    #Next, create a CustomNonbondedForce with LJ and Coulomb terms
    ONE_4PI_EPS0 = 138.935456
    energy_expression  = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*r)/r;"
    energy_expression += "epsilon = epsilon1*epsilon2;"
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
      [charge, sigma, epsilon] = original_nbforce.getParticleParameters(index)
      custom_nonbonded_force.addParticle([charge, sigma, epsilon])

    print('adding exceptions:')
    print('Number of exceptions before:', custom_nonbonded_force.getNumExclusions())
    print('Number of nb exceptions to add:', original_nbforce.getNumExceptions())
    for index in range(original_nbforce.getNumExceptions()):
      idx, jdx, c, s, eps = original_nbforce.getExceptionParameters(index)
      custom_nonbonded_force.addExclusion(idx, jdx)
    print('Number of exceptions after:', custom_nonbonded_force.getNumExclusions())

    return custom_nonbonded_force



filename = './step5_assembly.pdb'
output_directory = './'
topology.Topology.loadBondDefinitions('./dppc.xml')
topology.Topology.loadBondDefinitions('./chl1.xml')

if __name__ == '__main__':

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  #get an openmm system from the pdb file:
  system, pdb = setup_system(filename)

  ###now fetch the original nonbonded force:
  forces = { force.__class__.__name__ : force for force in system.getForces() }
  nbforce = forces['NonbondedForce']

  ###replicate the nonbondedforce object using CustomNonbondedForce's
  aqueous_cnb = create_cnb(nbforce)
  lipid_cnb = create_cnb(nbforce)

  ###now get indices of the three sets of atoms.
  chol = [int(i.index) for i in pdb.topology.atoms() if i.residue.name=='CHL1']
  dppc = [int(i.index) for i in pdb.topology.atoms() if not i.residue.name=='DPPC']
  solvent = [int(i.index) for i in pdb.topology.atoms() if i.residue.name not in ['DPPC', 'CHL1']]

  ###here's the complicated bit - adding interaction groups.
  ###We want water to interact with itself normally:
  aqueous_cnb.addInteractionGroup(solvent, solvent)
  ###We want water to interact with DPPC normally:
  aqueous_cnb.addInteractionGroup(solvent, dppc)
  ###We want water to interact with CHOL normally:
  aqueous_cnb.addInteractionGroup(solvent, chol)
  ###Now the final interaction between atom sets is the one we will bias:
  #now this uses lipid_cnb, NOT aqueous_cnb.
  lipid_cnb.addInteractionGroup(chol, dppc)

  ###Finally, delete the original nbforce,
  for count, f in enumerate(system.getForces()):
    if isinstance(f, simtk.openmm.openmm.NonbondedForce):
        system.removeForce(count)
  ###add the aqueous_cnb:
  system.addForce(aqueous_cnb)
  ###add the lipid_cnb:
  system.addForce(lipid_cnb)
  ###and set the forcegroup of the lipid_cnb to 2,
  ###so that we can scale it with simulated tempering.
  forces = [i for i in system.getForces()]
  forces[-1].setForceGroup(2)

  ###as a sanity check, print all the forces out with their
  ###force groups.
  for force in system.getForces():
    print(force.__class__.__name__, force.getForceGroup(), force.__class__)

  #carry on as usual:
  integrator = grestIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
  #integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
  simulation = setup_simulation(system, pdb, integrator)

  #simulation.reporters.append(DCDReporter(output_directory+'lipid_gst_equilibration.dcd', 1000))
  simulation.reporters.append(DCDReporter(output_directory+'lipid_gst_equilibration.dcd', 25000))

  #Just equilibrating the box first:
  for i in range(10):
    simulation.step(1000)
    print(i)
    print(simulation.context.getState(getEnergy=True, groups={2}).getPotentialEnergy())
    print(simulation.context.getState().getPeriodicBoxVectors())
  ###Then, equilibrate the weights:
  #The simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroup=2,
                              cutoff=1e-8,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=exchangeInterval,
                              reportFile='./lipid_gst_temp_equilibration.dat',
                             )

  stepsDone=0
  st.currentTemperature=numTemps-1
  while st._weightUpdateFactor>st.cutoff:
    simulation.step(exchangeInterval)
    stepsDone+=exchangeInterval
    print(st._weightUpdateFactor, st.currentTemperature, st.scalingFactors[st.currentTemperature])


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
                              reportInterval=dcdstride,
                              reportFile='./lipid_gst_temp.dat',
                             )
  #new st objects assume they need to equilibrate first. We don't.
  #so set the weights to what we determined, and turn off updating.
  st._weights = list(weights_store)
  st._updateWeights = False

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'lipid_gst_traj.dcd', dcdstride))
  simulation.reporters.append(StateDataReporter(stdout, 5000, step=True,totalSteps=number_steps+stepsDone,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.step(number_steps)
