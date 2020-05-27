import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
from simtk import unit

from sys import stdout
import numpy as np

####
###Setup
####
stride = 1000
number_steps = 1000000



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

def create_cnb(original_nbforce, solute):
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
    energy_expression  = "select(condition, scalingFactor, 1)*all;"
    energy_expression += "condition = soluteFlag2*soluteFlag2;" #solute must have flag int(1)
    energy_expression += "all=4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*r)/r;"
    energy_expression += "epsilon = epsilon1*epsilon2;"
    energy_expression += "sigma = 0.5*(sigma1+sigma2);"
    energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)  # already in OpenMM units
    energy_expression += "chargeprod = charge1*charge2;"
    energy_expression += "alpha_ewald = {:f};".format(alpha_ewald.value_in_unit_system(unit.md_unit_system))
    custom_nonbonded_force = CustomNonbondedForce(energy_expression)
    custom_nonbonded_force.addGlobalParameter('scalingFactor', 1)
    custom_nonbonded_force.addEnergyParameterDerivative('scalingFactor')
    custom_nonbonded_force.addPerParticleParameter('soluteFlag')
    custom_nonbonded_force.addPerParticleParameter('charge')
    custom_nonbonded_force.addPerParticleParameter('sigma')
    custom_nonbonded_force.addPerParticleParameter('epsilon')
    # Configure force
    custom_nonbonded_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom_nonbonded_force.setCutoffDistance(1*nanometer)
    custom_nonbonded_force.setUseLongRangeCorrection(False)
    custom_nonbonded_force.setUseSwitchingFunction(True)
    custom_nonbonded_force.setSwitchingDistance(cutoff_distance - 1.0*unit.angstroms)

    #Now we must add a bond force to handle exceptions.
    #(exceptions are altered interactions - but not completely
    #removed! Also PME is not required now.
    energy_expression  = "(4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r);"
    energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)  # already in OpenMM units
    custom_bond_force = openmm.CustomBondForce(energy_expression)
    custom_bond_force.addPerBondParameter('chargeprod')
    custom_bond_force.addPerBondParameter('sigma')
    custom_bond_force.addPerBondParameter('epsilon')

    print('adding particles to custom force')
    for index in range(system.getNumParticles()):
      [charge, sigma, epsilon] = original_nbforce.getParticleParameters(index)
      soluteFlag = 1 if index in solute else 0
      custom_nonbonded_force.addParticle([soluteFlag, charge, sigma, epsilon])

    print('adding exceptions:')
    print('Number of exceptions before:', custom_nonbonded_force.getNumExclusions())
    print('Number of nb exceptions to add:', original_nbforce.getNumExceptions())
    for index in range(original_nbforce.getNumExceptions()):
      idx, jdx, c, s, eps = original_nbforce.getExceptionParameters(index)
      custom_nonbonded_force.addExclusion(idx, jdx)

      #solute has some altered (not removed) exceptions,
      #so these are added to the bond_force:
      if (idx in solute) and (jdx in solute):
#        c_value = c/unit.elementary_charge**2
#        eps_value = eps/(unit.kilojoule/unit.mole)
#        if c_value != 0 or eps_value!=0:
        custom_bond_force.addBond(idx, jdx, [c, s, eps])
    print('Number of exceptions after:', custom_nonbonded_force.getNumExclusions())

    return custom_nonbonded_force, custom_bond_force


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
  all_indices = [int(i.index) for i in pdb.topology.atoms()]
  solvent = [int(i.index) for i in pdb.topology.atoms() if i.residue.name in ['HOH', 'Cl']]
  solute = [int(i.index) for i in pdb.topology.atoms() if not (i.residue.name in ['HOH', 'Cl'])]

  ###replicate the nonbondedforce object using CustomNonbondedForces
  custom_nonbonded_force, custom_bond_force = create_cnb(nbforce, solute)

  ####We want water to interact with itself normally:
  #custom_nonbonded_force.addInteractionGroup(solvent, solvent)
  #####We want protein to interact with itself normally:
  #custom_nonbonded_force.addInteractionGroup(solute, solute)
  #####We want solvent and solute to interact normally:
  #custom_nonbonded_force.addInteractionGroup(solute, solvent)

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
    print(force.__class__.__name__, force.getForceGroup(), force.__class__)
    
  #carry on as usual:
  integrator = LangevinIntegrator(323*kelvin, 1/picosecond, 0.002*picoseconds)
  simulation = setup_simulation(system, pdb, integrator)
  simulation.reporters.append(DCDReporter(output_directory+'reconstructed.dcd', stride))
  simulation.reporters.append(StateDataReporter(output_directory+'reconstructed.dat',
                                  stride, step=True,totalSteps=number_steps,remainingTime=True,
                                  potentialEnergy=True, density=True, speed=True))

  f = open('scaled_PEs.dat', 'w', buffering=1)
  for scalingFactor in [1,0.9,0.8,0.7]:
    simulation.context.setParameter('scalingFactor', scalingFactor)
    for i in range(50):
      simulation.step(1000)
      state = simulation.context.getState(getParameterDerivatives=True, getEnergy=True, groups={2})
      partial_energy = state.getEnergyParameterDerivatives()['scalingFactor']/scalingFactor
      total_energy = state.getPotentialEnergy()
      f.write(str(scalingFactor)+'\t'+str(partial_energy)+'\t'+str(total_energy)+'\n')
      
    
#  for i in range(int(number_steps/1000)):
#    simulation.step(1000)

#    print(state.getEnergyParameterDerivatives()['scalingFactor']/0.75) 
    #print(simulation.context.getState(getParameterDerivatives=True, groups=2**1).getEnergyParameterDerivatives()['lambda_vdw']]
