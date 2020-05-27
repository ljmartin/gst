import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
from simtk import unit

import numpy as np

def replace_nonbonded_force(system, original_nbforce, solute):
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
    energy_expression  = "scalingFactor*all;"
    energy_expression += "all=4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r;"
    energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)  # already in OpenMM units
    custom_bond_force = openmm.CustomBondForce(energy_expression)
    custom_bond_force.addGlobalParameter('scalingFactor', 1)
    custom_bond_force.addEnergyParameterDerivative('scalingFactor')
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
        c_value = c/unit.elementary_charge**2
        eps_value = eps/(unit.kilojoule/unit.mole)
        if c_value != 0 or eps_value!=0:
          custom_bond_force.addBond(idx, jdx, [c, s, eps])
    print('Number of exceptions after:', custom_nonbonded_force.getNumExclusions())

    return custom_nonbonded_force, custom_bond_force
