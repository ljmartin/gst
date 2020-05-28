import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk import unit

from sys import stdout
sys.path.append('../')
#from gst.gst import grestIntegrator,SimulatedSoluteTempering
import numpy as np

####
###Setup
####
#after equilibration, standardMD and GST-MD will run for this many seconds:
number_ns = 50
one_ns = int(5e5)
number_steps = number_ns*one_ns

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=310
maxTemp=700
numTemps=7
exchangeInterval=500

#for metadynamics, you need to set the number of bins and the time to run for. 
num_grid_pts = 25
bias_factor = 3 #simulation will behave as if at temp*bias_factor
number_metad_steps = 200000


def setup_system_implicit(filename, barostat=False):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,)
  if barostat:
    system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))
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
  simulation = Simulation(pdb.topology, system, integrator, platform, prop)
  simulation.context.setPositions(pdb.positions)
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(310*kelvin)
  print('Created simulation')
  return simulation


filename ='./alanine-dipeptide-implicit.pdb'
output_directory = './'

if __name__ == '__main__':
  print('Running standard MD.')
  #########################
  # Run standardMD first. #
  #########################
  system, pdb = setup_system_implicit(filename)

  forces = { force.__class__.__name__ : force for force in system.getForces() }
  torsionforce = forces['PeriodicTorsionForce']
  #choose the solute atoms:
  solute = [int(i.index) for i in pdb.topology.atoms()]

  energy_expression = ""
  #energy_expression  = "select(condition, scalingFactor, 1)*all;"
  #energy_expression += "condition = soluteFlag;"
  energy_expression += "all;all=k*(1+cos(periodicity*theta-theta0))"
  #energy_expression += "all;all=0.5*k*(1-cos(theta-theta0))"
  custom_torsion_force = CustomTorsionForce(energy_expression)
  #custom_torsion_force.addGlobalParameter("scalingFactor", 1)
  #custom_torsion_force.addEnergyParameterDerivative('scalingFactor')
  #custom_torsion_force.addPerTorsionParameter("soluteFlag")
  custom_torsion_force.addPerTorsionParameter("k");
  custom_torsion_force.addPerTorsionParameter("periodicity")
  custom_torsion_force.addPerTorsionParameter("theta0");

  for torsion_index in range(torsionforce.getNumTorsions()):
    a0, a1, a2, a3, periodicity, phase, k = torsionforce.getTorsionParameters(torsion_index)
    #condition = (a0 in solute) | (a1 in solute) | (a2 in solute) | (a3 in solute)
    #soluteFlag = 1 if condition else 0
    #custom_torsion_force.addTorsion(a0, a1, a2, a3, [soluteFlag, k, periodicity,phase])
    custom_torsion_force.addTorsion(a0, a1, a2, a3, [k, periodicity,phase])
    #custom_torsion_force.addTorsion(a0, a1, a2, a3, [k,phase])

  for count, force in enumerate(system.getForces()):
    if isinstance(force, PeriodicTorsionForce):
      system.removeForce(count)
  system.addForce(custom_torsion_force)

  ###as a sanity check, print all the forces out with their
  ###force groups.
  #for force in system.getForces():
  #  if isinstance(force, CustomTorsionForce):
  #    force.setForceGroup(2)
  #  if isinstance(force, PeriodicTorsionForce):
  #    force.setForceGroup(2)
  #  print(force.__class__.__name__, force.getForceGroup(), force.__class__)



  cv1 = CustomTorsionForce('theta')
  cv1.addTorsion(1, 6, 8, 14)
  phi = BiasVariable(cv1, -np.pi, np.pi, 0.5, True,)
  cv2 = CustomTorsionForce('theta')
  cv2.addTorsion(6, 8, 14, 16)
  psi = BiasVariable(cv2, -np.pi, np.pi, 0.5,  True)

  meta = Metadynamics(system, [phi, psi], 310*kelvin, bias_factor, 1.0*kilojoules_per_mole, 50)


  integrator = LangevinIntegrator(310*kelvin, 91/picosecond, 0.002*picoseconds)
  simulation = setup_simulation(system, pdb, integrator)
  simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(310*kelvin)  
  meta.step(simulation, 350000)
  fes = meta.getFreeEnergy()
  np.save('fes.npy', fes)

#  #Instantiate reporters
#  simulation.reporters.append(DCDReporter(output_directory+'diala_standardmd_traj.dcd', 2500))
#  simulation.reporters.append(StateDataReporter(stdout, 5000, step=True,
#                                              potentialEnergy=True, density=True, speed=True))#
#
#  f = open('odd_dihedrals.dat', 'w', buffering=1)
#  for i in range(1000):
#    simulation.step(100)
#    state = simulation.context.getState(getEnergy=True,getParameterDerivatives=True,groups={2})
#    value = str(state.getPotentialEnergy()/unit.kilojoule*unit.mole)
#    f.write(value+'\n')
#    f.flush()
#    print(value)
  #simulation.step(number_steps)


