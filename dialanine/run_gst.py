import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
sys.path.append('../')
from gst.gst import grestIntegrator,SimulatedSoluteTempering

def setup_system(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds) 
  #barostat = MonteCarloBarostat(1*bar, 310*kelvin)
  #system.addForce(barostat)
  set_dihedral_force_group(system)
  print('Created system')
  return system, pdb

def setup_system_vacuum(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
  system = forcefield.createSystem(pdb.topology, constraints=AllBonds, hydrogenMass=3*amu)
  #barostat = MonteCarloBarostat(1*bar, 310*kelvin)
  #system.addForce(barostat)
  set_dihedral_force_group(system)
  print('Created system')
  return system, pdb

def setup_system_implicit(filename):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,)
  #barostat = MonteCarloBarostat(1*bar, 310*kelvin)
  #system.addForce(barostat)
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
  platform = Platform.getPlatformByName('CPU')
  #prop = {'CudaPrecision':'single'}
  simulation = Simulation(pdb.topology, system, integrator, platform)
  simulation.context.setPositions(pdb.positions)
  simulation.context.setPeriodicBoxVectors(Vec3(2.846979248047194, 0.0, 0.0), 
                                           Vec3(-1.423489624023597, 2.8460968399199644, 0.0), 
                                           Vec3(0.0, 0.0, 2.767338609857905))
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(310*kelvin)
  print('Created simulation')
  return simulation




filename ='./alanine-dipeptide-implicit.pdb'
output_directory = './'

system, pdb = setup_system_implicit(filename)
integrator = grestIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
simulation = setup_simulation(system, pdb, integrator)

###Instantiate reporters
simulation.reporters.append(DCDReporter(output_directory+'diala_gst_traj2.dcd', 2500))

##The simulated tempering object:
st = SimulatedSoluteTempering(simulation,
                              forceGroup=2,
                              cutoff=1e-8,
                              numTemperatures=15,
                              tempChangeInterval=250,
                              minTemperature=310*kelvin,
                              maxTemperature=1500*kelvin,
                              reportInterval=500,
                              reportFile=output_directory+'diala_gst_temp2.dat',
                             )

one_ns = int(5e5)
number_ns = 100
number_steps = number_ns*one_ns
simulation.reporters.append(StateDataReporter(stdout, 2500, step=True,
        potentialEnergy=True, temperature=True,progress=True, totalSteps=number_steps,remainingTime=True))

simulation.step(number_steps)
