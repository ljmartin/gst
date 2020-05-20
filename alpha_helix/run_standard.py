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
number_ns = 1500
one_ns = int(5e5)
number_steps = number_ns*one_ns
dcdstride = 50000

def setup_system_implicit(filename, barostat=False):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,
                                   implicitSolvent=OBC2,implicitSolventKappa=1.0/nanometer)
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


filename = './input.pdb'
filename_implicit = './step1_pdbreader.pdb'
output_directory = './'

if __name__ == '__main__':
  print('Running standard MD.')
  #########################
  # Run standardMD first. #
  #########################
  system, pdb = setup_system_implicit(filename_implicit)
  integrator = grestIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
  simulation = setup_simulation(system, pdb, integrator)
  
  #Instantiate reporters
  simulation.reporters.append(DCDReporter(output_directory+'villin_standardmd_traj.dcd', dcdstride))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True, totalSteps=number_steps,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.step(number_steps)

  
