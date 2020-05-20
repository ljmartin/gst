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
number_ns = 1
one_ns = int(5e5)
number_steps = number_ns*one_ns
dcdstride = 50000


def setup_system_implicit(filename, barostat=False):
  """Creates a 'system' object given a pdb filename"""
  pdb = PDBFile(filename)
  forcefield = app.ForceField('amber99sb.xml', 'amber99_obc.xml')
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


filename = './TC5b_linear.pdb'
filename_implicit = './trp-min.pdb'
output_directory = './'

if __name__ == '__main__':
  print('Setting up an unfolded villin headpiece simulation')
  print('Saving minimized pdb structure')
  pdb = PDBFile(filename)
  #got to remove the waters and ions because we are running implicit:
  modeller = Modeller(pdb.topology, pdb.positions)
  to_delete = list()
  for atom in pdb.topology.atoms():
    if atom.residue.name in ['HOH', 'Cl']:
      to_delete.append(atom)
  modeller.delete(to_delete)

  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic,
                                   constraints=HBonds,implicitSolvent=OBC2,implicitSolventKappa=1.0/nanometer)
  integrator = grestIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)

  simulation = Simulation(modeller.topology, system, integrator)
  simulation.context.setPositions(modeller.positions)

  simulation.minimizeEnergy(50)
  
  PDBFile.writeFile(modeller.topology, 
                  simulation.context.getState(getPositions=True).getPositions(), 
                  open(filename_implicit, 'w'))
  
