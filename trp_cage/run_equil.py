import simtk
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
sys.path.append('../')
from gst.gst import grestIntegrator,SimulatedSoluteTempering
import numpy as np


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
  
