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

#for simulated tempering, you need min temp, max temp, and number of temps total
minTemp=310
maxTemp=600
numTemps=6
exchangeInterval=250

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


filename = './input.pdb'
filename_implicit = './villin-implicit.pdb'
output_directory = './'

if __name__ == '__main__':
  print('Setting up an unfolded villin headpiece simulation')
  print('Running at 650K')
  pdb = PDBFile(filename)
  #got to remove the waters and ions because we are running implicit:
  modeller = Modeller(pdb.topology, pdb.positions)
  to_delete = list()
  for atom in pdb.topology.atoms():
    if atom.residue.name in ['HOH', 'Cl']:
      to_delete.append(atom)
  modeller.delete(to_delete)

  forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
  system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.CutoffNonPeriodic, constraints=HBonds,)
  integrator = grestIntegrator(650*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)

  simulation = Simulation(modeller.topology, system, integrator)
  simulation.context.setPositions(modeller.positions)

  simulation.minimizeEnergy(50)
  simulation.context.setVelocitiesToTemperature(650*kelvin)
  simulation.reporters.append(DCDReporter('./villin/villin_equil.dcd', 500)) 
  simulation.reporters.append(StateDataReporter(stdout, 5000, step=True,
        potentialEnergy=True, density=True, speed=True))

  simulation.step(500000)
  PDBFile.writeFile(modeller.topology, modeller.positions, open(filename_implicit, 'w'))

  ######################
  # Run generalizedST. #
  ######################
  print('Equilibrating GST weights')
  system, pdb = setup_system_implicit(filename_implicit)
  integrator = grestIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds, 2, 1)
  simulation = setup_simulation(system, pdb, integrator)

  ###First, equilibrate the weights:
  #The simulated tempering object:
  st = SimulatedSoluteTempering(simulation,
                              forceGroup=2,
                              cutoff=1e-8,
                              numTemperatures=numTemps,
                              tempChangeInterval=exchangeInterval,
                              minTemperature=minTemp*kelvin,
                              maxTemperature=maxTemp*kelvin,
                              reportInterval=100,
                              reportFile='./villin_gst_temp_equilibration.dat',
                             )

  stepsDone=0
  while st._weightUpdateFactor>st.cutoff:
    simulation.step(250)
    stepsDone+=250
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
                              reportFile='./villin_gst_temp.dat',
                             )
  #new st objects assume they need to equilibrate first. We don't.
  #so set the weights to what we determined, and turn off updating.
  st._weights = list(weights_store)
  st._updateWeights = False

  #now add the normal reporters and run!
  simulation.reporters.append(DCDReporter(output_directory+'villin_dip_gst_traj.dcd', dcdstride))
  simulation.reporters.append(StateDataReporter(stdout, 50000, step=True,totalSteps=number_steps+stepsDone,remainingTime=True,
                                              potentialEnergy=True, density=True, speed=True))
  simulation.step(number_steps)
  







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

  
