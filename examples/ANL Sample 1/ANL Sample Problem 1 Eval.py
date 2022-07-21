#%%
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 11:47:15 2022

@author: mjmullin
"""

from srlife import receiver, solverparams, library, thermal, structural, system, damage, managers

# Load the receiver we previously saved
directory = ''
modelName = 'sample1_noCloud'
model = receiver.Receiver.load((directory+modelName+".hdf5"))


#Cut down on run time for now by making the tube analyses 2D
for panel in model.panels.values():
    for tube in panel.tubes.values():
        tube.make_2D(tube.dim[2])

# Choose the material models
fluid_mat = library.load_fluid("salt", "base") # Generic chloride salt  model
material = '740H'
deformation_mat_string = 'base' #elastic-model
thermal_mat, deformation_mat, damage_mat = library.load_material(material, "base", deformation_mat_string, "base")

# Setup some solver parameters
params = solverparams.ParameterSet()
params['progress_bars'] =True # Print a progress bar to the screen as we solve
params['nthreads'] = 8 # Solve will run in multithreaded mode, set to number of available cores
params['system']['atol'] = 1.0e-4 # During the standby very little happens, lower the atol to accept this result
# params['structural']['atol'] = 1.0e-3
# params['system']['atol'] = 1.0e-3
# params['system']['miter'] = 20
# params['structural']['miter'] = 50
# params['structural']['rtol'] = 1.0e-5

# Choose the solvers, i.e. how we are going to solve the thermal,
# single tube, structural system, and damage calculation problems.
# Right now there is only one option for each
# Define the thermal solver to use in solving the heat transfer problem
thermal_solver = thermal.FiniteDifferenceImplicitThermalSolver(
    params["thermal"])
# Define the structural solver to use in solving the individual tube problems
structural_solver = structural.PythonTubeSolver(params["structural"])
# Define the system solver to use in solving the coupled structural system
system_solver = system.SpringSystemSolver(params["system"])
# Damage model to use in calculating life
damage_model = damage.TimeFractionInteractionDamage(params["damage"])



# The solution manager
solver = managers.SolutionManager(model, thermal_solver, thermal_mat, fluid_mat, 
                                  structural_solver, deformation_mat, damage_mat, system_solver, damage_model, pset = params)

# Reset the temperatures each night
solver.add_heuristic(managers.CycleResetHeuristic())

life = solver.solve_life()
print("Best estimate life: %f daily cycles" % life)

#%%
# # Save the tube data out for additional visualization
# for pi, panel in model.panels.items():
#   for ti, tube in panel.tubes.items():
#     tube.write_vtk("tube-%s-%s" % (pi, ti))

# %%
