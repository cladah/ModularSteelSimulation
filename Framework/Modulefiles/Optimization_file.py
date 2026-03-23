from .ModuleStructure_file_new import CalcModule
import numpy as np
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from Framework.HelpFile import read_input
from tkinter import filedialog
from pathlib import Path
import shutil
import tkinter as tk
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import json

class OptimizationModule(CalcModule):

    def __init__(self, infile, modulenr, simulation_method):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Optimization", infile, modulenr)

        initial_dir = Path.cwd() / "Experimental"
        if not initial_dir.exists():
            initial_dir = Path.cwd()
        root = tk.Tk()
        root.withdraw()
        selected_dir = filedialog.askopenfilename(
            initialdir=str(initial_dir),
            title="Select experimental data for optimization"
        )
        root.destroy()
        df = pd.read_csv(selected_dir)
        self.goal_data = df.to_numpy()
        self.goal_file = str(selected_dir)
        self.goal_par = self.minput["OptPar"]
        self.par_type = list()
        self.par_files = list()
        self.initial_guess = list()
        self.opt_path = Path.cwd() / "Inputs" / "TmpOpt"
        self.run_simulation = simulation_method

    def run(self):
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Optimization module: " + self.inputfile,
                  "File used for optimization: " + self.goal_file,
                  "Model: " + str(self.minput["OptType"]) + "\n",
                  "Parameters used for optimization:"]
        for par in self.minput["Parameters"]:
            tmpfile = self.minput["Parameters"][par][0]
            self.par_files.append(tmpfile)
            self.par_type.append(par)
            self.initial_guess.append(self.minput["Parameters"][par][1])
            outstr.append(f"{par} in {tmpfile}")

        outstr.append(f"\n{self.goal_par} used as objective function")
        outstr.append("\n---------------------------------------------------------------------\n")

        for line in outstr:
            self.writeoutput(line)
            print(line)

        self.new_input_directory()

        initial_guess = np.array(self.initial_guess)

        result = minimize(
            self.objective_function,
            initial_guess,
            method=self.minput["OptType"],
            options={'xatol': 1e-2,
                     'fatol': 1e-5,
                     'disp': True}
        )

        if result.success:
            print(f"Convergence reached. Optimized parameters: {result.x}")
        else:
            print("Optimization failed.")


    def objective_function(self, params):

        x_exp = self.goal_data[:,0]
        y_exp = self.goal_data[:, 1]

        self.change_input(params)
        self.run_simulation(True, self.opt_path)

        # Extracting data from Datastream file
        x_sim = getaxisvalues("nodes")
        y_sim = getaxisvalues(self.goal_par)

        # Interpolation due to missmatch between simulation and experiment
        interp_func = interp1d(x_sim, y_sim, kind='linear', fill_value="extrapolate")
        y_sim_interp = interp_func(x_exp)

        # Mean square error
        mse = np.mean((y_sim_interp - y_exp) ** 2)
        return mse

    def change_input(self, params):
        for i in range(len(params)):
            filename = self.opt_path / self.par_files[i]
            if filename:
                with open(filename, 'r+') as f:
                    data = json.load(f)
                    data[self.par_type[i]] = params[i]
                    f.seek(0)
                    json.dump(data, f, indent=4)
                    f.truncate()

    def new_input_directory(self):
        base_path = Path.cwd() / "Inputs"
        source_dir = base_path / self.ginput["InputDirectory"]
        target_dir = self.opt_path

        try:
            if not source_dir.exists():
                print(f"Error: Source directory {source_dir} does not exist.")
                return

            # If the target directory already exists, remove it to ensure a clean copy
            if target_dir.exists():
                shutil.rmtree(target_dir)

            # Copy the entire directory tree
            shutil.copytree(source_dir, target_dir)

            print(f"Successfully copied:\n{source_dir}\nto\n{target_dir}")

        except Exception as e:
            print(f"An error occurred: {e}")

        return target_dir


def optimize(self, exp_data, e_range=(1e5, 1e7), nu_range=(0.2, 0.4)):
    """Monte Carlo Wrapper."""
    best_error = float('inf')
    best_params = None

    print(f"Starting Monte Carlo optimization ({self.n_samples} iterations)...")

    for i in range(self.n_samples):
        # Random Sampling
        test_E = np.random.uniform(*e_range)
        test_nu = np.random.uniform(*nu_range)

        try:
            sim_result = self.run_fenics_sim(test_E, test_nu)
            error = self.objective_function(sim_result, exp_data)

            if error < best_error:
                best_error = error
                best_params = (test_E, test_nu)
                print(f"Iter {i}: New Best! E={test_E:.2e}, nu={test_nu:.3f}, Error={error:.2e}")
        except Exception as e:
            continue

    return best_params, best_error


def run_parallel_optimization(n_samples, exp_stresses, sensor_pts):
    # Prepare random parameters
    E_samples = np.random.uniform(1e5, 1e7, n_samples)
    nu_samples = np.random.uniform(0.2, 0.4, n_samples)

    # Bundle arguments for the executor
    tasks = [(E_samples[i], nu_samples[i], sensor_pts, exp_stresses) for i in range(n_samples)]

    # Use all available cores minus one to keep the OS responsive
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    print(f"Launching {n_samples} simulations across {num_workers} cores...")

    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # map distributes the tasks and yields results as they finish
        for result in executor.map(run_single_simulation, tasks):
            results.append(result)
            if len(results) % 10 == 0:
                print(f"Progress: {len(results)}/{n_samples} completed.")

    # Find the best result
    best_result = min(results, key=lambda x: x[2])
    return best_result