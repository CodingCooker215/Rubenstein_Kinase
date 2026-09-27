import MDAnalysis as mda
import os
from MDAnalysis.analysis import align, rms, hbonds, rms, align
from MDAnalysis.analysis.bat import BAT
import pandas as pd
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from MDAnalysis.coordinates.DCD import DCDWriter
import csv
from MDAnalysis.analysis import distances
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# ============================
# GLOBAL SETTINGS
# ============================

timestep_ns = 0.000002
report_interval = 25000

SSD_PATH = "/Volumes/Jason/Simulations/"
def ssd_file(filename):
    return os.path.join(SSD_PATH, filename)

# ============================
# GLOBAL VARIABLES (set per job)
# ============================

pdb = None
traj = None
reference = None
currname = None
trajnow = None
refer = None

# ============================
# ORIGINAL FUNCTIONS (UNCHANGED)
# ============================
def plot(t, values, name):
    plt.plot(t, values)
    plt.xlabel("Time (NS)")
    plt.ylabel("RMSD")
    plt.savefig(f"{name}_rmsd.png")
    plt.close()

def RMSD_visualizer_tot():
    n_frames = trajnow.n_frames
    steps = np.arange(n_frames) * report_interval
    times = steps * timestep_ns

    atom_indices = trajnow.topology.select("backbone")
    traj_aligned = trajnow.superpose(refer, atom_indices=atom_indices)
    rmsd_vals = md.rmsd(traj_aligned, refer, atom_indices=atom_indices, frame=0)

    minima_idx = argrelextrema(rmsd_vals, np.less)[0]
    print("Local minima indices:", minima_idx)

    # plt.plot(times, rmsd_vals)
    # plt.xlabel("Time (NS)")
    # plt.ylabel("RMSD")
    # plt.savefig(f"{currname}_rmsd.png")
    # plt.close()
    
    #plot(times, rmsd_vals, currname)

    csv_file = ssd_file(f"{currname}_rmsd.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Time_ns', 'RMSD_A'])
        for t, r in zip(times, rmsd_vals):
            writer.writerow([t, r])

def RMSD_visualizer(start_residue, end_residue):
    atom_indices = trajnow.topology.select(f"resid {start_residue} to {end_residue}")

    n_frames = trajnow.n_frames
    steps = np.arange(n_frames) * report_interval
    times = steps * timestep_ns

    traj_aligned = trajnow.superpose(refer, atom_indices=atom_indices)
    rmsd_vals = md.rmsd(traj_aligned, refer, atom_indices=atom_indices, frame=0)

    # plt.plot(times, rmsd_vals)
    # plt.xlabel("Time")
    # plt.ylabel("RMSD")
    # plt.savefig(f"{currname}_{start_residue}_{end_residue}.png")
    # plt.close()

    # name=f"{currname}_{start_residue}_{end_residue}"
    # plot(times, rmsd_vals,name)

    csv_file = ssd_file(f"{currname}_{start_residue}_{end_residue}_rmsd.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Time_ns', 'RMSD_A'])
        writer.writerows(zip(times, rmsd_vals))

def distance_visualizer(res1, res2, selection):
    u = mda.Universe(ssd_file(pdb), ssd_file(traj))

    sel1 = u.select_atoms(f"resid {res1} and name CA")
    sel2 = u.select_atoms(f"resid {res2} and name CA")

    distances_list = []
    times = []

    for ts in u.trajectory:
        d = np.linalg.norm(sel1.center_of_mass() - sel2.center_of_mass())
        time_ns = ts.frame * report_interval * timestep_ns
        distances_list.append(d)
        times.append(time_ns)

    # plt.plot(times, distances_list)
    # plt.xlabel("Time")
    # plt.ylabel("Distance")
    # plt.savefig(f"{currname}_{res1}_{res2}.png")
    # plt.close()

    # names=f"{currname}_{res1}_{res2}"
    # plot(times, distances_list, names)

    csv_file = ssd_file(f"{currname}_{res1}_{res2}_distance.csv")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Time_ns', 'Distance_A'])
        writer.writerows(zip(times, distances_list))

# ============================
# PARALLEL JOB RUNNER
# ============================

def run_job(job):
    global pdb, traj, reference, currname, trajnow, refer

    pdb = job["pdb"]
    traj = job["traj"]
    reference = job["reference"]
    currname = job["name"]

    print(f"Running: {currname}")

    trajnow = md.load(ssd_file(traj), top=ssd_file(pdb))
    refer = md.load(ssd_file(reference))

    # === YOUR PIPELINE ===
    RMSD_visualizer(142, 144)
    distance_visualizer(30, 49, 0)
    RMSD_visualizer(160, 180)
    #only run rmsd_visualizer_tot when reference is itself
    #RMSD_visualizer_tot()   

    return {
        "name": currname,
        "traj": traj,
        "reference": reference
    }

# ============================
# JOB DEFINITIONS
# ============================

jobs = [
    {
        "name": "t232a_73_v_r1inactive_run2",
        "pdb": "t232a_73.pdb",
        "traj": "sim2/trajectory_t232a_73_run2.dcd",
        "reference": "t232d_3_007_1eqst.pdb"
    }, 
    {
        "name": "t232a_94_v_r1self_run2",
        "pdb": "t232a_94.pdb",
        "traj": "sim2/trajectory_t232a_94_run2.dcd",
        "reference": "t232a_94_1steq.pdb"
    }, 
    {
        "name": "t232a_94_v_r1active_run2",
        "pdb": "t232a_94.pdb",
        "traj": "sim2/trajectory_t232a_94_run2.dcd",
        "reference": "wt_3264_r46_2_008_1eqst.pdb"
    }, 
    {
        "name": "t232a_94_v_r1inactive_run2",
        "pdb": "t232a_94.pdb",
        "traj": "sim2/trajectory_t232a_94_run2.dcd",
        "reference": "t232d_3_007_1eqst.pdb"
    }
]
# jobs = [
#     {
#         "name": "t232a_94_self_run2",
#         "pdb": "t232a_94.pdb",
#         "traj": "sim2/trajectory_t232a_94_run2.dcd",
#         "reference": "sim2/t232a_94_1steq_r2.pdb"
#     }
# ]
# ============================
# PARALLEL EXECUTION
# ============================

# def run_all_jobs(jobs):
#     results = []
#     n_workers = 1

#     with ProcessPoolExecutor(max_workers=n_workers) as executor:
#         futures = [executor.submit(run_job, job) for job in jobs]

#         for future in as_completed(futures):
#             res = future.result()
#             results.append(res)

#     return results

def run_all_jobs(jobs):
    results = []
    for job in jobs:
        res = run_job(job)
        results.append(res)
    return results

# ============================
# REPORT GENERATION
# ============================
import datetime
def write_report(results):
    report_path = ssd_file("parallel_run_report.json")

    # load existing report
    if os.path.exists(report_path):
        with open(report_path, "r") as f:
            try:
                existing = json.load(f)
            except json.JSONDecodeError:
                existing = []
    else:
        existing = []

    # timestamp
    now = datetime.datetime.now().isoformat()

    # add timestamp to each result
    for r in results:
        r["timestamp"] = now
        existing.append(r)

    # write back
    with open(report_path, "w") as f:
        json.dump(existing, f, indent=4)

    print(f"Report appended to {report_path}")

# ============================
# MAIN
# ============================

if __name__ == "__main__":
    results = run_all_jobs(jobs)
    write_report(results)