import MDAnalysis as mda

#u = mda.Universe("/Volumes/Jason/Simulations/AukB_k215r_t232d_64_128_19_11.pdb", "/Volumes/Jason/Simulations/trajectory_4af3.dcd")

u=mda.Universe("/Volumes/Jason/Simulations/AukB_wt_3264_r46_2_008.pdb", "/Volumes/Jason/Simulations/trajectory_wt_080825.dcd")
n_frames = len(u.trajectory)
dt = u.trajectory.dt   # in picoseconds
total_time_ps = (n_frames - 1) * dt
total_time_ns = total_time_ps / 1000

print("Frames:", n_frames)
print("dt (ps):", dt)
print("Total time (ns):", total_time_ns)
