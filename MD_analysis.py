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
#from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
#run with python /Volumes/Jason/MD_analysis.py

pdb_inactive="AurB_t232d_3_007.pdb"
#trajectory_t232d.dcd
#trajectory_t232d
#t232d_3_007_1eqst.pdb
#sim2/t232d_run2_1eq.pdb
#/Volumes/Jason/Simulations/sim2/trajectory_t232d_3007_run2.dcd
pdb_wt="AukB_wt_3264_r46_2_008.pdb"
#trajectory_wt_080825.dcd
#wt_3264_r46_2_008_1eqst
#sim2/aukb_wt_3264_run2.pdb
def kt():
	pdb="AukB_k215r_t232d_64_128_19_11.pdb"
	traj="trajectory_AukB_k215r_t232d_64_128_19_11_resume_1.dcd"
	currname="trajectory_AukB_k215r_t232d_64_128_19_11_rem01_active"
	#reference="k215r_t232d_1eqst.pdb"
def k():
	#/Volumes/Jason/Simulations/k215r_1eqst.pdb
	pdb="k215r_64128_r88.pdb"
	traj="trajectory_k215r_64128_r88.dcd"
	currname="trajectory_k215r_64128_r88_rem0_active"
def t94():
	#t232a_94_1steq.pdb
	pdb="t232a_94.pdb"
	traj="trajectory_t232a_94.dcd" 
	#/Volumes/Jason/Simulations/trajectory_t232a_94_resume_0.dcd
def t73():
	#t232a_73_1steqst.pdb
	pdb="t232a_73.pdb"
	traj="trajectory_t232a_73.dcd" 

reference="t232d_3_007_1eqst.pdb"
pdb="t232a_94.pdb"
traj="trajectory_t232a_94_resume_0.dcd"
currname="trajectory_t232a_94_v_inactive_c0"
#trajectory_t232a_73_resume_0_c1
#trajectory_t232a_94_resume_0_c1

#traj1="trajectory_t232a_73_resume_0_c1.dcd"
#combtraj="trajectory_AukB_k215r_t232d_64_128_19_11_c0.dcd"
#comb_name="trajectory_t232a_73_resume_0_c0_c1.dcd"
timestep_ns = 0.000002
report_interval = 25000

SSD_PATH = "/Volumes/Jason/Simulations/"
def ssd_file(filename):
    return os.path.join(SSD_PATH, filename)

def __combiner__(pdb_file, traj0_file, traj1_file, comb_name):
	def dcdcombiner_mod():
		u = mda.Universe(
			ssd_file(pdb_file),
			[ssd_file(traj0_file), ssd_file(traj1_file)]
		)

		# Get timestep from first trajectory
		dt = u.trajectory.dt   # ps
		print(f"Detected dt = {dt} ps")

		n_atoms = u.atoms.n_atoms
		with DCDWriter(ssd_file(comb_name), n_atoms=n_atoms, dt=dt) as W:
			for ts in u.trajectory:
				W.write(u.atoms)
		print("Combined trajectory written successfully.")
	def comb_verify():
		u2 = mda.Universe(ssd_file(pdb_file), ssd_file(comb_name))
		times = [ts.time for ts in u2.trajectory]
		print("First 5 times:", times[:5])
		print("Last 5 times:", times[-5:])
		deltas = [round(times[i+1] - times[i], 5)
				for i in range(len(times)-1)]
		print("Δt unique:", set(deltas))
	return dcdcombiner_mod, comb_verify
#combiner, verifier = __combiner__(pdb, traj, traj1, comb_name)
#combiner()     
#verifier()  
trajnow = md.load(ssd_file(traj), top=ssd_file(pdb))
refer = md.load(ssd_file(reference))
def RMSD_visualizer_tot():
	n_frames = trajnow.n_frames

	steps = np.arange(n_frames) * report_interval
	times = steps * timestep_ns
	# Calculate RMSD (automatically aligns)
	atom_indices = trajnow.topology.select("backbone")
	traj_aligned = trajnow.superpose(refer, atom_indices=atom_indices)
	rmsd = md.rmsd(traj_aligned, refer, atom_indices=atom_indices,frame=0)  # 0 means align to first frame
	
	minima_idx = argrelextrema(rmsd, np.less)[0]
	print("Local minima indices:", minima_idx)
	print("Local minima RMSD values:", rmsd[minima_idx])

	plt.plot(times, rmsd, color="blue", label="RMSD")
	plt.xlabel("Time (NS)")
	plt.ylabel("RMSD")
	plt.title("RMSD fluctuations")
	plt.legend()
	plt.savefig(f"{currname}_rmsd.png")

	# Save results
	np.savetxt(f'{currname}.csv', rmsd)
	csv_file = ssd_file(f'{currname}_rmsd.csv')
	with open(csv_file, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(['Time_ns', 'RMSD_A'])  # header
		for t, r in zip(times, rmsd):
			writer.writerow([t, r])
	print("RMSD Done")
##New Memory safe code
def RMSD_visualizer_tot_safe(traj_file, top_file, currname, report_interval=1, timestep_ns=0.002):
    # Load trajectory
    trajnow = md.load(traj_file, top=top_file, stride=10)  # downsample every 10 frames
    refer = trajnow[0]  # first frame as reference

    # Only use backbone atoms
    atom_indices = trajnow.topology.select("backbone")

    # Compute RMSD frame by frame
    rmsd = md.rmsd(trajnow, refer, atom_indices=atom_indices)  # aligns internally

    # Time axis
    n_frames = trajnow.n_frames
    steps = np.arange(n_frames) * report_interval * 10  # multiply by stride
    times = steps * timestep_ns  # ns

    # Local minima
    minima_idx = argrelextrema(rmsd, np.less)[0]
    print("Local minima indices:", minima_idx)
    print("Local minima RMSD values:", rmsd[minima_idx])

    # Plot RMSD
    plt.figure(figsize=(7,4))
    plt.plot(times, rmsd, color="blue", label="RMSD")
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title(f"{currname} RMSD fluctuations")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{currname}_rmsd.png", dpi=300)
    plt.close()

    # Save RMSD to CSV
    csv_file = f"{currname}_rmsd.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Time_ns', 'RMSD_A'])
        for t, r in zip(times, rmsd):
            writer.writerow([t, r])

    print("RMSD Done")
def RMSF_visualizer_trial():
    print("RMSF")
    # Align to first frame instead of average structure
    align.AlignTraj(u, u, 
                    select='protein and resid 1:150 and backbone',
                    ref_frame=0).run()

    c_alphas = u.select_atoms('protein and resid 1:150 and name CA')
    u.trajectory.rewind()

    R = rms.RMSF(c_alphas).run()

    rmsf_values = R.results.rmsf
    residue_ids = c_alphas.resids

    plt.figure(figsize=(7,4))
    plt.plot(residue_ids, rmsf_values)
    plt.xlabel('Residue number')
    plt.ylabel('RMSF (Å)')

    plt.axvspan(142, 144, alpha=0.2, label='DFG')
    plt.axvspan(8, 13, alpha=0.2, label='G-loop')
    plt.axvspan(43, 53, alpha=0.2, label="αC helix")

    plt.legend()
    plt.tight_layout()

    outpng = ssd_file(f"{currname}_rmsf.png")
    plt.savefig(outpng, dpi=300)
    plt.close()

    df_to_save = pd.DataFrame({
        'Residue': residue_ids,
        'RMSF_Angstrom': rmsf_values
    })

    csv_path = ssd_file(f"{currname}_RMSF.csv")
    df_to_save.to_csv(csv_path, index=False)

    print("RMSF Done")
def RMSF_visualizer():
	print("RMSF")
	u = mda.Universe(ssd_file(pdb), ssd_file(traj))
	average = align.AverageStructure(u, u, select='protein and resid 1:250 and backbone', ref_frame=0).run()
	ref = average.results.universe
	align.AlignTraj(u, ref, select='protein and resid 1:250 and backbone', in_memory=True).run()
	c_alphas = u.select_atoms('protein and resid 1:250 and name CA')
	u.trajectory.rewind()
	#protein backbone and name CA
	R = rms.RMSF(c_alphas).run()
	rmsf_values = R.results.rmsf      # Å in MDAnalysis
	residue_ids = c_alphas.resids

	#normal: 142-144 for DFG; 8-13 for GL; 43-53 for aC
	plt.figure(figsize=(7,4))
	plt.plot(c_alphas.resids, R.results.rmsf)
	plt.xlabel('Residue number')
	plt.ylabel('RMSF (Å)')
	plt.axvspan(142, 144, alpha=0.2, label='DFG', color="lime")
	plt.axvspan(8, 13, alpha=0.2, label='G-loop', color="red")
	plt.axvspan(43, 53, alpha=0.2, label="αC helix", color="yellow")

	plt.legend()
	plt.tight_layout()

	outpng = ssd_file(f"{currname}_rmsf.png")
	plt.savefig(outpng, dpi=300)
	plt.close()
	csv_path = ssd_file(f"{currname}_RMSF.csv")

	df_to_save = pd.DataFrame({
		'Residue': residue_ids,
		'RMSF_Angstrom': rmsf_values
	})

	df_to_save.to_csv(csv_path, index=False)
	print("RMSF Done")
def RMSD_visualizer(start_residue, end_residue):
	atom_indices = trajnow.topology.select(f"resid {start_residue} to {end_residue}")
	n_frames = trajnow.n_frames
	steps = np.arange(n_frames) * report_interval
	times = steps * timestep_ns

	traj_aligned = trajnow.superpose(refer, atom_indices=atom_indices)
	rmsd = md.rmsd(traj_aligned, refer, atom_indices=atom_indices,frame=0)  # 0 means align to first frame
	
	minima_idx = argrelextrema(rmsd, np.less)[0]
	print("Local minima indices:", minima_idx)
	print("Local minima RMSD values:", rmsd[minima_idx])

	plt.plot(times, rmsd, color="blue", label=f"RMSD (Residues {start_residue}-{end_residue})")
	plt.xlabel("Time (ps)")
	plt.ylabel("RMSD (Å)")
	plt.legend()
	plt.savefig(f"{currname}_{start_residue}_{end_residue}_rmsd.png")
	
	#np.savetxt(f"{currname}_{start_residue}_{end_residue}_rmsd.dat", rmsd)
	csv_file = ssd_file(f"{currname}_{start_residue}_{end_residue}_rmsd.csv")
	with open(csv_file, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(['Time_ns', 'RMSD_A'])  # header
		for t, r in zip(times, rmsd):
			writer.writerow([t, r])
	print(f"RMSD for residues {start_residue}-{end_residue} done.")
def distance_visualizer(res1, res2, selection):
	"""
	Compute and visualize distance between two residues (default: CA atoms)
	across trajectory using MDAnalysis.
	"""
	##usually atom_name is CA
	#salt bridge
	u = mda.Universe(ssd_file(pdb), ssd_file(traj))
	if selection==0:
		sel1 = u.select_atoms(f"resid {res1} and name CA")
		sel2 = u.select_atoms(f"resid {res2} and name CA")
		#nz and OE1/2
	#cation-pi proxy ##double check
	if selection==1:	
		sel1 = u.select_atoms(f"resid {res1} and name NZ")
		sel2 = u.select_atoms(f"resid {res2} and (name OE1 OE2)")
		sel2 = u.select_atoms(f"resid {res2} and (name CG CZ)")

	if len(sel1) == 0 or len(sel2) == 0:
		raise ValueError(f"Could not find atom atoms in residues {res1} or {res2}.")

	# Compute distances for each frame
	distances_list = []
	times = []
	steps= []
	for ts in u.trajectory:
		#d = distances.distance_array(sel1.positions, sel2.positions)[0, 0]

		#center of mass
		d = np.linalg.norm(sel1.center_of_mass() - sel2.center_of_mass())

		frame = ts.frame
		step = frame * report_interval
		time_ps = step * timestep_ns
	
		distances_list.append(d)
		times.append(time_ps)  # picoseconds if time info is present
		steps.append(step)

	distances_list = np.array(distances_list)
	times = np.array(times)
	steps=np.array(steps)

	# Identify local minima (optional)
	minima_idx = argrelextrema(distances_list, np.less)[0]
	print("Local minima indices:", minima_idx)
	print("Local minima distances (Å):", distances_list[minima_idx])

	# Plot
	plt.plot(times, distances_list, color="green", label=f"Distance: resid {res1}-{res2} (atoms)")
	plt.xlabel("Time (ps)")
	plt.ylabel("Distance (Å)")
	plt.legend()
	plt.tight_layout()
	plt.savefig(f"{currname}_{res1}_{res2}.png")

	# Save data
	#np.savetxt(f"{currname}_{res1}_{res2}_dist.dat", distances_list)
	csv_file = ssd_file(f"{currname}_{res1}_res{res2}_dist.csv")
	print("About to write CSV:", csv_file)
	with open(csv_file, 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(['Time_ns', 'Distance_A'])
		for t, d in zip(times, distances_list):
			writer.writerow([t, d])
	print(f"Distance analysis for residues {res1}-{res2} done.")
def grapher():
	df = pd.read_csv(f'/Volumes/Jason/Simulations/{currname}_rmsd.csv')
	time_col = df['Time_ns']
	# Plot
	plt.figure(figsize=(10,5))
	plt.plot(time_col, df['RMSD_A'])
	plt.xlabel('Time (ns)')
	plt.ylabel('RMSD (Å)')
	plt.title('RMSD vs Time')
	plt.show()

#test the one below
def distance_visualizer_test_test(res1, res2, mode="min"):
    """
    mode:
        "ca"  → CA–CA distance
        "min" → minimum heavy atom distance (recommended)
    """

    u = mda.Universe(ssd_file(pdb), ssd_file(traj))

    if mode == "ca":
        sel1 = u.select_atoms(f"resid {res1} and name CA")
        sel2 = u.select_atoms(f"resid {res2} and name CA")
    else:
        sel1 = u.select_atoms(f"resid {res1} and not name H*")
        sel2 = u.select_atoms(f"resid {res2} and not name H*")

    if len(sel1) == 0 or len(sel2) == 0:
        raise ValueError(f"Selection failed for residues {res1}, {res2}")

    distances_list = []
    times = []

    for ts in u.trajectory:
        box = ts.dimensions  # PBC box

        d_array = distances.distance_array(
            sel1.positions,
            sel2.positions,
            box=box
        )

        if mode == "ca":
            d = d_array[0, 0]
        else:
            d = d_array.min()

        time_ns = ts.frame * report_interval * timestep_ns

        distances_list.append(d)
        times.append(time_ns)

    distances_list = np.array(distances_list)
    times = np.array(times)

    # Plot
    plt.figure(figsize=(7,4))
    plt.plot(times, distances_list, label=f"{res1}-{res2} ({mode})")
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance (Å)")
    plt.title(f"Residue {res1}-{res2} Distance")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{currname}_{res1}_{res2}_{mode}.png", dpi=300)
    plt.close()

    print(f"Distance ({mode}) analysis done for {res1}-{res2}")
RMSD_visualizer(142, 144)
distance_visualizer(30, 49, 0)
RMSD_visualizer(160, 180)
RMSD_visualizer_tot()
#RMSF_visualizer()

print(reference)

#from sklearn.decomposition import PCA

def conformational_landscape():
    u = md.load(ssd_file(traj), top=ssd_file(pdb))

    atom_indices = u.topology.select("backbone")
    u.superpose(u, 0, atom_indices=atom_indices)

    X = u.xyz[:, atom_indices, :].reshape(u.n_frames, -1)

    pca = PCA(n_components=2)
    proj = pca.fit_transform(X)

    plt.hexbin(proj[:,0], proj[:,1], gridsize=50, bins='log')
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Free Energy Landscape")
    plt.colorbar(label="log density")
    plt.savefig(f"{currname}_FEL.png")

    print("Conformational landscape done")
def pocket_rgyr(res_start=63, res_end=144):
    u = mda.Universe(ssd_file(pdb), ssd_file(traj))

    pocket = u.select_atoms(f"resid {res_start}:{res_end}")

    rg = []
    times = []

    for ts in u.trajectory:
        rg.append(pocket.radius_of_gyration())
        times.append(ts.frame * report_interval * timestep_ns)

    plt.plot(times, rg)
    plt.xlabel("Time (ns)")
    plt.ylabel("Rg (Å)")
    plt.title("Binding Pocket Compactness")
    plt.savefig(f"{currname}_pocket_rg.png")

    print("Pocket Rg done")

#conformational_landscape()
#pocket_rgyr()
#grapher()




###fix
def radius_gyration():
	u=mda.Universe("AukB_k215r_t232d_64_128_19_11.pdb", "trajectory_AukB_k215r_t232d_64_128_19_11.dcd")
	rgyr = []
	time = []
	protein = u.select_atoms("protein")
	for ts in u.trajectory:
	    time.append(u.trajectory.time)
	    rgyr.append(protein.radius_of_gyration())
	rgyr_df = pd.DataFrame(rgyr, columns=['Radius of gyration (A)'], index=time)
	rgyr_df.index.name = 'Time (ps)'
	rgyr_df.head()
	rgyr_df.plot(title='Radius of gyration')
	plt.savefig("radius of Gyration k215q")
def HBondanal():
	u=mda.Universe("AukB_k215q_rank46_3264.pdb", "trajectory_k215q.dcd")
	# Set up hydrogen bond analysis
	hbonds = HBA(universe=u,
		    donors_sel="name O or name N or name OH",  # Common donors
		    hydrogens_sel="name H*",                   # All hydrogens
		    acceptors_sel="name O or name N",          # Common acceptors
		    d_h_cutoff=3.0,                           # Donor-Hydrogen distance cutoff (Å)
		    d_a_cutoff=3.0,                           # Donor-Acceptor distance cutoff (Å)
		    d_h_a_angle_cutoff=150,                   # D-H-A angle cutoff (degrees)
		    update_selections=True)                    # Update selections for each frame

	# Run analysis
	hbonds.run()

	# Generate results
	results = hbonds.results.hbonds

	# Plot time evolution of number of H-bonds
	plt.plot(hbonds.results.times, hbonds.results.count_by_time())
	plt.xlabel("Time (ps)")
	plt.ylabel("Number of H-bonds")
	plt.title("Hydrogen Bond Dynamics")
	plt.savefig("HBond k215q")
	
	np.savetxt("hbonds_all.dat", hbonds.results.hbonds, 
           fmt="%8d %8d %8d %8d %8.3f %8.3f %8.3f",
           header="Frame DonorID HydrogenID AcceptorID Distance Angle")

	# Get unique hydrogen bonds (persistent throughout simulation)
	unique_hbonds = hbonds.count_by_ids()
	print(f"Found {len(unique_hbonds)} unique H-bond pairs")

	# Plot hydrogen bond lifetime
	lifetimes = hbonds.lifetime()
	plt.bar(range(len(lifetimes)), lifetimes)
	plt.xlabel("H-bond ID")
	plt.ylabel("Lifetime (frames)")
	plt.show()
	plt.savefig("unique_bond_k215q")
	
	network = hbonds.network()

	# Get most common H-bonds
	from collections import Counter
	common_bonds = Counter((bond[1], bond[3]) for bond in hbonds.results.hbonds)  # (donor, acceptor)
	print("Most common H-bonds:", common_bonds.most_common(5))

	# Plot donor-acceptor distance for specific pair
	specific_pair = hbonds.results.hbonds[
	    (hbonds.results.hbonds[:,1] == donor_id) & 
	    (hbonds.results.hbonds[:,3] == acceptor_id)]  # Use actual atom IDs

	plt.plot(specific_pair[:,0], specific_pair[:,4])  # Frame vs Distance
	plt.xlabel("Frame")
	plt.ylabel("Donor-Acceptor Distance (Å)")
	plt.savefig("hbond network k215q")