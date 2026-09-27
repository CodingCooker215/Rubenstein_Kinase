import pandas as pd
import matplotlib.pyplot as plt

# Load OpenMM log file
file1="/Volumes/Jason/old_data/log_equilibrium_k215r_64128_r88.pdb_resume_0.txt"
file2="/Volumes/Jason/Simulations/00_log_equilibrium_k215r_64128_r88.txt"
merge="/Volumes/Jason/Simulations/00_log_equilibrium_k215r_64128_r88_r2.txt"
range="k215r"
def merge_logs(file1, file2, merged_file):
    with open(file1) as f1, open(file2) as f2, open(merged_file, "w") as out:
        out.write(f1.read())

        next(f2) 
        out.write("\n") 
        out.write(f2.read())
def plot(file):
    df = pd.read_csv(file, sep="\t")
    print(df.columns)

    # Extract data
    steps = df['#"Step"']
    potential_energy = df["Potential Energy (kJ/mole)"]
    dt_ps = 0.002
    time_ns = steps * dt_ps * 0.001

    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(time_ns, potential_energy)
    plt.xlabel("Time (ns)")
    plt.ylabel("Potential Energy (kJ/mol)")
    plt.title(f"Potential Energy vs Simulation Step for {range}")
    plt.tight_layout()
    plt.savefig(f"/Volumes/Jason/Simulations/equil_{range}.png", dpi=300)
    plt.show()

#merge_logs(file1, file2, merge)
plot(merge)