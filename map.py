import pandas as pd
import matplotlib.pyplot as plt

ref="active_r2"
range="total"

# kt="/Users/jasonhwang/Documents/Brown/Rubenstein/AukB_k215r_t232d_64_128_19_11_rmsf.csv"
# k="/Users/jasonhwang/Documents/Brown/Rubenstein/AukB_k215q_rank46_3264_rmsf.csv"
# inact="/Users/jasonhwang/Documents/Brown/Rubenstein/AurB_t232d.csv"
# wt="/Users/jasonhwang/Documents/Brown/Rubenstein/AurB_4af3.csv"
# ta_94="/Volumes/Jason/Simulations/trajectory_t232a_94_v_active_r0_t3_RMSF.csv"
# ta_73="/Volumes/Jason/Simulations/trajectory_t232a_73v_active_r0_t3_RMSF.csv"

kt="/Volumes/Jason/Simulations/4March2026/total/done/00_trajectory_k215r_t232d_total_self_dist_total.csv"
k="/Volumes/Jason/Simulations/run2/trajectory_k215r_r88_run2_v_active_rmsd.csv"
inact="/Volumes/Jason/Simulations/run2/trajectory_t232d_3007_run2_v_active_rmsd.csv"
wt="/Volumes/Jason/Simulations/run2/trajectory_wt_3264_r46_run2_v_active_rmsd.csv"
ta_94="/Volumes/Jason/Simulations/run2/trajectory_t232a_94_run2_v_active_rmsd.csv"
ta_73="/Volumes/Jason/Simulations/run2/trajectory_t232a_73_run2_v_active_rmsd.csv"
#file="/Volumes/Jason/Simulations/4March2026/cook"

#run1="/Volumes/Jason/Simulations/4March2026/total/done/00_trajectory_k215r_total_inact_dist_total.csv"
#run2="/Volumes/Jason/Simulations/run2/trajectory_k215r_r88_run2_v_inactive_rmsd.csv"
file="/Volumes/Jason/Simulations/sim2/data"
#f1 = pd.read_csv(f"{run1}")
#f2 = pd.read_csv(f"{run2}")
f1 = pd.read_csv(f"{inact}")
f2 = pd.read_csv(f"{wt}")
f3 = pd.read_csv(f"{kt}")
f4 = pd.read_csv(f"{k}")
f5 = pd.read_csv(f"{ta_94}")
f6 = pd.read_csv(f"{ta_73}")


#f1 = pd.read_csv("/Volumes/Jason/Simulations/RMSF/AurB_t232d_3_007_reeun_remsf_RMSF.csv")
#f2 = pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_4af3_inactive_RMSF.csv")
# f3 = pd.read_csv("/Volumes/Jason/Simulations/RMSF/AukB_k215r_t232d_64_128_19_11_c0_RMSF.csv")
# f4 = pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_k215r_64128_r88_rem02_RMSF.csv")

plt.figure(figsize=(10,6))
def rmsd():
    time_col = f3.columns[0]
    value_col = f3.columns[1]

    print(f1.columns)
    print(f3.columns)

    # Zero each time axis
    t1 = (f1[time_col] - f1[time_col].iloc[0])
    t2 = (f2[time_col] - f2[time_col].iloc[0])
    #t3 = (f3[time_col] - f3[time_col].iloc[0])
    t4 = (f4[time_col] - f4[time_col].iloc[0])
    t5 = (f5[time_col] - f5[time_col].iloc[0])
    t6 = (f6[time_col] - f6[time_col].iloc[0])
    plt.plot(t1, f1[value_col], label='AukB_T232D')
    plt.plot(t2, f2[value_col], label='AukB_wt')
    #plt.plot(t1, f1[value_col], label='run1')
    #plt.plot(t2, f2[value_col], label='run2')
    #plt.plot(t3, f3[value_col], label='AukB_k215r_t232d')
    plt.plot(t4, f4[value_col], label='k215r_64128_r88')
    plt.plot(t5, f5[value_col], label='t232a_94')
    plt.plot(t6, f6[value_col], label='t232a_73')

    plt.xlabel("Time (ns, aligned)")
    plt.ylabel(f"{range} RMSD (Å)")
    plt.title("Fluctuation Comparison (Time-Aligned)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{file}/{range}_rmsd_{ref}.png", dpi=300)
    #_wo_t232a
def rmsd_wo_t232a():
    time_col = f3.columns[0]
    value_col = f3.columns[1]

    print(f1.columns)
    print(f3.columns)

    # Zero each time axis
    t1 = (f1[time_col] - f1[time_col].iloc[0])
    t2 = (f2[time_col] - f2[time_col].iloc[0])
    t3 = (f3[time_col] - f3[time_col].iloc[0])
    t4 = (f4[time_col] - f4[time_col].iloc[0])
    plt.plot(t1, f1[value_col], label='AukB_T232D')
    plt.plot(t2, f2[value_col], label='AukB_wt')
    plt.plot(t3, f3[value_col], label='AukB_k215r_t232d')
    plt.plot(t4, f4[value_col], label='k215r_64128_r88')

    plt.xlabel("Time (ns, aligned)")
    plt.ylabel(f"{range} RMSD (Å)")
    plt.title("Fluctuation Comparison (Time-Aligned)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{file}/{range}_rmsd_{ref}_wo_t232a.png", dpi=300)
def rmsf_test_Test():
    plt.plot(t1, f1[value_col], label='inactive')
    plt.plot(t2, f2[value_col], label='active')
    plt.plot(t3, f3[value_col], label='k215r_t232d')
    plt.plot(t4, f4[value_col], label='k215r')
    plt.plot(t5, f5[value_col], label='t232a_94')
    plt.plot(t6, f6[value_col], label='t232a_73')

    plt.axvspan(142, 144, alpha=0.2, label='DFG', color="lime")
    plt.axvspan(8, 13, alpha=0.2, label='G-loop', color="red")
    plt.axvspan(43, 53, alpha=0.2, label="αC helix", color="yellow")

    plt.xlabel('Residue number')
    plt.ylabel('RMSF (Å)')
    plt.title("Aurora Kinase RMSF Comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{range}_rmsf_{ref}.png", dpi=300)
def standardize(df):
    df.columns = ["Residue", "RMSF"]
    return df
f1 = standardize(f1)
f2 = standardize(f2)
f3 = standardize(f3)
f4 = standardize(f4)
f5 = standardize(f5)
f6 = standardize(f6)
def rmsf():
    plt.figure()

    plt.plot(f1["Residue"], f1["RMSF"], label='T232D')
    plt.plot(f2["Residue"], f2["RMSF"], label='WT')
    plt.plot(f3["Residue"], f3["RMSF"], label='k215r_t232d')
    plt.plot(f4["Residue"], f4["RMSF"], label='k215r')
    plt.plot(f5["Residue"], f5["RMSF"], label='t232a_94')
    plt.plot(f6["Residue"], f6["RMSF"], label='t232a_73')

    plt.xlim(0, 150)
    plt.axvspan(142, 144, alpha=0.2, color="lime", label='DFG')
    plt.axvspan(8, 13, alpha=0.2, color="red", label='G-loop')
    plt.axvspan(43, 53, alpha=0.2, color="yellow", label="αC helix")

    plt.xlabel('Residue number')
    plt.ylabel('RMSF (Å)')
    plt.title("Aurora Kinase RMSF Comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{range}_rmsf_{ref}.png", dpi=300)

rmsd()
#rmsd_wo_t232a()
#naming= range, what, reference
plt.show()
