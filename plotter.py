import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

#data in format of [location, name]
data_total_single = [
    {
        "name": "wt_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/0_wt_npt_run1_RMSF.csv"
    },
    {
        "name": "wt_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/0_wt_npt_run2_RMSF.csv"
    }, 
    {
        "name": "t232d_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/1_t232d_npt_run1_RMSF.csv"
    }, 
    {
        "name": "t232d_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/trajectory_t232a_94_self_RMSF.csv"
    }, 
    {
        "name": "k215r_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/2_k215r_npt_run1_RMSF.csv"
    }, 
    {
        "name": "k215r_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/2_k215r_npt_run2_RMSF.csv"
    }, 
    {
        "name": "k215r_t232d_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run1_RMSF.csv"
    }, 
    {
        "name": "k215r_t232d_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run2_RMSF.csv"
    }, 
    {
        "name": "t232a_73_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/4_t232a_73_npt_run1_RMSF.csv"
    }, 
    {
        "name": "t232a_73_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/4_t232a_73_npt_run2_RMSF.csv"
    }, 
    {
        "name": "t232a_94_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/5_t232a_94_npt_run1_RMSF.csv"
    }, 
    {
        "name": "t232a_94_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/5_t232a_94_npt_run2_RMSF.csv"
    }
]
data_total_average = [
    {
        "name": "wt_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/wt_npt_run1_average_RMSF.csv"
    },
    {
        "name": "wt_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/wt_npt_run2_average_RMSF.csv"
    }, 
    {
        "name": "t232d_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232d_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "t232d_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232d_npt_run2_average_RMSF.csv"
    }, 
    {
        "name": "k215r_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "k215r_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_npt_run2_average_RMSF.csv"
    }, 
    {
        "name": "k215r_t232d_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_t232d_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "k215r_t232d_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_t232d_npt_run2_average_RMSF.csv"
    }, 
    {
        "name": "t232a_73_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232a_73_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "t232a_73_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232a_73_npt_run2_average_RMSF.csv"
    }, 
    {
        "name": "t232a_94_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232a_94_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "t232a_94_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/t232a_94_npt_run2_average_RMSF.csv"
    }
]
data = [
    {
        "name": "t232d_k215r_run1_first",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run1_RMSF.csv"
    }, 
    {
        "name": "t232d_k215r_run2_first",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run2_RMSF.csv"
    },
    {
        "name": "t232d_k215r_run1_average",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_t232d_npt_run1_average_RMSF.csv"
    }, 
    {
        "name": "t232d_k215r_run1_average",
        "file": "/Volumes/Jason/Simulations/RMSF/0_rmsf_average/k215r_t232d_npt_run2_average_RMSF.csv"
    }
]
data_residue = [
    {
        "name": "wt_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/wt_r1t2_vs_AukB_wt_3264_r46_2_008_30NZ-49OE1OE2.csv"
    },
    {
        "name": "wt_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/wt_r2t2_vs_AukB_wt_3264_r46_2_008_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232d_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232d_r1t2_vs_AurB_t232d_3_007_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232d_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232d_r2t2_vs_AurB_t232d_3_007_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "k215r_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/k215r_r1t2_vs_k215r_64128_r88_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "k215r_run_3",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/k215r_r3t2_vs_k215r_64128_r88_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "k215r_t232d_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/k215r_t232d_r1t2_vs_AukB_k215r_t232d_64_128_19_11_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "k215r_t232d_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/k215r_t232d_r2t2_vs_AukB_k215r_t232d_64_128_19_11_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "k215r_t232d_run_3",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/k215r_t232d_r3t2_vs_AukB_k215r_t232d_64_128_19_11_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232a_73_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232a_73_r1t2_vs_t232a_73_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232a_73_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232a_73_r2t2_vs_t232a_73_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232a_94_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232a_94_r1t2_vs_t232a_94_30NZ-49OE1OE2.csv"
    }, 
    {
        "name": "t232a_94_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/t232a_94_r2t2_vs_t232a_94_30NZ-49OE1OE2.csv"
    },
    {
        "name": "4af3_incenp_run_1",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/4af3_incenp_r1t2_vs_4af3_incenp_box_106NZ-125OE1OE2.csv"
    }, 
    {
        "name": "4af3_incenp_run_2",
        "file": "/Volumes/Jason/Simulations/Data/csv/term_amide_K106_term_carboxyl_E125/4af3_incenp_r2t2_vs_4af3_incenp_box_106NZ-125OE1OE2.csv"
    }
]
#print(data_residue[1]["name"])

complex_data=[
    # {
    #     "name": "4af3_pure_r1t2",
    #     "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r1t2_vs_4af3_pure_box_77_338.csv"
    # },
    # {
    #     "name": "4af3_pure_r2t2",
    #     "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r2t2_vs_4af3_pure_box_77_338.csv"
    # }, 
    # {
    #     "name": "4af3_pure_r3t2",
    #     "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r3t2_vs_4af3_pure_box_77_338.csv"
    # }
    {
        "name": "2bfx_incenp_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r1t2_vs_2bfx_box_78_357.csv"
    },
    {
        "name": "2bfx_incenp_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r2t2_vs_2bfx_box_78_357.csv"
    },
    {
        "name": "2bfx_incenp_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r3t2_vs_2bfx_box_78_357.csv"
    }
]

af3_files=[
    {
        "name": "4af3_AurB_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r1t2_vs_4af3_pure_box_77_338.csv"
    },
    {
        "name": "4af3_AurB_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r2t2_vs_4af3_pure_box_77_338.csv"
    }, 
    {
        "name": "4af3_AurB_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/4af3_pure_r3t2_vs_4af3_pure_box_77_338.csv"
    }, 
    {
        "name": "4af3_incenp_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r1t2_vs_4af3_incenp_box_77_338.csv"
    },
    {
        "name": "4af3_incenp_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r2t2_vs_4af3_incenp_box_77_338.csv"
    },
    {
        "name": "4af3_incenp_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r3t2_vs_4af3_incenp_box_77_338.csv"
    }
]

bfx_files=[
    {
        "name": "2bfx_AurB_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r1t2_vs_2bfx_box_78_357.csv"
    },
    {
        "name": "2bfx_AurB_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r2t2_vs_2bfx_box_78_357.csv"
    }, 
    {
        "name": "2bfx_AurB_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/whole/2bfx_incenp_r3t2_vs_2bfx_box_78_357.csv"
    },
    {
        "name": "2bfx_incenp_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r1t2_vs_2bfx_AurB_box_78_357.csv"
    },
    {
        "name": "2bfx_incenp_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r2t2_vs_2bfx_AurB_box_78_357.csv"
    },
    {
        "name": "2bfx_incenp_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r3t2_vs_2bfx_AurB_box_78_357.csv"
    }
]

complex_data_1=[
    {
        "name": "2bfx_pure_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r1t2_vs_2bfx_AurB_box_78_357.csv"
    },
    {
        "name": "2bfx_pure_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r2t2_vs_2bfx_AurB_box_78_357.csv"
    }, 
    {
        "name": "2bfx_pure_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/2bfx_AurB_r3t2_vs_2bfx_AurB_box_78_357.csv"
    },
    {
        "name": "4af3_incenp_r1t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r1t2_vs_4af3_incenp_box_77_338.csv"
    },
    {
        "name": "4af3_incenp_r2t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r2t2_vs_4af3_incenp_box_77_338.csv"
    },
    {
        "name": "4af3_incenp_r3t2",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/complex_07302026/4af3_incenp_r3t2_vs_4af3_incenp_box_77_338.csv"
    }
]
averages=[
    {
        "name": "2bfx_average_complex",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/DF_218_219/2bfx_incenp_average.csv"
    }, 
    {
        "name": "4af3_average",
        "file": "/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/DF_218_219/4af3_pure_average.csv"
    }
]
RESIDUE_COL = "Time_ns"
#RMSF_COL    = "Distance_A"
RMSF_COL = "RMSD_A"
#RMSF_COL="mean"
#arjun=pd.read_csv("/Volumes/Jason/Simulations/sim2/trial2/complex/07202026/DF_218_219/2bfx_incenp_r1t2_vs_2bfx_box_234_235.csv")
#print(arjun)

# # ── Plot ──────────────────────────────────────────────────────────────────────

n = len(data_residue)

colors = (
    cm.tab20(np.linspace(0, 1, n))
    if n <= 20
    else cm.turbo(np.linspace(0, 1, n))
)
fig, ax = plt.subplots(figsize=(15, 8))

def sum():
    arjun_ace_phos_1 = arjun[["Residue", 'Rep1_Mut']].copy()
    arjun_ace_phos_1.columns = [RESIDUE_COL, RMSF_COL]
    arjun_ace_phos_1[RESIDUE_COL] = arjun_ace_phos_1[RESIDUE_COL] + 76
    arjun_ace_phos_1 = arjun_ace_phos_1[arjun_ace_phos_1[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_ace_phos_1[RESIDUE_COL], arjun_ace_phos_1[RMSF_COL], linewidth=1.4, label="k215r_t232d_run1_charm", color="burlywood")

    arjun_ace_phos_2 = arjun[["Residue", 'Rep2_Mut']].copy()
    arjun_ace_phos_2.columns = [RESIDUE_COL, RMSF_COL]
    arjun_ace_phos_2[RESIDUE_COL] = arjun_ace_phos_2[RESIDUE_COL] + 76
    arjun_ace_phos_2 = arjun_ace_phos_2[arjun_ace_phos_2[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_ace_phos_2[RESIDUE_COL], arjun_ace_phos_2[RMSF_COL], linewidth=1.4, label="k215r_t232d_run2_charm", color="navajowhite")

    arjun_ace_phos_3 = arjun[["Residue", 'Rep3_Mut']].copy()
    arjun_ace_phos_3.columns = [RESIDUE_COL, RMSF_COL]
    arjun_ace_phos_3[RESIDUE_COL] = arjun_ace_phos_3[RESIDUE_COL] + 76
    arjun_ace_phos_3 = arjun_ace_phos_3[arjun_ace_phos_3[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_ace_phos_3[RESIDUE_COL], arjun_ace_phos_3[RMSF_COL], linewidth=1.4, label="k215r_t232d_run3_charm", color="moccasin")


    arjun_phosp_1 = arjun[["Residue", 'Rep1_WT']].copy()
    arjun_phosp_1.columns = [RESIDUE_COL, RMSF_COL]
    arjun_phosp_1[RESIDUE_COL] = arjun_phosp_1[RESIDUE_COL] + 76
    arjun_phosp_1 = arjun_phosp_1[arjun_phosp_1[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_phosp_1[RESIDUE_COL], arjun_phosp_1[RMSF_COL], linewidth=1.4, label="t232d_run1_charm", color="orange")

    arjun_phosp_2 = arjun[["Residue", 'Rep2_WT']].copy()
    arjun_phosp_2.columns = [RESIDUE_COL, RMSF_COL]
    arjun_phosp_2[RESIDUE_COL] = arjun_phosp_2[RESIDUE_COL] + 76
    arjun_phosp_2 = arjun_phosp_2[arjun_phosp_2[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_phosp_2[RESIDUE_COL], arjun_phosp_2[RMSF_COL], linewidth=1.4, label="t232d_run2_charm", color="green")

    arjun_phosp_3 = arjun[["Residue", 'Rep3_WT']].copy()
    arjun_phosp_3.columns = [RESIDUE_COL, RMSF_COL]
    arjun_phosp_3[RESIDUE_COL] = arjun_phosp_3[RESIDUE_COL] + 76
    arjun_phosp_3 = arjun_phosp_3[arjun_phosp_3[RESIDUE_COL] <= 250+76]
    ax.plot(arjun_phosp_3[RESIDUE_COL], arjun_phosp_3[RMSF_COL], linewidth=1.4, label="k215r_t232d_run3_charm", color="red")
import pandas as pd

def compute_group_average(entries, time_col=RESIDUE_COL, val_col=RMSF_COL):
    dfs = []
    for entry in entries:
        df = pd.read_csv(entry["file"])[[time_col, val_col]]
        dfs.append(df.set_index(time_col)[val_col])
    # align on Time_ns and average across replicates row-wise
    merged = pd.concat(dfs, axis=1)
    avg = merged.mean(axis=1).reset_index()
    avg.columns = [time_col, val_col]
    return avg

# # split complex_data by system prefix
groups = {}
for entry in complex_data_1:
    key = "2bfx" if entry["name"].startswith("2bfx") else "4af3"
    groups.setdefault(key, []).append(entry)

avg_2bfx = compute_group_average(groups["2bfx"])
avg_4af3 = compute_group_average(groups["4af3"])

# # ── Plot ──────────────────────────────────────────────────────────────────────

n = len(data_residue)

colors = (
    cm.tab20(np.linspace(0, 1, n))
    if n <= 20
    else cm.turbo(np.linspace(0, 1, n))
)
fig, ax = plt.subplots(figsize=(15, 8))

# #add 76 for FastConformation stuff
# for entry, color in zip(complex_data_1, colors):
#     df = pd.read_csv(entry["file"])
#     print(entry["file"])
#     print(df)
#     print(entry["name"], df[RMSF_COL].mean())
#     ax.plot(df[RESIDUE_COL], df[RMSF_COL], linewidth=1.4,label=entry["name"], color=color)
# ax.plot(avg_2bfx[RESIDUE_COL], avg_2bfx[RMSF_COL], linewidth=2, label="2bfx_average", color="black", alpha=0.5)
# ax.plot(avg_4af3[RESIDUE_COL], avg_4af3[RMSF_COL], linewidth=2, label="4af3_average", color="purple",  alpha=0.5)

for entry, color in zip(af3_files, colors):
    df = pd.read_csv(entry["file"])
    print(entry["file"])
    print(df)
    print(entry["name"], df[RMSF_COL].mean())
    ax.plot(df[RESIDUE_COL], df[RMSF_COL], linewidth=1.4,label=entry["name"], color=color)

ax.set_xlabel("time (ns)", fontsize=12)
ax.set_ylabel("RMSD (Å)", fontsize=12)
ax.set_title("Whole", fontsize=14, fontweight="bold")
ax.legend(fontsize=10, framealpha=0.7)
ax.grid(True, linestyle="--", alpha=0.4)

#print("2bfx average:", avg_2bfx[RMSF_COL].mean())
#print("4af3 average:", avg_4af3[RMSF_COL].mean())

# plt.axvspan(142+76, 144+76, alpha=0.2, color="lime", label='DFG')
# plt.axvspan(8+76, 13+76, alpha=0.2, color="red", label='G-loop')
# plt.axvspan(43+76, 53+76, alpha=0.2, color="yellow", label="αC helix")
fig.tight_layout()
plt.show()
