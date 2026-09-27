import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

data = [
    {
        "name": "k215r_t232d_run1",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run1_RMSF.csv"
    }, 
    {
        "name": "k215r_t232d_run2",
        "file": "/Volumes/Jason/Simulations/RMSF/0_curr/3_k215r_t232d_npt_run2_RMSF.csv"
    }
]

RESIDUE_COL = "Residue"
RMSF_COL    = "RMSF_Angstrom"

arjun = pd.read_csv("/Volumes/Jason/Simulations/RMSF/arjun/Arjun_data_combined.csv")

# ── Build Arjun's combined DataFrame ─────────────────────────────────────────
arjun_runs = {}
for rep, color in zip(['Rep1_Mut', 'Rep2_Mut', 'Rep3_Mut'],
                      ['burlywood', 'navajowhite', 'moccasin']):
    df = arjun[["Residue", rep]].copy()
    df.columns = [RESIDUE_COL, RMSF_COL]
    df[RESIDUE_COL] = df[RESIDUE_COL] + 76
    arjun_runs[rep] = df

# Stack all Arjun RMSF values per residue and compute average
arjun_stacked = pd.concat(list(arjun_runs.values()))
arjun_avg = arjun_stacked.groupby(RESIDUE_COL)[RMSF_COL].mean().reset_index()

# ── Load your data ────────────────────────────────────────────────────────────
n      = len(data)
colors = cm.tab10(np.linspace(0, 1, min(n, 10)))

my_dfs = []
for entry in data:
    df = pd.read_csv(entry["file"])
    df[RESIDUE_COL] = df[RESIDUE_COL] + 76
    my_dfs.append(df)

my_stacked = pd.concat(my_dfs)
my_avg = my_stacked.groupby(RESIDUE_COL)[RMSF_COL].mean().reset_index()

# ── Compute min/max across ALL values (arjun + mine) per residue for shading ─
all_data = pd.concat([arjun_stacked, my_stacked])
all_min  = all_data.groupby(RESIDUE_COL)[RMSF_COL].min()
all_max  = all_data.groupby(RESIDUE_COL)[RMSF_COL].max()
residues = all_min.index

# Per-group min/max for individual highlights
arjun_min_val = arjun_avg[RMSF_COL].min()
arjun_max_val = arjun_avg[RMSF_COL].max()
arjun_min_res = arjun_avg.loc[arjun_avg[RMSF_COL].idxmin(), RESIDUE_COL]
arjun_max_res = arjun_avg.loc[arjun_avg[RMSF_COL].idxmax(), RESIDUE_COL]

my_min_val = my_avg[RMSF_COL].min()
my_max_val = my_avg[RMSF_COL].max()
my_min_res = my_avg.loc[my_avg[RMSF_COL].idxmin(), RESIDUE_COL]
my_max_res = my_avg.loc[my_avg[RMSF_COL].idxmax(), RESIDUE_COL]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, (ax, ax_avg) = plt.subplots(2, 1, figsize=(16, 14), sharex=True,
                                 gridspec_kw={'hspace': 0.35})

# ── Top panel: individual runs ────────────────────────────────────────────────
# Shaded envelope (global min–max across all runs)
ax.fill_between(residues, all_min, all_max, alpha=0.12, color='gray',
                label='Min–Max envelope (all runs)')

rep_labels = ['k215r_t232d_run1_charm', 'k215r_t232d_run2_charm', 'k215r_t232d_run3_charm']
rep_colors = ['burlywood', 'navajowhite', 'moccasin']
for (rep, df), lbl, col in zip(arjun_runs.items(), rep_labels, rep_colors):
    ax.plot(df[RESIDUE_COL], df[RMSF_COL], linewidth=1.2, label=lbl, color=col)

for entry, color in zip(data, colors):
    df = pd.read_csv(entry["file"])
    df[RESIDUE_COL] = df[RESIDUE_COL] + 76
    ax.plot(df[RESIDUE_COL], df[RMSF_COL], linewidth=1.4,
            label=entry["name"], color=color)

ax.set_ylabel("RMSF (Å)", fontsize=12)
ax.set_title("Per-Residue RMSF — Individual Runs", fontsize=14, fontweight="bold")
ax.legend(fontsize=9, framealpha=0.7, ncol=2)
ax.grid(True, linestyle="--", alpha=0.4)
ax.axvspan(142+76, 144+76, alpha=0.2, color="lime")
ax.axvspan(8+76,   13+76,  alpha=0.2, color="red")
ax.axvspan(43+76,  53+76,  alpha=0.2, color="yellow")

# ── Bottom panel: averages with min/max highlights ────────────────────────────
ax_avg.plot(arjun_avg[RESIDUE_COL], arjun_avg[RMSF_COL],
            linewidth=2, color='saddlebrown', label="Arjun avg (CHARMM)")
ax_avg.plot(my_avg[RESIDUE_COL], my_avg[RMSF_COL],
            linewidth=2, color='steelblue', label="My avg (AMBER/GROMACS)")

# Shaded region between the two averages (difference)
common_res = np.intersect1d(arjun_avg[RESIDUE_COL], my_avg[RESIDUE_COL])
a_interp = arjun_avg.set_index(RESIDUE_COL).reindex(common_res)[RMSF_COL]
m_interp = my_avg.set_index(RESIDUE_COL).reindex(common_res)[RMSF_COL]
ax_avg.fill_between(common_res, a_interp, m_interp,
                    where=(a_interp >= m_interp), alpha=0.2, color='saddlebrown',
                    label='Arjun > Mine')
ax_avg.fill_between(common_res, a_interp, m_interp,
                    where=(m_interp > a_interp),  alpha=0.2, color='steelblue',
                    label='Mine > Arjun')

# ── Arjun min/max markers ─────────────────────────────────────────────────────
ax_avg.scatter(arjun_min_res, arjun_min_val, zorder=5, s=90,
               color='saddlebrown', marker='v',
               label=f"Arjun min: {arjun_min_val:.2f} Å @ res {arjun_min_res}")
ax_avg.scatter(arjun_max_res, arjun_max_val, zorder=5, s=90,
               color='saddlebrown', marker='^',
               label=f"Arjun max: {arjun_max_val:.2f} Å @ res {arjun_max_res}")
ax_avg.annotate(f"min {arjun_min_val:.2f}", xy=(arjun_min_res, arjun_min_val),
                xytext=(5, -14), textcoords='offset points',
                fontsize=8, color='saddlebrown')
ax_avg.annotate(f"max {arjun_max_val:.2f}", xy=(arjun_max_res, arjun_max_val),
                xytext=(5, 6), textcoords='offset points',
                fontsize=8, color='saddlebrown')

# ── My min/max markers ────────────────────────────────────────────────────────
ax_avg.scatter(my_min_res, my_min_val, zorder=5, s=90,
               color='steelblue', marker='v',
               label=f"My min: {my_min_val:.2f} Å @ res {my_min_res}")
ax_avg.scatter(my_max_res, my_max_val, zorder=5, s=90,
               color='steelblue', marker='^',
               label=f"My max: {my_max_val:.2f} Å @ res {my_max_res}")
ax_avg.annotate(f"min {my_min_val:.2f}", xy=(my_min_res, my_min_val),
                xytext=(5, -14), textcoords='offset points',
                fontsize=8, color='steelblue')
ax_avg.annotate(f"max {my_max_val:.2f}", xy=(my_max_res, my_max_val),
                xytext=(5, 6), textcoords='offset points',
                fontsize=8, color='steelblue')

# ── Annotation box: overall average RMSF comparison ──────────────────────────
arjun_global_avg = arjun_avg[RMSF_COL].mean()
my_global_avg    = my_avg[RMSF_COL].mean()
diff             = my_global_avg - arjun_global_avg
textstr = (f"Global mean RMSF\n"
           f"  Arjun : {arjun_global_avg:.3f} Å\n"
           f"  Mine  : {my_global_avg:.3f} Å\n"
           f"  Δ     : {diff:+.3f} Å")
props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
ax_avg.text(0.02, 0.97, textstr, transform=ax_avg.transAxes,
            fontsize=10, verticalalignment='top', bbox=props)

ax_avg.set_xlabel("Residue", fontsize=12)
ax_avg.set_ylabel("RMSF (Å)", fontsize=12)
ax_avg.set_title("Average RMSF Comparison: Arjun vs Mine", fontsize=14, fontweight="bold")
ax_avg.legend(fontsize=8, framealpha=0.7, ncol=2)
ax_avg.grid(True, linestyle="--", alpha=0.4)
ax_avg.axvspan(142+76, 144+76, alpha=0.2, color="lime",   label='DFG')
ax_avg.axvspan(8+76,   13+76,  alpha=0.2, color="red",    label='G-loop')
ax_avg.axvspan(43+76,  53+76,  alpha=0.2, color="yellow", label="αC helix")

fig.tight_layout()
plt.savefig("/Volumes/Jason/Simulations/RMSF/rmsf_comparison.png", dpi=150, bbox_inches='tight')
plt.show()