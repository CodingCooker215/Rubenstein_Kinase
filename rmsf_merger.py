import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df1=pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_k215r_self_RMSF.csv")
df2=pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_t232d_self_RMSF.csv")
df3=pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_k215rt_t232d_self_RMSF.csv")
df4=pd.read_csv("/Volumes/Jason/Simulations/RMSF/trajectory_t232a_73_self_RMSF.csv")
#df3=pd.read_csv("/Volumes/Jason/Simulations/4March2026/ttrajectory_t232a_73_resume_0_v_self_r1_99_RMSF.csv")
#df4=pd.read_csv("/Volumes/Jason/Simulations/4March2026/trajectory_t232a_94_resume_0_v_self_r1_99_RMSF.csv")


def standardize(df):
    df.columns = df.columns.str.strip()

    if "Residue" in df.columns:
        df = df.rename(columns={"Residue":"X_values"})
    if "RMSF" in df.columns:
        df = df.rename(columns={"RMSF":"Y_values"})
    if "RMSF_Angstrom" in df.columns:
        df = df.rename(columns={"RMSF_Angstrom":"Y_values"})

    return df


df1 = standardize(df1)
df2 = standardize(df2)
df3 = standardize(df3)
df4 = standardize(df4)


def annotate(offset=76):
    plt.axvspan(142+offset,144+offset,alpha=0.2,color='orange',label='DFG')
    plt.axvspan(8+offset,13+offset,alpha=0.2,color='green',label='G-loop')
    plt.axvspan(43+offset,53+offset,alpha=0.2,color='yellow',label="αC helix")


def calculate_sem(df):
    return df.groupby("X_values")["Y_values"].agg(['mean','sem']).reset_index()


sem1 = calculate_sem(df1[df1["X_values"]<=250])
sem2 = calculate_sem(df2[df2["X_values"]<=250])
sem3 = calculate_sem(df3[df3["X_values"]<=250])
sem4 = calculate_sem(df4[df4["X_values"]<=250])


def plot_sem(data, color, label, offset=76):
    x = data["X_values"] + offset
    mean = data["mean"]
    sem = data["sem"]

    plt.plot(x, mean, color=color, label=label)
    plt.fill_between(x, mean-sem, mean+sem, color=color, alpha=0.25)


plt.figure(figsize=(10,5))

plot_sem(sem1, "blue", "K215R")
plot_sem(sem2, "red", "T232D")
plot_sem(sem3, "green", "K215R+T232D")
plot_sem(sem4, "purple", "T232A_73")

annotate()

plt.xlabel("Residue Number")
plt.ylabel("RMSF (Å)")
plt.legend()
plt.tight_layout()

plt.show()