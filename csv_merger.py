import pandas as pd
# file_list=["/Volumes/Jason/Simulations/vs_self/30_49/recorded/trajectory_k215r_64128_r88_v_self_rem0_30_res49_dist.csv",
#            "/Volumes/Jason/Simulations/vs_self/30_49/recorded/trajectory_k215r_64128_r88_v_self_rem1_30_res49_dist.csv", 
#            "/Volumes/Jason/Simulations/vs_self/30_49/recorded/trajectory_k215r_64128_r88_v_self_rem2_30_res49_dist.csv", 
#            "/Volumes/Jason/Simulations/vs_self/30_49/recorded/trajectory_k215r_64128_r88_v_self_rem3_30_res49_dist.csv"]

#/Volumes/Jason/Simulations/Feb2026/self/trajectory_AukB_k215r_t232d_64_128_19_11_vs_self_rem_0_30_res49_dist.csv

file_list=["/Volumes/Jason/Simulations/4March2026/total/done/trajectory_t232a_94_resume_0_v_self_r1_99_rmsd.csv", 
           "/Volumes/Jason/Simulations/4March2026/total/trajectory_t232a_94_v_self_r2_rmsd.csv"]
#/Volumes/Jason/Simulations/4March2026/160180/trajectory_t232a_94_v_active_r2_99_160_180_rmsd.csv
output_path = "/Volumes/Jason/Simulations/4March2026/total/00_trajectory_t232a_94_v_self_dist_total.csv"

def merge_recursive(file_list, output_path, time_col="Time_ns", idx=0, accumulated_df=None):

    # Base case: first file
    if idx == 0:
        accumulated_df = pd.read_csv(file_list[0]).copy()
        return merge_recursive(file_list, output_path, time_col, idx=1, accumulated_df=accumulated_df)

    # Base case: finished recursion
    if idx >= len(file_list):
        accumulated_df.to_csv(output_path, index=False)
        print("Final merged file written successfully.")
        print("Total frames:", len(accumulated_df))
        return accumulated_df

    # Recursive step
    print(f"Merging file {idx+1}/{len(file_list)}")

    next_df = pd.read_csv(file_list[idx]).copy()

    # Remove overlap frames
    accumulated_df = accumulated_df.iloc[:-1].copy()   # drop last
    next_df = next_df.iloc[1:].copy()                   # drop first

    # Shift time
    last_time = accumulated_df[time_col].iloc[-1]
    next_df.loc[:, time_col] += last_time

    # Merge
    accumulated_df = pd.concat([accumulated_df, next_df], ignore_index=True)

    # Recurse
    return merge_recursive(file_list, output_path, time_col, idx+1, accumulated_df)

merge_recursive(file_list, output_path)



