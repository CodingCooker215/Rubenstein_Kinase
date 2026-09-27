import csv
import os
from pymol import cmd

csv_folder = "/Volumes/Jason/Simulations/sim2/trial2/complex/"

# Color scheme per chain: A=protein, B=ATP, C=compound, D=Mg2+
chain_colors = {
    "A": "orange",
    "B": "yellow",
    "C": "green",
    "D": "purple",
}

pdb_objects = cmd.get_names("objects")
cmd.set("cartoon_transparency", 0.85)

for obj_name in pdb_objects:
    file_name = f"{obj_name}_binding_residues.csv"
    csv_path = os.path.join(csv_folder, file_name)

    if not os.path.exists(csv_path):
        print(f"Skipping {obj_name} (CSV not found)")
        continue

    # Group residue numbers by chain so we don't mix up e.g. A/12 with D/12
    chain_residues = {}
    with open(csv_path, mode="r") as file:
        reader = csv.reader(file)
        next(reader)  # skip header row: Side, Chain, Residue_Code, Residue_Number
        for row in reader:
            side, chain_id, res_code, res_num = row
            chain_residues.setdefault(chain_id, []).append(res_num)

    if not chain_residues:
        print(f"No residues found in {file_name}")
        continue

    for chain_id, resnums in chain_residues.items():
        joined_numbers = "+".join(resnums)
        sele_name = f"{obj_name}_{chain_id}_contact"
        cmd.select(sele_name, f"{obj_name} and chain {chain_id} and resi {joined_numbers}")
        cmd.show("sticks", sele_name)
        cmd.set("stick_transparency", 0, sele_name)
        cmd.color(chain_colors.get(chain_id, "white"), sele_name)

    print(f"{obj_name}: highlighted chains {list(chain_residues.keys())}")