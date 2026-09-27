import csv
from pymol import cmd

# Dictionary to map 3-letter residue names to 1-letter codes
aa_mapping = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
    'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
    'SEC':'U', 'PYL':'O'
}

# 1. Get the names of all loaded PDB objects
pdb_objects = cmd.get_object_list("all")

# 2. Loop through each PDB structure
for pdb in pdb_objects:

    # --- Side 1: Chain A residues within 4A of Chain D ---
    site_A = f"{pdb}_site_A"
    cmd.select(site_A, f"byres ({pdb} and chain A within 3 of chain D)")
    residues_A = set()
    cmd.iterate(site_A, "residues_A.add((chain, resn, resi))", space={'residues_A': residues_A})

    # --- Side 2: Chain D residues within 4A of Chain A ---
    site_D = f"{pdb}_site_D"
    cmd.select(site_D, f"byres ({pdb} and chain D within 3 of chain A)")
    residues_D = set()
    cmd.iterate(site_D, "residues_D.add((chain, resn, resi))", space={'residues_D': residues_D})

    if residues_A or residues_D:
        filename = f"/Volumes/Jason/Simulations/sim2/trial2/complex/{pdb}_binding_residues.csv"
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Side", "Chain", "Residue_Code", "Residue_Number"])

            # Write Chain A side, sorted by residue number
            for chain_id, res_name, res_num in sorted(residues_A, key=lambda x: int(x[2])):
                one_letter_code = aa_mapping.get(res_name.upper(), res_name)
                writer.writerow(["A_contacts_D", chain_id, one_letter_code, res_num])

            # Write Chain D side, sorted by residue number
            for chain_id, res_name, res_num in sorted(residues_D, key=lambda x: int(x[2])):
                one_letter_code = aa_mapping.get(res_name.upper(), res_name)
                writer.writerow(["D_contacts_A", chain_id, one_letter_code, res_num])

        print(f"Successfully created: {filename}")
    else:
        print(f"No binding residues found for object: {pdb}")

    # Clean up selections to avoid clutter across loop iterations
    #cmd.delete(site_A)
    #cmd.delete(site_D)