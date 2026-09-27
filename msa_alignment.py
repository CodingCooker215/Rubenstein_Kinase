"""
Align Aurora Kinase (2bfx chain B / 4af3 chain A) and INCENP (2bfx chain C / 4af3 chain D)
sequences pairwise, then map each structure's contact-residue annotations (from the
*_binding_residues.csv files) onto the alignment so the two structures' binding-site
residues can be compared column-by-column.

Why a "monotonic offset" residue mapper is needed:
    PDB residue numbers usually do NOT correspond 1:1 with the index of the sequence in
    the *_box.csv files. Crystal structures almost always have a few residues that are
    disordered / not resolved in the density, so the author residue numbering has small
    jumps (gaps) that don't exist in the full construct sequence used for AlphaFold /
    SEQRES. This script reconstructs the residue-number -> sequence-index mapping
    empirically per chain from the contact list itself, assuming the offset
    (resnum - index) only ever increases (never decreases) as you move along the chain
    -- true so long as the PDB numbering has no insertions relative to the sequence,
    only deletions (missing/unresolved residues).
"""

import csv
import os
from itertools import groupby

from Bio.Align import PairwiseAligner, substitution_matrices

DATA_DIR = "/Volumes/Jason/Simulations/sim2/trial2/complex"
OUT_DIR = "/Volumes/Jason/Simulations/sim2/trial2/complex"
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Load full-length sequences from the *_box.csv files
# ---------------------------------------------------------------------------
def load_sequences(path):
    seqs = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            seqs[row["Object"]] = row["Sequence"].strip()
    return seqs


# ---------------------------------------------------------------------------
# 2. Load contact residues from the *_binding_residues.csv files, grouped by chain
# ---------------------------------------------------------------------------
def load_contacts_by_chain(path):
    by_chain = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            chain = row["Chain"]
            resnum = int(row["Residue_Number"])
            code = row["Residue_Code"]
            side = row["Side"]
            by_chain.setdefault(chain, []).append((resnum, code, side))
    for chain in by_chain:
        by_chain[chain].sort(key=lambda x: x[0])
    return by_chain


# ---------------------------------------------------------------------------
# 3. Map PDB residue numbers -> 0-based index in the box.csv sequence
#    using a monotonic-offset greedy matcher (see module docstring)
# ---------------------------------------------------------------------------
def map_resnum_to_index(contacts, seq):
    """
    contacts: list of (resnum, code, side), sorted ascending by resnum
    seq: full-length sequence string
    returns: dict resnum -> 0-based index into seq (only for resolvable contacts)
    """
    mapping = {}
    current_offset = None  # offset = resnum - index

    for resnum, code, side in contacts:
        candidates = [i for i, c in enumerate(seq) if c == code]
        if not candidates:
            continue

        if current_offset is None:
            # first contact: take the smallest plausible offset
            chosen = min(candidates, key=lambda i: resnum - i)
        else:
            # prefer the candidate whose offset is >= current_offset and closest to it
            valid = [i for i in candidates if (resnum - i) >= current_offset]
            if valid:
                chosen = min(valid, key=lambda i: (resnum - i) - current_offset)
            else:
                # fallback: closest candidate overall (shouldn't normally happen)
                chosen = min(candidates, key=lambda i: abs((resnum - i) - current_offset))

        current_offset = resnum - chosen
        mapping[resnum] = chosen

    return mapping


# ---------------------------------------------------------------------------
# 4. Pairwise global alignment (Biopython), BLOSUM62 + affine gaps
# ---------------------------------------------------------------------------
def align_sequences(seq_a, seq_b):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    alignment = aligner.align(seq_a, seq_b)[0]
    return alignment


def alignment_to_columns(alignment):
    """
    Returns two parallel lists (aligned_a, aligned_b) of characters,
    '-' for gaps, plus index maps: seq-index -> alignment-column for each side.
    """
    aligned_a, aligned_b = alignment.aligned, None
    a_str, b_str = str(alignment).split("\n")[0], str(alignment).split("\n")[2]
    # Bio.Align's default __str__ wraps; safer to rebuild from alignment coordinates
    idx_a, idx_b = 0, 0
    col_a, col_b = [], []
    idx_to_col_a, idx_to_col_b = {}, {}
    for col in range(alignment.shape[1]):
        # use alignment.indices for robust mapping (Bio >=1.80)
        pass
    return idx_to_col_a, idx_to_col_b


def build_index_to_column_maps(alignment):
    """Use alignment.indices (2 x n_columns array) to map sequence-index -> column."""
    indices = alignment.indices  # shape (2, n_columns); -1 marks a gap
    idx_to_col_a, idx_to_col_b = {}, {}
    for col in range(indices.shape[1]):
        a_idx, b_idx = indices[0][col], indices[1][col]
        if a_idx != -1:
            idx_to_col_a[a_idx] = col
        if b_idx != -1:
            idx_to_col_b[b_idx] = col
    return idx_to_col_a, idx_to_col_b


# ---------------------------------------------------------------------------
# 5. Build the combined annotated alignment CSV for one protein pair
# ---------------------------------------------------------------------------
def build_annotated_alignment(
    label,
    obj_a, seq_a, contacts_a,
    obj_b, seq_b, contacts_b,
    out_path,
):
    alignment = align_sequences(seq_a, seq_b)
    idx_to_col_a, idx_to_col_b = build_index_to_column_maps(alignment)
    indices = alignment.indices
    n_cols = indices.shape[1]

    resmap_a = map_resnum_to_index(contacts_a, seq_a)
    resmap_b = map_resnum_to_index(contacts_b, seq_b)

    # invert: column -> (resnum, code, side) if that column is a contact for that object
    col_annot_a = {}
    for resnum, idx in resmap_a.items():
        if idx in idx_to_col_a:
            col = idx_to_col_a[idx]
            code = next(c for r, c, s in contacts_a if r == resnum)
            side = next(s for r, c, s in contacts_a if r == resnum)
            col_annot_a[col] = (resnum, code, side)

    col_annot_b = {}
    for resnum, idx in resmap_b.items():
        if idx in idx_to_col_b:
            col = idx_to_col_b[idx]
            code = next(c for r, c, s in contacts_b if r == resnum)
            side = next(s for r, c, s in contacts_b if r == resnum)
            col_annot_b[col] = (resnum, code, side)

    rows = []
    for col in range(n_cols):
        a_idx, b_idx = indices[0][col], indices[1][col]
        a_char = seq_a[a_idx] if a_idx != -1 else "-"
        b_char = seq_b[b_idx] if b_idx != -1 else "-"
        a_resnum, a_side = ("", "")
        b_resnum, b_side = ("", "")
        if col in col_annot_a:
            a_resnum, _, a_side = col_annot_a[col]
        if col in col_annot_b:
            b_resnum, _, b_side = col_annot_b[col]

        conserved = (a_char == b_char and a_char != "-")
        either_contact = (col in col_annot_a) or (col in col_annot_b)

        rows.append({
            "Column": col + 1,
            f"{obj_a}_Residue": a_char,
            f"{obj_a}_ResNum": a_resnum,
            f"{obj_a}_ContactSide": a_side,
            f"{obj_b}_Residue": b_char,
            f"{obj_b}_ResNum": b_resnum,
            f"{obj_b}_ContactSide": b_side,
            "Conserved": "Y" if conserved else "",
            "ContactInEitherStructure": "Y" if either_contact else "",
        })

    with open(out_path, "w", newline="") as f:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    n_contact_cols = sum(1 for r in rows if r["ContactInEitherStructure"] == "Y")
    n_shared_contact_cols = sum(
        1 for r in rows if r[f"{obj_a}_ContactSide"] and r[f"{obj_b}_ContactSide"]
    )
    print(f"[{label}] alignment length: {n_cols} columns")
    print(f"[{label}] columns with a contact residue in either structure: {n_contact_cols}")
    print(f"[{label}] columns where BOTH structures contact at that position: {n_shared_contact_cols}")
    print(f"[{label}] wrote: {out_path}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    seqs_2bfx = load_sequences(os.path.join(DATA_DIR, "2bfx_box.csv"))
    seqs_4af3 = load_sequences(os.path.join(DATA_DIR, "4af3_incenp.csv"))

    contacts_2bfx = load_contacts_by_chain(os.path.join(DATA_DIR, "2bfx_box_binding_residues.csv"))
    contacts_4af3 = load_contacts_by_chain(os.path.join(DATA_DIR, "4af3_incenp_box_binding_residues.csv"))

    # --- Aurora Kinase: 2bfx chain B  vs  4af3 chain A ---
    build_annotated_alignment(
        label="Aurora Kinase (2bfx_B vs 4af3_A)",
        obj_a="2bfx_B", seq_a=seqs_2bfx["2bfx_box_B"], contacts_a=contacts_2bfx.get("B", []),
        obj_b="4af3_A", seq_b=seqs_4af3["4af3_incenp_box_A"], contacts_b=contacts_4af3.get("A", []),
        out_path=os.path.join(OUT_DIR, "aurora_kinase_alignment_annotated.csv"),
    )

    # --- INCENP IN-box: 2bfx chain C  vs  4af3 chain D ---
    build_annotated_alignment(
        label="INCENP IN-box (2bfx_C vs 4af3_D)",
        obj_a="2bfx_C", seq_a=seqs_2bfx["2bfx_box_C"], contacts_a=contacts_2bfx.get("C", []),
        obj_b="4af3_D", seq_b=seqs_4af3["4af3_incenp_box_D"], contacts_b=contacts_4af3.get("D", []),
        out_path=os.path.join(OUT_DIR, "incenp_inbox_alignment_annotated.csv"),
    )


if __name__ == "__main__":
    main()