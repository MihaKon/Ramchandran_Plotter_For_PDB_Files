import pathlib
import argparse

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

from Bio.PDB import PDBParser, PPBuilder, calc_dihedral, Structure, Model, Residue


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path", type=pathlib.Path)
    return parser.parse_args()


def read_pdb(file_path: pathlib.Path, file_name: str) -> Structure.Structure:
    parser = PDBParser()
    structure = parser.get_structure(file_name, file_path)
    return structure


def calculate_phi_psi(model: Model.Model) -> None:
    ppb = PPBuilder()
    phi_psi_angles = []

    for polypeptide in ppb.build_peptides(model):
        residues = polypeptide
        for i in range(1, len(residues) - 1):
            try:
                prev_res: Residue.Residue = residues[i - 1]
                curr_res: Residue.Residue = residues[i]
                next_res: Residue.Residue = residues[i + 1]

                n = curr_res["N"].get_vector()
                ca = curr_res["CA"].get_vector()
                c = curr_res["C"].get_vector()

                prev_c = prev_res["C"].get_vector()

                next_n = next_res["N"].get_vector()

                phi = calc_dihedral(prev_c, n, ca, c)
                psi = calc_dihedral(n, ca, c, next_n)

                phi_psi_angles.append((np.degrees(phi), np.degrees(psi)))
            except KeyError:
                phi_psi_angles.append((None, None))

    return phi_psi_angles


def visualize_ramachandran_plot(
    file_name: str, model_id: str, angles: list[tuple[np.float32, np.float32]]
) -> plt:
    angles_df = pd.DataFrame(angles, columns=["Phi", "Psi"])
    plt.figure(figsize=(10, 8))
    sns.kdeplot(
        data=angles_df,
        x="Phi",
        y="Psi",
        cmap="coolwarm",
        fill=True,
    )
    sns.scatterplot(data=angles_df, x="Phi", y="Psi", color="black", alpha=0.3, s=10)
    plt.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.title(f"Ramachandran Plot: {file_name} | {model_id}", fontsize=16)
    plt.xlabel("Phi (°)", fontsize=12)
    plt.ylabel("Psi (°)", fontsize=12)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(visible=True, linestyle="--", alpha=0.6)
    return plt


def main():
    args = parse_args()
    file_name = args.file_path.name.split(".")[0].upper()
    structure = read_pdb(args.file_path, file_name)

    models = list(structure.get_models())
    if not models:
        print("No models found in the PDB file.")
        return

    ramachandran_plots = {}
    for model in models:
        phi_psi_angles = calculate_phi_psi(model)
        ramachandran_plots[model.id] = phi_psi_angles

    for model_id, angles in ramachandran_plots.items():
        plt = visualize_ramachandran_plot(file_name, model_id, angles)
        plt.savefig(f"{file_name}_model_{model_id}_ramachandran_plot.png")


if __name__ == "__main__":
    main()
