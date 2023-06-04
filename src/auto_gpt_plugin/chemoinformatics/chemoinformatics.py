import os
import subprocess
from pathlib import Path
from typing import Any, List, Tuple, Union

from rdkit.Chem import AllChem as Chem
from rdkit.DataStructs import TanimotoSimilarity


def calculate_similarity(smiles1: str, smiles2: str) -> float:
    """
    Calculate the similarity between two chemical structures given their SMILES strings using RDKit.

    This function computes the Tanimoto coefficient between the Morgan fingerprints of the input SMILES strings.

    Args:
        smiles1 (str): The SMILES string of the first chemical structure.
        smiles2 (str): The SMILES string of the second chemical structure.

    Returns:
        float: The Tanimoto coefficient between the input chemical structures, representing their similarity.
    """
    # SMILES文字列からRDKitの分子オブジェクトを作成
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # 分子オブジェクトからMorganフィンガープリントを計算
    fp1 = Chem.GetMorganFingerprintAsBitVect(mol1, 2)
    fp2 = Chem.GetMorganFingerprintAsBitVect(mol2, 2)

    # タニモト係数を計算して返す
    similarity = TanimotoSimilarity(fp1, fp2)
    return similarity


def predict_CBL_B_activity(smiles: str) -> float:
    """
    Predict CBL_B activity of chemical structure fiven its SMILES string using machine learning model.

    Args:
        smiles (str): The SMILES string of the chemical structure.

    Returns:
        float: Predicted activity value.
    """
    # TODO: 予測モデルをビルドして予測する
    return 0.0


def validate_smiles_string(
    smiles: Union[List[str], str]
) -> Union[List[bool], bool]:
    """
    Evaluate whether SMILES strings are valid and returns True if so

    Args:
        smiles (List[str] or str): The SMILES string(s) of the chemical structure.

    Returns:
        List[bool] or bool: is valid or not
    """

    def _validate(smiles: str):
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            return False
        return mol is not None

    if type(smiles) == list:
        return [_validate(s) for s in smiles]
    elif type(smiles) == str:
        return _validate(smiles)
    else:
        return False


def scaffold_hopping(smiles: str) -> str:
    pass


def target_protein_preparation(setting_json: str) -> Tuple[str]:
    """
    Do target protein preparation with preparation setting json

    Args:
        setting_json (str): Setting json absolute path for target protein preparation

    Return:
        receptor_path, logfile
    """
    import json
    import os
    from pathlib import Path

    dockstream_path = os.environ["DOCKSTREAM_PATH"]
    dockstream_env = os.environ["DOCKSTREAM_ENV"]

    target_preparator = Path(dockstream_path) / "target_preparator.py"
    python_path = Path(dockstream_env) / "bin" / "python3"

    workdir = os.environ["WORKDIR"]
    os.chdir(workdir)

    command = [
        str(python_path),
        str(target_preparator),
        "-conf",
        setting_json,
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error: {stderr.decode()}")
    else:
        print(f"Output: {stdout.decode()}")

    with open(setting_json, "r") as f:
        setting = json.load(f)

    receptor_path = setting["target_preparation"]["runs"][0]["output"][
        "receptor_path"
    ]
    log_file = setting["target_preparation"]["header"]["logging"]["logfile"]

    return receptor_path, log_file


# def dock_ligand_compounds_and_target_protein(json_setting: str) -> None:
#     """
#     Dock ligand compounds and target protein with setting json
#     """
#     subprocess


def make_apo_protein_pdb(pdb: str) -> str:
    """
    Make apo protein (only protein, not including ligand compounds) PDB

    Args:
        pdb str: Input PDB file path

    Returns:
        str: Apo PDB file path
    """
    import os
    from Bio.PDB import PDBIO, PDBParser, Select
    from pathlib import Path

    class SelectProtein(Select):
        def accept_residue(self, residue):
            if (
                residue.get_full_id()[3][0] == " "
                or residue.get_resname() == "HOH"
            ):
                return 1
            return 0

    workdir = os.environ["WORKDIR"]
    os.chdir(workdir)

    name = Path(pdb).stem

    parser = PDBParser()
    structure = parser.get_structure("protein_structure", pdb)

    io = PDBIO()
    io.set_structure(structure)

    output_name = name + "_apo.pdb"
    io.save(output_name, select=SelectProtein())
    return output_name


def make_only_ligand_compound_pdb(pdb: str) -> str:
    """
    Make only ligand compounds PDB

    Args:
        pdb str: Input PDB file path

    Returns:
        str: only ligand compound PDB file path
    """
    import os
    from Bio.PDB import PDBIO, PDBParser, Select
    from pathlib import Path

    class SelectLigand(Select):
        def __init__(self, ligand_resname):
            self.ligand_resname = ligand_resname

        def accept_residue(self, residue):
            if (
                residue.get_full_id()[3][0] == " "
                or residue.get_resname() == "HOH"
                or residue.get_resname() != self.ligand_resname
            ):
                return 0
            return 1

    workdir = os.environ["WORKDIR"]
    os.chdir(workdir)

    parser = PDBParser()
    structure = parser.get_structure("protein_structure", pdb)

    io = PDBIO()
    io.set_structure(structure)

    ligand_resnames = set(
        residue.get_resname()
        for residue in structure.get_residues()
        if residue.get_full_id()[3][0] != " "
    )
    ligand_resnames.discard("HOH")

    for ligand_resname in ligand_resnames:
        io.save(
            ligand_resname + ".pdb",
            SelectLigand(ligand_resname),
        )
    return [ligand_resname + ".pdb" for ligand_resname in ligand_resnames]


def run_reinvent_with_docking(setting_json: str) -> None:
    """
    Run reinvnet with docking by DockStream

    Args:
        setting_json (str): Setting json path for reinvent

    Return:
        None
    """
    import json
    import os
    from pathlib import Path

    reinvent_path = os.environ["REINVENT_PATH"]
    reinvent_env = os.environ["REINVENT_ENV"]

    input_ = Path(reinvent_path) / "input.py"
    python_path = Path(reinvent_env) / "bin" / "python"

    workdir = os.environ["WORKDIR"]
    os.chdir(workdir)

    command = [
        str(python_path),
        str(input_),
        setting_json,
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error: {stderr.decode()}")
    else:
        print(f"Output: {stdout.decode()}")

    return


def compare_smiles_with_sdf(csv_file: str, sdf_file: str, output_file: str):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit import DataStructs

    """
    Compare smiles in text file with SDF and output smiles which not included
    in SDF and dissimilarity top 10

    Args:
        csv_file (str): smiles csv
        sdf_file (str): sdf file
        output_file (str): output file name for result csv
    
    Return:
        output_file name (str)
    """

    def read_sdf_to_mols(file_path):
        suppl = Chem.SDMolSupplier(file_path)
        mols = [mol for mol in suppl if mol is not None]
        return mols

    def read_csv_to_mols(file_path):
        df = pd.read_csv(file_path)
        mols = [Chem.MolFromSmiles(smiles) for smiles in df["SMILES"]]
        return mols

    def compute_avg_tanimoto(csv_mols, sdf_mols):
        fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator()
        csv_fps = [fp_gen.GetFingerprint(mol) for mol in csv_mols]
        sdf_fps = [fp_gen.GetFingerprint(mol) for mol in sdf_mols]

        avg_tanimoto = []
        for csv_fp in csv_fps:
            tanimoto_scores = [
                DataStructs.FingerprintSimilarity(csv_fp, sdf_fp)
                for sdf_fp in sdf_fps
            ]
            avg_tanimoto.append(sum(tanimoto_scores) / len(tanimoto_scores))

        return avg_tanimoto

    sdf_mols = read_sdf_to_mols(sdf_file)
    csv_mols = read_csv_to_mols(csv_file)

    # Convert sdf_mols to a set for faster lookup
    sdf_mols_set = set(Chem.MolToSmiles(mol) for mol in sdf_mols)
    not_included_mols = [
        mol for mol in csv_mols if Chem.MolToSmiles(mol) not in sdf_mols_set
    ]

    not_included_smiles = [Chem.MolToSmiles(mol) for mol in not_included_mols]

    avg_tanimoto = compute_avg_tanimoto(not_included_mols, sdf_mols)
    df_not_included = pd.DataFrame(
        {
            "SMILES": not_included_smiles[:10],
            "Average Tanimoto Similarity": avg_tanimoto[:10],
        }
    )
    df_not_included.to_csv(output_file, index=False)
    return output_file
