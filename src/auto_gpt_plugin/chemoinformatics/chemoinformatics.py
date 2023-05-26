import os
import subprocess
from pathlib import Path
from typing import Any, List, Union

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


def generate_molecule() -> List[Any]:
    pass


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


# def generate_target_protein_preparation_setting_json() -> None:
#     """
#     Generate setting json for preparation of target protein before docking. An example of a setting file is shown below.

#     ```json
#     {
#         "target_preparation": {
#             "header": {                                                                                                                     # general settings
#                 "logging": {                                                                                                                # logging settings (e.g. which file to write to)
#                     "logfile": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_target_prep.log"
#                 }
#             },
#             "input_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/DockStreamCommunity/data/1UYD/1UYD_apo.pdb",                        # this should be an absolute path
#             "fixer": {
#                 "enabled": true,
#                 "standardize": True,                                                                                                        # enables standardization of residues
#                 "remove_heterogens": True,                                                                                                  # remove hetero-entries
#                 "fix_missing_heavy_atoms": True,                                                                                            # if possible, fix missing heavy atoms
#                 "fix_missing_hydrogens": True,                                                                                              # add hydrogens, which are usually not present in PDB files
#                 "fix_missing_loops": False,                                                                                                 # add missing loops; CAUTION: the result is usually not sufficient
#                 "add_water_box": False,                                                                                                     # if you want to put the receptor into a box of water molecules
#                 "fixed_pdb_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_fixed_target.pdb"                              # if specified and not "None", the fixed PDB file will be stored here
#             },
#             "runs": [                                                                                                                       # "runs" holds a list of backend runs; at least one is required
#                 {
#                     "backend": "AutoDockVina",                                                                                              # Use AutoDockVina
#                     "output": {
#                         "receptor_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_receptor.pdbqt"                         # the generated receptor file will be saved to this location
#                     },
#                     "parameters": {
#                         "pH": 7.4,                                                                                                          # sets the protonation states (NOT used in Vina)
#                         "extract_box": {                                                                                                    # in order to extract the coordinates of the pocket (see text)
#                             "reference_ligand_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/DockStreamCommunity/data/1UYD/PU8.pdb",  # path to the reference ligand
#                             "reference_ligand_format": "PDB"                                                                                # format of the reference ligand
#                         }
#                     }
#                 }
#             ]
#         }
#     }
#     ```
#     """


def target_protein_preparation(setting_json: str) -> None:
    """
    Do target protein preparation with preparation setting json

    Args:
        setting_json (str): Setting json absolute path for target protein preparation

    Return:
        None: Output file path is assigned in setting file
    """
    import os
    from pathlib import Path

    dockstream_path = os.environ["DOCKSTREAM_PATH"]
    dockstream_env = os.environ["DOCKSTREAM_ENV"]

    target_preparator = Path(dockstream_path) / "target_preparator.py"
    python_path = Path(dockstream_env) / "bin" / "python3"

    workdir = os.environ["WORKDIR"]
    setting_json = str((Path(workdir) / setting_json).resolve())

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


# def generate_docking_setting_json() -> None:
#     """
#     Generate setting json for docking target protein and ligands. An example of a setting file is shown below.

#     ```json
#     {
#       "docking": {
#         "header": {                                                                                                                         # general settings
#           "logging": {                                                                                                                      # logging settings (e.g. which file to write to)
#             "logfile": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_docking.log"
#           }
#         },
#         "ligand_preparation": {                                                                                                             # the ligand preparation part, defines how to build the pool
#           "embedding_pools": [
#             {
#               "pool_id": "RDkit",                                                                                                           # Use RDkit
#               "type": "RDkit",                                                                                                              # Use RDkit
#               "parameters": {
#                 "prefix_execution": ""                                                                                                      # only required, if a module needs to be loaded to execute "Corina"
#               },
#               "input": {
#                 "standardize_smiles": false,
#                 "type": "smi",
#                 "input_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/DockStreamCommunity/notebooks/../data/1UYD/ligands_smiles.txt"
#               },
#               "output": {                                                                                                                   # the conformers can be written to a file, but "output" is not required as the ligands are forwarded internally
#                 "conformer_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_embedded_ligands.sdf",
#                 "format": "sdf"
#               }
#             }
#           ]
#         },
#         "docking_runs": [
#           {
#             "backend": "AutoDockVina",
#             "run_id": "AutoDockVina",
#             "input_pools": [
#               "RDkit"
#             ],
#             "parameters": {
#               "binary_location": "/home/vscode/micromamba/envs/DockStream/bin/",  # fixed
#               "parallelization": {
#                 "number_cores": 4
#               },
#               "seed": 42,                                                                                                                   # use this "seed" to generate reproducible results; if varied, slightly different results will be produced paths to the receptor files number of poses to be generated search space (cavity definition); see text
#               "receptor_pdbqt_path": [
#                 "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_receptor.pdbqt"                                                  # The target protein to be assigned must be pretreated.
#               ],
#               "number_poses": 2,
#               "search_space": {
#                 "--center_x": 3.3,
#                 "--center_y": 11.5,
#                 "--center_z": 24.8,
#                 "--size_x": 15,
#                 "--size_y": 10,
#                 "--size_z": 10
#               }
#             },
#             "output": {
#               "poses": {
#                 "poses_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_ligands_docked.sdf"
#               },
#               "scores": {
#                 "scores_path": "/workspaces/Auto-GPT-Plugin-Chemoinformatics/tutorial/ADV_scores.csv"
#               }
#             }
#           }
#         ]
#       }
#     }
#     ```
#     """


# def dock_ligand_compounds_and_target_protein(json_setting: str) -> None:
#     """
#     Dock ligand compounds and target protein with setting json
#     """
#     subprocess


# def dock_compound_to_CBL_B(smiles: str) -> float:
#     """
#     Evaluate whether SMILES strings are valid and returns True if so

#     Args:
#         smiles (List[str] or str): The SMILES string(s) of the chemical structure.

#     Returns:
#         List[bool] or bool: is valid or not
#     """
