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


def dock_compound_to_CBL_B(smiles: str) -> float:
    pass
