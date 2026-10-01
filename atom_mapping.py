from rdkit import Chem
from rxnmapper import RXNMapper


def mol_to_tracked_smiles(mol: Chem.Mol) -> str:
    # Hydrogens must already be explicit here.  Keeping this assertion close
    # to reaction-SMILES construction prevents silently sending an unmapped H
    # population to RXNMapper.
    if any(atom.GetNumImplicitHs() for atom in mol.GetAtoms()):
        raise ValueError(
            "Reaction SMILES requires explicit hydrogens; call "
            "add_explicit_hydrogens() before make_mapping()"
        )
    tracked = Chem.Mol(mol)
    for atom in tracked.GetAtoms():
        atom.SetIsotope(atom.GetIdx() + 1)
    return Chem.MolToSmiles(tracked, canonical=False, isomericSmiles=True)


def parse_mapped_reaction(mapped_rxn: str):
    reactant_smiles, product_smiles = mapped_rxn.split(">>")
    reactants = [Chem.MolFromSmiles(s) for s in reactant_smiles.split(".") if s]
    products = [Chem.MolFromSmiles(s) for s in product_smiles.split(".") if s]
    if any(mol is None for mol in reactants + products):
        raise RuntimeError("Failed to parse RXNMapper output")
    return reactants, products


def _tracked_maps(mols):
    result = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            rxn_map = atom.GetAtomMapNum()
            if not rxn_map:
                continue
            isotope = atom.GetIsotope()
            if not isotope:
                raise RuntimeError("RXNMapper output lost original atom isotope")
            if rxn_map in result:
                raise RuntimeError(f"Duplicate RXN atom map: {rxn_map}")
            result[rxn_map] = isotope - 1
    return result


def make_mapping(reactant_mol: Chem.Mol, product_mol: Chem.Mol):
    reaction = f"{mol_to_tracked_smiles(reactant_mol)}>>{mol_to_tracked_smiles(product_mol)}"
    print("\nTracked reaction SMILES:")
    print(reaction, "\n")
    print("Running RXNMapper...")
    result = RXNMapper().get_attention_guided_atom_maps(
        [reaction], canonicalize_rxns=False
    )[0]
    mapped_rxn = result["mapped_rxn"]
    confidence = result["confidence"]
    print("Mapped reaction:")
    print(mapped_rxn, "\n")
    print(f"RXNMapper confidence: {confidence:.6f}\n")
    mapped_reactants, mapped_products = parse_mapped_reaction(mapped_rxn)
    reactant_maps = _tracked_maps(mapped_reactants)
    product_maps = _tracked_maps(mapped_products)
    mapping = {
        r_idx: product_maps[rxn_map]
        for rxn_map, r_idx in reactant_maps.items()
        if rxn_map in product_maps
    }
    return mapping, confidence, mapped_rxn


def validate_mapping(reactant_mol, product_mol, mapping):
    for r_idx, p_idx in mapping.items():
        if reactant_mol.GetAtomWithIdx(r_idx).GetSymbol() != product_mol.GetAtomWithIdx(p_idx).GetSymbol():
            raise RuntimeError(f"Element mismatch: reactant {r_idx}, product {p_idx}")
    print("Element consistency check: OK")
