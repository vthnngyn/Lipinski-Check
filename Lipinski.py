#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
from rdkit.Chem import Lipinski as lip
from rdkit.Chem import Descriptors as des
from rdkit.Chem import Crippen
from rdkit.Chem import Draw
import requests
import os

H_DONERS_LIMIT = 5
H_ACCEPTORS_LIMIT = 10
MOLECULAR_MASS_LIMIT = 500
CLOGP_LIMIT = 5
TOTAL_RULES = 4

#command for loading molecule from input
def load_molecule(user_input: str):
    mol_data = user_input.strip()

    #extracting from a weblink
    if mol_data.startswith('http'):
        print("Is a URL")
        try:
            response = requests.get(mol_data)
            response.raise_for_status()
            mol_data = response.text
        except Exception as e:
            print(f"Cannot read URL: {e}")
            return None

        mol = Chem.MolFromMolBlock(mol_data)
        if mol:
            print("Is a MOL/SDF")
            return mol
        else:
            print("URL is not a MOL/SDF")
            return None

    #extracting from local file
    elif os.path.exists(mol_data):
        print("Is a local file")
        try:
            with open(mol_data, 'r') as f:
                mol_block = f.read()
            mol = Chem.MolFromMolBlock(mol_block)
            if mol:
                print("Is a MOL/SDF")
                return mol
            else:
                print("Is not a MOL/SDF")
                return None
        except Exception as e:
            print(f"Cannot read local file: {e}")
            return None

    #extracting from rich text
    mol = Chem.MolFromMolBlock(mol_data)
    if mol:
        print("Is a MOL/SDF")
        return mol

    #extracting from SMILES
    mol = Chem.MolFromSmiles(mol_data)
    if mol:
        print("Is a SMILES")
        return mol

    print("Cannot read input")
    return None

def main():
    #prompt for input
    user_input = input("Input URL, local file path, or formatted text of MOL/SDF, or SMILES").strip()
    mol = load_molecule(user_input)

    #counter for Lipinski's rules
    passed = 0

    #checking each of Lipinski's rules and whether inputted molecule passes or not
    if mol:
        print("Loading molecule")

        h_donors = lip.NumHDonors(mol)
        if h_donors <= H_DONERS_LIMIT:
            print(f"Checking if {H_DONERS_LIMIT} or less H-donors... HAS {h_donors} H-DONORS")
            passed += 1
        else:
            print(f"Checking if {H_DONERS_LIMIT} or less H-donors... RULE VIOLATION: HAS {h_donors} H-DONORS")

        h_acceptors = lip.NumHAcceptors(mol)
        if h_acceptors <= H_ACCEPTORS_LIMIT:
            print(f"Checking if {H_ACCEPTORS_LIMIT} or less H-acceptors... HAS {h_acceptors} H-ACCEPTORS")
            passed += 1
        else:
            print(f"Checking if {H_ACCEPTORS_LIMIT} or less H-acceptors... RULE VIOLATION: HAS {h_acceptors} H-ACCEPTORS")

        mass = round(des.MolWt(mol),2)
        if mass <= MOLECULAR_MASS_LIMIT:
            print(f"Checking if molecular mass is {MOLECULAR_MASS_LIMIT} or less daltons... {mass} DALTONS")
            passed += 1
        else:
            print(f"Checking if molecular mass is {MOLECULAR_MASS_LIMIT} or less daltons... RULE VIOLATION: {mass} DALTONS")

        cLogP = round(Crippen.MolLogP(mol),2)
        if cLogP <= CLOGP_LIMIT:
            print(f"Checking if computational partition coefficient is {CLOGP_LIMIT} or less... IS {cLogP}")
            passed += 1
        else:
            print(f"Checking if computational partition coefficient is {CLOGP_LIMIT} or less... RULE VIOLATION: IS {cLogP}")

        print(f"{passed} out of {TOTAL_RULES}")
    else:
        print("Molecule could not be loaded")

    if passed < TOTAL_RULES:
        print("Failed Lipinski check")
    else:
        print("Passed Lipinski check")

    #draw molecule
    img = Draw.MolToImage(mol, size=(300, 300))
    img.show()

if __name__ == "__main__":
    main()