import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

# Titre de l'application
st.title("Green solvent")

# Demander à l'utilisateur d'entrer le SMILES du produit
smiles_input = st.text_input("Enter the SMILES of the reaction product", "")

if smiles_input:
    # Convertir le SMILES en molécule
    mol = Chem.MolFromSmiles(smiles_input)
    
#    if mol:
#        # Dessiner la molécule
#        img = Draw.MolToImage(mol)
#        st.image(img, caption="Chemical structure of the product")
        
        # Afficher des informations supplémentaires sur la molécule
#        st.write("Molecule Information:")
        
        # Exemple : afficher la masse molaire
#        mol_weight = Chem.Descriptors.MolWt(mol)
#        st.write(f"Molar weight : {mol_weight:.2f} g/mol")
        
#    else:
#        st.write("The entered SMILES is not valid.")
