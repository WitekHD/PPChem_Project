import pandas as pd
from striprtf.striprtf import rtf_to_text
from io import StringIO
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import re
import os

# Fonction pour nettoyer les SMILES (supprime les annotations de type :xx)
def clean_smiles(smiles):
    return re.sub(r':\d+', '', smiles)

# Chemins des fichiers
rtf_path = "/Users/lealombard/PPChem_Project/Raw_ReactionsA.rtf"
csv_output = "/Users/lealombard/PPChem_Project/Reactions_extraites.csv"
txt_output = "/Users/lealombard/PPChem_Project/Reactions_extraites.txt"
img_dir = "/Users/lealombard/PPChem_Project/Reactions_Images"

# Créer le dossier pour les images s'il n'existe pas
os.makedirs(img_dir, exist_ok=True)

# Lire le contenu du fichier .rtf
with open(rtf_path, 'r', encoding='utf-8') as file:
    rtf_content = file.read()

# Convertir en texte brut
plain_text = rtf_to_text(rtf_content)

# Parser le texte en DataFrame
data_io = StringIO(plain_text)
df = pd.read_csv(data_io)

# Sauvegarder en CSV
df.to_csv(csv_output, index=False)
print(f"✅ Fichier CSV sauvegardé : {csv_output}")

print("🔍 Colonnes détectées dans le fichier :")
print(df.columns.tolist())


# Sauvegarder les réactions dans un fichier texte
with open(txt_output, "w", encoding="utf-8") as f:
    for idx, row in df.iterrows():
        f.write(f"--- Réaction {idx + 1} ---\n")
        f.write(f"ID du brevet       : {row['patentID']}\n")
        f.write(f"Classe de réaction : {row['yrxn_Class']}\n")
        f.write(f"SMILES réaction    :\n{row['rxn_Smiles']}\n")
        f.write(f"Reactant Set       : {row['reactantSet']}\n\n")
print(f"✅ Fichier texte lisible sauvegardé : {txt_output}")

# Générer les images de réactions chimiques
for idx, row in df.iterrows():
    raw_smiles = row['rxn_Smiles']
    cleaned_smiles = clean_smiles(raw_smiles)

    try:
        rxn = AllChem.ReactionFromSmarts(cleaned_smiles, useSmiles=True)
        if rxn is not None:
            img = Draw.ReactionToImage(rxn, subImgSize=(300, 300))
            img_path = os.path.join(img_dir, f"reaction_{idx + 1}.png")
            img.save(img_path)
            print(f"✅ Image de la réaction {idx + 1} sauvegardée : {img_path}")
        else:
            print(f"❌ Impossible d'interpréter la réaction {idx + 1}")
    except Exception as e:
        print(f"❌ Erreur pour la réaction {idx + 1} : {raw_smiles}")
        print(f"    → Détail : {e}")