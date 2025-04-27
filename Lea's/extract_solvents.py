import pandas as pd
import os

def extract_solvents_from_smiles(file_path, smiles_column='rxn_Smiles'):
    # Chargement du fichier CSV
    df = pd.read_csv(file_path)
    
    # Liste des solvants communs et leurs SMILES
    known_solvents = [
        'CO',      # Méthanol (CH3OH)
        'CCO',     # Éthanol (C2H5OH)
        'CC(C)=O', # Acétone (C3H6O)
        'ClCH2Cl', # Dichlorométhane (CH2Cl2)
        'C1CCO1',  # Tétrahydrofuran (THF) (C4H8O)
        'CC#N',    # Acétonitrile (CH3CN)
        'CS[O]',   # DMSO (Diméthylsulfoxyde, C2H6OS)
        'CN(C)C=O',# DMF (Diméthylformamide, C3H7NO)
        'Cc1ccccc1', # Toluène (C6H5CH3)
        'ClC(Cl)Cl', # Chloroforme (CHCl3)
        'C1CCCCC1',  # Cyclohexane (C6H12)
        'CCCCCC',    # Hexane (C6H14)
        'C1=CC=CC=C1', # Benzène (C6H6)
        'CCCCCCO',   # Hexanol (C6H14O)
        'CC(=O)O'    # Acide acétique (C2H4O2)
    ]
    
    # Fonction pour extraire les solvants des réactifs dans la réaction SMILES
    def get_solvents(rxn_smiles):
        # Vérifier la partie avant '>>' qui correspond aux réactifs
        reactants = rxn_smiles.split('>>')[0]
        
        solvents_found = []  # Liste pour stocker les solvants trouvés
        
        # Vérifier si un solvant connu est présent parmi les réactifs
        for solvent in known_solvents:
            if solvent in reactants:
                solvents_found.append(solvent)  # Ajouter le solvant à la liste
        
        # Retourner une chaîne de solvants séparés par des virgules, ou None si aucun solvant n'est trouvé
        return ', '.join(solvents_found) if solvents_found else None
    
    # Appliquer la fonction pour extraire les solvants
    df['solvents'] = df[smiles_column].apply(get_solvents)
    
    # Créer le chemin pour le fichier modifié avec un suffixe "_with_solvents"
    base_name, ext = os.path.splitext(file_path)
    modified_file_path = f"{base_name}_with_solvents{ext}"
    
    # Sauvegarder le fichier modifié avec le nouveau nom
    df.to_csv(modified_file_path, index=False)
    
    print(f"Le fichier a été modifié et sauvegardé sous : {modified_file_path}")

# Exemple d'utilisation :
file_path = "/Users/lealombard/PPChem_Project/Lea's/extract.csv"  # Remplace par le chemin réel de ton fichier
extract_solvents_from_smiles(file_path)
