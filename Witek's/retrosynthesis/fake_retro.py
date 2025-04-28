import pandas as pd

def manual_retrosynthesis(target_smiles, reaction_db_path):
    reaction_df = pd.read_csv(reaction_db_path)
    results = []
    
    def search_reactions(product_smiles, step_prefix):
        # Find reactions producing the product_smiles
        matches = reaction_df[reaction_df['Reaction_SMILES'].str.split(">>").str[1] == product_smiles]
        if matches.empty:
            return  # No reactions produce this molecule
        
        for i, row in matches.iterrows():
            reaction_smiles = row['Reaction_SMILES']
            step_number = f"{step_prefix}"
            results.append({
                "Step": step_number,
                "Reaction_SMILES": reaction_smiles
            })

            # Get reactants from this reaction
            reactants = reaction_smiles.split(">>")[0].split(".")
            
            # Recursively search for reactions producing each reactant
            for j, reactant in enumerate(reactants, start=1):
                sub_step_prefix = f"{step_number}.{j}"
                search_reactions(reactant, sub_step_prefix)

    # Start the search from the target molecule
    search_reactions(target_smiles, "1")

    return pd.DataFrame(results)
