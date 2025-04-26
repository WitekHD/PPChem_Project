import pandas as pd

from aizynthfinder.aizynthfinder import AiZynthFinder
from aizynthfinder.context.config import Configuration




def retrosynthesis_reaction_smiles(smiles: str, config_path: str = "config.yml") -> pd.DataFrame:
    """
    Perform retrosynthesis and return a table of forward-ordered one-step Reaction SMILES.

    Args:
        smiles (str): Target molecule in SMILES format.
        config_path (str): Path to AiZynthFinder's config.yml file.

    Returns:
        pd.DataFrame: Table with step number, reactants, product, and Reaction SMILES.
    """
    # Load AiZynthFinder
    config = Configuration(config_path)
    finder = AiZynthFinder()
    finder.config = config
    finder.target_smiles = smiles
    finder.tree_search()
    
    # Run retrosynthesis
    
    if not finder.routes:
        print("No routes found.")
        return []

    # Extract best route
    route = finder.routes.best

    if not route:
        raise ValueError("No synthesis route found for the given molecule.")

    # Collect reactions (retrosynthesis order)
    reactions = []
    for node in route.nodes:
        if node.is_reaction_node:
            reactants = ".".join(mol.smiles for mol in node.reactants)
            product = node.smiles
            reaction_smiles = f"{reactants}>>{product}"
            reactions.append((reactants, product, reaction_smiles))

    # Reverse to forward synthesis order
    reactions = list(reversed(reactions))

    # Create a DataFrame
    df = pd.DataFrame(
        [(i + 1, r[0], r[1], r[2]) for i, r in enumerate(reactions)],
        columns=["Step", "Reactants (SMILES)", "Product (SMILES)", "Reaction SMILES"]
    )

    return df


