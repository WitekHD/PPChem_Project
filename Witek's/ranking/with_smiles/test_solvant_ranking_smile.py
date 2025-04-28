import pandas as pd
import numpy as np
from tabulate import tabulate
from solvant_ranking_smile import rank_similar_solvents


# Get the results
results = rank_similar_solvents('CCO')

# Print each table with a header
print("\n=== Target Solvent Properties ===")
# Convert dictionary to DataFrame for better display
target_props_df = pd.DataFrame([results['target_solvent_properties']]).T
print(tabulate(target_props_df, headers=['Value'], tablefmt='pretty'))

print("\n=== Ranked by Similarity ===")
print(tabulate(results['by_similarity'], headers='keys', tablefmt='pretty', floatfmt='.3f'))

print("\n=== Ranked by Environmental Impact ===")
print(tabulate(results['by_environment'], headers='keys', tablefmt='pretty', floatfmt='.2f'))

print("\n=== Ranked by Health Impact ===")
print(tabulate(results['by_health'], headers='keys', tablefmt='pretty', floatfmt='.2f'))

print("\n=== Ranked by Safety ===")
print(tabulate(results['by_safety'], headers='keys', tablefmt='pretty', floatfmt='.2f'))

print("\n=== Ranked by Overall Score ===")
print(tabulate(results['by_overall_ranking'], headers='keys', tablefmt='pretty', floatfmt='.2f'))

