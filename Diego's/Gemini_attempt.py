import os
from dotenv import load_dotenv
from google import genai
from google.genai import types
load_dotenv()

rxn_name="Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters" #Change 
solvent_list= [
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
        'CC(=O)O',    # Acide acétique (C2H4O2)
        'CCCCO', # Butan-1-ol (C4H10O)
        'O', # Water
    ]

client = genai.Client(api_key=os.environ.get("GEMINI_API_KEY")) 
# API key is supposed to be in the .env and .env has to be in gitignore
prompt=f""" 
1. Main goal and context: 
You are an expert in assigning solvents to reactions and know which solvent can be used for a given reactant. 
You are part of a retro-synthesis program which will be used in mostly in solvent prediction (this is your job). 
After your answer, the rest of the code will evaluate the “greeness” (how sustainable the possible solvents are).
I hence need you to output only the Simplified Molecular Input Line Entry System (SMILES) of up to three possible 
solvents which can be used for the given reaction name (or class if the name cannot be extracted for you, 
of course, the class is quite general therefore you can be more general for your answer too). 

The solvent you must find is to do with the following reaction name/type: {rxn_name}

2. Constraints and examples
The solvent you propose must be part of this list: {solvent_list}

If you cannot find two solvents, one will do. 
You must output at least one solvent and everything you output must be in the list 
and in smiles format!
You MUST ONLY output the smiles of the solvents in the following format: 
“solventsmiles_1, solventsmiles_2, solventsmiles_3”

This is an example output for you to visualise with the SMILES: "CN(C)C=O, ClCCl, CS(C)=O"

3. Problems
If you are given a reactionn name or type which you do now know how to answer, you MUST simply reply with "nan"

"""
#use f string
#tell him exactly what to do
#tell him what not to do
#fix errors - ask for a specific format to allow the usage of the answer

response = client.models.generate_content(model="gemini-2.0-flash", contents=[prompt])
response_stripped=response.text.strip()
print(response_stripped)