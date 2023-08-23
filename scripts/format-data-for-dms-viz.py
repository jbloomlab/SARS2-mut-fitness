"""Parse the mutation data for dms-viz."""

import os
import Bio.PDB
from Bio.SeqUtils import seq1
import pandas as pd
import requests
import warnings
from io import StringIO

# Functions 
def get_structure(pdb_input):
    """
    Fetch a PDB structure from the RCSB PDB web service or load it from a local file.

    This function takes a string as input, which should either be a 4-character PDB ID or
    a path to a local PDB file. The function fetches the structure with the specified PDB ID
    from the RCSB PDB web service, or reads the structure from the specified local PDB file,
    and returns a Bio.PDB structure object.

    Parameters
    ----------
    pdb_input : str
        A string that is either a 4-character PDB ID or a path to a local .pdb file.

    Returns
    -------
    structure : Bio.PDB.Structure.Structure
        A Bio.PDB structure object.

    Raises
    ------
    ValueError
        If the pdb_input is neither a valid PDB ID nor a local PDB file path.
        If there was an error reading the local PDB file or parsing the PDB content.
        If there was an error downloading the PDB file from the RCSB PDB web service.

    """

    # Check if the input is a local file path
    if os.path.isfile(pdb_input) and pdb_input.endswith(".pdb"):
        try:
            # Ignore warnings about discontinuous chains
            with warnings.catch_warnings():
                warnings.simplefilter(
                    "ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning
                )
                structure = Bio.PDB.PDBParser().get_structure(pdb_input[:-4], pdb_input)
        except Exception as e:
            raise ValueError(f"Error reading PDB file {pdb_input}: {e}")
    elif len(pdb_input) == 4 and pdb_input.isalnum():  # Check for a valid PDB ID format
        # Try to fetch the structure from RCSB PDB
        response = requests.get(f"https://files.rcsb.org/download/{pdb_input}.cif")
        if response.status_code == 200:
            try:
                pdb_file_content = StringIO(response.text)
                # Ignore warnings about discontinuous chains
                with warnings.catch_warnings():
                    warnings.simplefilter(
                        "ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning
                    )
                    structure = Bio.PDB.MMCIFParser().get_structure(
                        pdb_input, pdb_file_content
                    )
            except Exception as e:
                raise ValueError(f"Error parsing PDB content for {pdb_input}: {e}")
        else:
            raise ValueError(
                f"Failed to download {pdb_input} from the RCSB database. Status code: {response.status_code}"
            )
    else:
        raise ValueError(
            f"Invalid input: {pdb_input}. Please provide a valid PDB ID or a local PDB file path."
        )

    return structure


def get_structure_sequence(structure):
    """Get the sequence and index for chains in a given structure"""
        
    # Get the sequence for each chain (assuming a single model)
    structure_seqs = {}
    for model in structure:
        for chain in model:
            sequence = ""
            indices = []

            for residue in chain:
                if Bio.PDB.is_aa(residue):  
                    sequence += seq1(residue.get_resname())
                    indices.append(residue.id[1])  
                    
            structure_seqs[chain.id] = [sequence, indices]
            
    return structure_seqs


# Inputs
fitness_effects = snakemake.input.aa_fitness
wildtype_residues = snakemake.input.clade_founder_aas
structure_info = snakemake.input.structure_info
protein = snakemake.wildcards.protein
# Outputs
fitness_output = snakemake.output.fitness_df
sitemap_output = snakemake.output.sitemap_df

# Read in the fitness effects for all genes
all_gene_fitness_df = (
    pd.read_csv(fitness_effects)
    .rename(columns={'aa_site': 'site', 'aa': 'mutant'})
    .drop(columns=['aa_differs_among_clade_founders', 'subset_of_ORF1ab'])
)

# Read in the 'Wildtype' amino acids from clade 19A ancestor
wildtype_aa_df = (
    pd.read_csv(wildtype_residues)
        .query("clade == '19A (B)'")[['gene', 'site', 'amino acid']]
        .rename(columns={'amino acid': 'wildtype'})
        .drop_duplicates()
        .sort_values(['gene', 'site'])
        .reset_index(drop=True)
)
all_gene_fitness_df['gene'] = all_gene_fitness_df['gene'].str.split(' ').str[0]
wildtype_aa_df['gene'] = wildtype_aa_df['gene'].str.split(' ').str[0]

# Make the 'fitness' dataframe for the gene
fitness_df = (all_gene_fitness_df
                    .merge(wildtype_aa_df, how='left', on=['gene', 'site'])
                    .query(f"gene == '{protein}'")
                    .drop(columns=['gene'])
                    )

# Make the 'sitemap' dataframe for the gene
sitemap = (fitness_df[['site']]
            .drop_duplicates()
            .reset_index(drop=True)
            .rename(columns={'site': 'reference_site'})
            .sort_values(['reference_site'])
            )
sitemap['sequential_site'] = sitemap.index + 1                 

# Some of the structures are not aligned to the reference sequence
structure_info = pd.read_csv(structure_info)
pdb_file = structure_info.query(f"selection == '{protein}'")['pdb'].values[0]
offset = structure_info.query(f"selection == '{protein}'")['offset'].values[0]
indexChain = structure_info.query(f"selection == '{protein}'")['indexChain'].values[0]

if not pd.isnull(offset):
    print(f"Offset for {protein} and chains {indexChain}: {offset}")
    # Load the structure into a BioPython object
    structure = get_structure(pdb_file)
    # Get the sequence and index for each chain in the PDB
    structure_sequences = get_structure_sequence(structure)
    # Get the current wildtype residue at each site
    protein_df = fitness_df[['site', 'wildtype']].drop_duplicates().reset_index(drop=True)

    # Only get the sequence for the included chains
    for chain, values in structure_sequences.items(): 
        if indexChain == "None" or chain == indexChain:
            indices = values[1]
            sequence = values[0]
            
    # Get the indicies and sequence for the data and structure
    structure_seq = dict(zip([i for i in indices], sequence))
    data_seq = dict(zip(protein_df['site'], protein_df['wildtype']))
    
    protein_sites = []
    after_offset = []
    # Use this index by updating the protein_site column of the sitemap dataframe
    for site, res in data_seq.items():
        index = site + offset
        if index in structure_seq.keys():
            after_offset.append(res == structure_seq[index])
            # Add to protein_site column 
            protein_sites.append(index)
        else:
            # Add a nonsense value
            protein_sites.append(-1000)
            
    print(f"{sum(after_offset) / len(after_offset) * 100}% match after offset.")

    sitemap['protein_site'] = protein_sites
 
# Write out the sitemap and fitness dataframes
sitemap.to_csv(sitemap_output, index=False)
fitness_df.to_csv(fitness_output, index=False)
