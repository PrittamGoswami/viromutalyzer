from Bio import SeqIO
from Bio import Phylo
import pandas as pd
import re
from .reference_based_mutation_detection import save_dict_to_pickle, find_mutations_wrt_subject_sequence

__all__ = ["detect_mutations_phylogeny_based"]

def _get_reconstructed_ancestral_sequences(ancestral_state_file_path):
    """
    Reconstructs ancestral node sequences from the ancestral state file.
    
    Parameters:
        - ancestral_state_file_path (str): The path to the ancestral state file.
    
    Returns:
        - dict: A dictionary containing the reconstructed ancestral sequences, where the keys are the ancestral node ids and the values are the reconstructed sequences.
    """

    # Look for the start line
    with open(ancestral_state_file_path, 'r') as file:
        for line_count, line in enumerate(file):
            if line.strip().startswith("Node"):
                break  # Return the line number where the data starts
    
    # Read the ancestral state file
    ancestral_states = pd.read_csv(ancestral_state_file_path, sep='\t', skiprows=line_count)

    # Concatenate the ancestral states for all the positions for a given ancestral node to reconstruct the sequence for that node
    ancestral_sequences = {}
    for node, group in ancestral_states.groupby('Node'):
        sequence = ''.join(group['State'])
        ancestral_sequences[node] = sequence
    
    return ancestral_sequences


def _find_ancestral_node(tree, node_name):
    """
    To find the nearest ancestral node to a particular node in the tree.

    Parameters:
        - tree: The phylogenetic tree to search in.
        - node_name: The name of the node for which we wish to find the nearest ancestor.

    Returns:
        - ancestral_node: The ancestral node if found, None otherwise.
    """
    # Search the target node in the tree
    target_node = next((clade for clade in tree.find_clades() if str(clade.name).split("/")[0] == node_name), None)
    if target_node is None:
        print(f"Node '{node_name}' not found.")
        return None

    # Traverse the tree to find the nearest ancestral node to the target node
    ancestral_node = None
    for clade in tree.find_clades():
        if target_node in clade:
            ancestral_node = clade
            break
    
    #return nearest ancestral node
    return ancestral_node

tree = Phylo.read("Sample_1.treefile", "newick")



def _get_mutation_data_for_ancestral_nodes(ancestral_sequences, msa_file_path, tree, outgroup_id):
    """ 
    To compute the mutations for the internal nodes of the tree.

    Args:
        ancestral_sequences (dict): Dictionary of reconstructed ancestral sequences.
        msa_file_path (str): The path to the mutiple sequence alignment file.
        tree (tree object): The phylogenetic tree.
        outgroup_id(str or None): The outgroup sequence ID 

    Returns:
         - denovo_mutations: A dictionary mapping each internal node to its denovo mutations.
         - ancestral_mutations: A dictionary mapping each internal node to its ancestral mutations.
         - ancestral_node_data: A dictionary mapping each internal node to its nearest ancestral node.

    """
    # Initialize the dictionaries to store the mutation data
    denovo_mutations = {}
    ancestral_mutations = {}
    ancestral_node_data = {}

    # Sort the ancestral nodes based on their node number so that older ancestral nodes comes first. This is done to detect the ancestral mutations that are inherited. 
    ancestral_nodes = list(ancestral_sequences.keys())
    ancestral_nodes = sorted(ancestral_nodes, key=lambda node: int(re.search(r'\d+', node).group()))
    
    # Sequence records in multifasta file
    records = list(SeqIO.parse(msa_file_path, "fasta"))

    # Iterate through each ancestral node
    for node in ancestral_nodes:
        ancestral_node = _find_ancestral_node(tree, node)
        # If it is the initial node, then the reference sequence is used as the ancestral sequence
        if node == "Node1":
            # If the outgroup_id is specified
            if outgroup_id:
                # Search for corresponding record in the multifasta file
                outgroup_record = next((record for record in records if record.id == str(outgroup_id)), None)
                # If it is not found raise and error
                if not outgroup_record:
                    raise ValueError(f"Outgroup sequence with ID '{outgroup_id}' not found in the input FASTA file.")
                else:
                    # Otherwise compute the mutations
                    outgroup_seq = str(outgroup_record.seq).upper()
                    denovo_mutations[node] = find_mutations_wrt_subject_sequence(outgroup_seq, ancestral_sequences[node])
                    ancestral_mutations[node] = []
                    ancestral_node_data["Node1"] = outgroup_seq    
            else:
                # If the outgroup_id is not specified then continue
                continue
        else:
            # If it is not the initial node, then we find the nearest ancestor and compute the mutations against it
            denovo_mutations[node] = find_mutations_wrt_subject_sequence(ancestral_sequences[str(ancestral_node).split("/")[0]], ancestral_sequences[node])
            ancestral_mutations[node] = (denovo_mutations[str(ancestral_node).split("/")[0]]
                                        + ancestral_mutations[str(ancestral_node).split("/")[0]]
                                        )
            ancestral_node_data[node] = ancestral_node

    return denovo_mutations, ancestral_mutations, ancestral_node_data



def _get_mutation_data_for_leaf_nodes(input_msa_fasta_path, ancestral_sequences, denovo_mutations, ancestral_mutations, tree, ancestral_node_data, outgroup_id):
    """
    To compute the mutations for the leaf nodes of the tree
    
    Parameters:
        - input_msa_fasta_path (str): Path to the multiple sequence alignment file .
        - ancestral_sequences (dict): Dictionary of reconstructed ancestral sequences.
        - denovo_mutations (dict): Dictionary to store denovo mutations. 
        - ancestral_mutations (dict): Dictionary to store ancestral mutations. 
        - tree: the Phylogenetic tree object.
        - ancestral_node_data (dict): Dictionary to store ancestral node data.

    Returns:
        - denovo_mutations (dict): Updated dictionary of denovo mutations.
        - ancestral_mutations (dict): Updated dictionary of ancestral mutations.
        - ancestral_node_data (dict): Updated dictionary of ancestral node data for each leaf node.
    """
    # Get the leaf node names from the alignment file, ignoring the outgroup sequence    
    leaf_nodes = [record.id for record in SeqIO.parse(input_msa_fasta_path, "fasta") if record.id != outgroup_id]

    #Compute the mutations for each leaf node
    for leaf_node in leaf_nodes:
        #find its nearest ancestral node data
        ancestral_node = _find_ancestral_node(tree, leaf_node)
        #store the ancestral node data for the leaf node
        ancestral_node_data[leaf_node] = ancestral_node
        #extract the sequence of this leaf node
        sequence = next((str(record.seq).upper() for record in SeqIO.parse(input_msa_fasta_path, "fasta") if record.id == leaf_node), None)
        #compute the mutations against the nearest ancestral sequence
        denovo_mutations[leaf_node] = find_mutations_wrt_subject_sequence(
            ancestral_sequences[str(ancestral_node).split("/")[0]], 
            sequence)
        #bootstrapping during tree construction can modify the names of the nodes in the tree by including the bootstrap results, hence we need to extract only the node name from the tree
        #and use it to save mutation data for that leafnode
        ancestral_mutations[leaf_node] = denovo_mutations[str(ancestral_node).split("/")[0]] + ancestral_mutations[str(ancestral_node).split("/")[0]]
    return denovo_mutations, ancestral_mutations, ancestral_node_data



def detect_mutations_phylogeny_based(input_msa_fasta_path, output_dir, tree_file_path, reconstructed_ancestral_state_file_path, outgroup_id, include_ancestral_node=False):
    """ 
    Detect mutations in a phylogeny-based manner using an input multiple sequence alignment (MSA), 
    a phylogenetic tree, and reconstructed ancestral state sequences.

    This function analyzes mutations for both internal (ancestral) and leaf nodes in a given tree 
    and outputs mutation data while considering an outgroup for rooting. If `include_ancestral_node` 
    is set to False, mutations from internal nodes are excluded.

    Args:
        input_msa_fasta_path (str or Path): Path to the input MSA FASTA file.
        output_dir (str or Path): Directory to store the output mutation files.
        tree_file_path (str or Path): Path to the phylogenetic tree file in Newick format.
        reconstructed_ancestral_state_file_path (str or Path): Path to the ancestral state reconstruction file.
        outgroup_id (str): Identifier for the outgroup sequence in the tree.
        include_ancestral_node (bool, optional): Whether to include mutation data for ancestral nodes. Defaults to False.

    Returns:
        tuple: 
            - denovo_mutations (dict): Dictionary mapping sequences to their de novo mutations.
            - ancestral_mutations (dict): Dictionary mapping sequences to their ancestral mutations.
    """
    
    # Read the tree file
    tree = Phylo.read(tree_file_path, "newick")

    # we now reconstruct the ancestral sequences of the tree from the ancestral state file
    print(f"Recontructing the ancestral sequences")
    ancestral_sequences=_get_reconstructed_ancestral_sequences(reconstructed_ancestral_state_file_path)

    # We now compute the mutations for the internal nodes of the tree
    print(f"Computing mutations for the ancestral nodes of the tree")
    denovo_mutations, ancestral_mutations,ancestral_node_data=_get_mutation_data_for_ancestral_nodes(ancestral_sequences, input_msa_fasta_path, tree, outgroup_id)

    #we now compute the mutations for the leaf nodes of the tree   
    print(f"Computing mutations for the leaf nodes of the tree")
    denovo_mutations, ancestral_mutations,ancestral_node_data=_get_mutation_data_for_leaf_nodes(input_msa_fasta_path, ancestral_sequences, denovo_mutations,
                                                                                                 ancestral_mutations, tree, ancestral_node_data, outgroup_id)

    # If include_ancestral_nodes flag is set to false, we remove the ancestral nodes mutation data
    if not include_ancestral_node:
        denovo_mutations = {key: value for key, value in denovo_mutations.items() if not key.startswith("Node")}
        ancestral_mutations = {key: value for key, value in denovo_mutations.items() if not key.startswith("Node")}


    print(f"Saving the mutation files to {output_dir}")
    # Save denovo mutations
    save_dict_to_pickle(denovo_mutations, f"{output_dir}/denovo_mutations.pkl")

    # Save ancestral mutations
    save_dict_to_pickle(ancestral_mutations, f"{output_dir}/ancestral_mutations.pkl")


    return denovo_mutations, ancestral_mutations