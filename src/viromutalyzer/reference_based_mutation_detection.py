from Bio import SeqIO
from itertools import product
import pickle
import multiprocessing

__all__ = ["mutation_detection_wrt_reference", "save_dict_to_pickle", "get_kmer_frequency", "find_mutations_wrt_subject_sequence"]


#Function to compute the trinucleotide frequency of a sequence to be used in normalization step while computing the trinucleotide mutation rate
def get_kmer_frequency(sequence, k = 3):
    """
    Calculate the frequency of k-mers in a given sequence.

    Parameters:
        - sequence (str): The input RNA sequence.
        - k (int): The length of the k-mer. Default is 3 

    Returns:
        - kmers_frequency_dict: A dictionary containing all possible k-mers as keys and their frequencies in the input sequence as values.
    """
    # Remove gaps from the sequence if any
    sequence = sequence.replace("-", "")
    #Total number of kmers in the sequence
    total_kmers = len(sequence) - k + 1

    # Declare and Initialize the dictionary to store the kmer frequency by assigning all possible k-mers with frequency 0
    nucleotides = ["C", "T", "G", "A"]
    kmers_frequency_dict = {''.join(combo): 0 for combo in product(nucleotides, repeat=k)}

    # Iterate through the sequence to find all k-mers
    for i in range(total_kmers):
        kmer = sequence[i:i + k]  # Extract the k-mer
        # Check if the k-mer consists only of A, T, G, or C and not any other characters
        if set(kmer).issubset(set(nucleotides)):
            # Increment the frequency of this k-mer in the kmer_frequency dictionary
            kmers_frequency_dict[kmer] += 1  
    # Return the kmer dictionary
    return kmers_frequency_dict


def find_mutations_wrt_subject_sequence(subject_seq, sample_seq):
    """
    To compute the mutations by comparing each sequence against its nearest reconstructed ancestral sequence

    Parameters:
        - subject_seq (str): The ancestor sequence.
        - sample_seq (str): The sample sequence.

    Returns:
        - mutations (list): A list of mutations, where each mutation is represented as a list containing the mutation, the trinucleotide context, and the trinucleotide frequency
    """
    # Initialize the list to store the mutations
    mutations = []
    # Get the trinucleotide mutation frequencies in the ancestor sequence
    trinucleotide_mutation_frequency = get_kmer_frequency(subject_seq, 3)
    # Iterate through each position and bases in the two sequences
    for i, (a, b) in enumerate(zip(subject_seq, sample_seq)):
        # If the base is different between the ancestor and sample sequences 
        if a != b:
            #position of the mutation in the alignment is index + 1
            position = i + 1  
            # Extract the trinucleotide context of the mutation
            trinucleotide = str(subject_seq[position-2:position+1]).upper()  
            # Base in the ancestor sequence
            ref_base = a.upper()
            # Base in the sample sequence
            sample_base = b.upper()
            # We ignore those mutations at the starting or end position of the alignment as their trinucleotide context cannot be captured
            # We also ignore those mutations in which the surrounding trinucleotide contains abnormal basses or the mutation itself contains abnormal bases
            if (
                any(base not in ["C", "T", "G", "A"] for base in trinucleotide) or
                any(base not in ["C", "T", "G", "A"] for base in [ref_base, sample_base]) or
                len(trinucleotide) != 3
            ):                
                continue
            else:
                # Record the mutation data which includes the mutation, its trinucleotide context and trinuleotide frequency 
                mutation = f"{ref_base}{position}{sample_base}".upper()
                mutations.append([mutation, trinucleotide, trinucleotide_mutation_frequency[trinucleotide]])
               
    return mutations





def save_dict_to_pickle(dictionary, out_dir):
    """
        Saves a dictionary to a file using pickle.
    
    Parameters:
        - dictionary(dictionary object): The dictionary to save
        - file_path(str): The filename to save the dictionary to
    """
    with open(f"{out_dir}/mutation.pkl", 'wb') as file:
        pickle.dump(dictionary, file)



def _detect_mutations_batchwise(batch, reference_sequence):
    """
    Process a batch of sequences to detect mutations in parallel.
    
    Parameters:
        - batch (list): A list of SeqRecord objects to be processed.
        - reference_sequence (str): The reference genome sequence.
    
    Returns:
        - batch_mutations_dict (dict): A dictionary mapping sequence IDs to their mutation data.
    """
    # Dictionary to record mutation data for this batch
    batch_mutations_dict = {}
    # Iterate through each record in this batch
    for record in batch:
        # Extract the sample sequence and ID
        sample_sequence = str(record.seq)
        sample_sequence_id = record.id
        # Compute the mutation
        batch_mutations_dict[sample_sequence_id] = find_mutations_wrt_subject_sequence(reference_sequence, sample_sequence)
    return batch_mutations_dict


def detect_mutations_reference_based(input_msa_fasta_path,out_dir, ref_seq_id=None, num_threads=1):
    """
    Detect mutations by comparing sequences against the reference genome using multiprocessing.

    Parameters:
        - aligned_file_path (str): The path to the MSA input file.
        - ref_seq_id (str): The ID of the reference genome.
        - num_threads (int): The number of processes to use.

    Returns:
        - mutations_dict (dict): A dictionary containing mutation data for each genome.
    """
    # Load sequences
    records = list(SeqIO.parse(input_msa_fasta_path, "fasta"))

    # Determine the reference sequence
    if ref_seq_id:
        reference_record = next((record for record in records if record.id == ref_seq_id), None)
        if not reference_record:
            raise ValueError(f"Reference sequence with ID '{ref_seq_id}' not found in the input FASTA file.")
        reference_sequence = str(reference_record.seq)
    else:
        reference_record = records[0]
        reference_sequence = str(reference_record.seq)

    # Split records into batches for multiprocessing
    batch_size = max(1, len(records) // num_threads)  # Avoid zero batch size
    batches = [records[i:i + batch_size] for i in range(0, len(records), batch_size)]

    # Use multiprocessing to process batches
    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.starmap(_detect_mutations_batchwise, [(batch, reference_sequence) for batch in batches])

    # Merge results into a final dictionary
    mutations_dict = {}
    for batch_result in results:
        mutations_dict.update(batch_result)

    save_dict_to_pickle(mutations_dict, out_dir)

    # Return the final mutation dictionary
    return mutations_dict

