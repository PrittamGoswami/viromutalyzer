from viromutalyzer import detect_mutations_reference_based, detect_mutations_phylogeny_based
import typer
from pathlib import Path

app = typer.Typer()

@app.command()
def analyze(
    input_msa_file_path: Path = typer.Option(..., "-i", "--input", help="Input Multiple Sequence ALigned FASTA file"),
    output_dir: Path = typer.Option(..., "-o", "--output", help="Output directory"),
    ref_based_mutation_detection: bool = typer.Option(
        True, 
        "--ref-based-mutation-detection", 
        help="Enable reference-based mutation detection (default: True)"
    ),
    phylogeny_based_mutation_detection: bool = typer.Option(
        False, 
        "--phylogeny-based-mutation-detection", 
        help="Enable phylogeny-based mutation detection (default: False)"
    ),
    outgroup_id: Path = typer.Option(None, "-og", "--outgroup_id", help="Outgroup Sequence ID"),
    tree_file_path: Path = typer.Option(None, "-t", "--tree_file", help="Path to the tree file"),
    ancestral_state_file_path: Path = typer.Option(None, "-as", "--ancestral_state", help="Path to the Ancestral state file"),
    include_ancestral_nodes: bool = typer.Option(
        False, 
        "--include_ancestral_nodes", 
        help="Include the ancestral nodes in the results (default: False)"
    ),
):
    """
    Analyze mutations in the given FASTA file and save results to the output directory.
    """
    # Check if the input MSA fasta file exits otherwise raise an error
    if not input_msa_file_path.exists():
        typer.echo("Error: Input file does not exist!", err=True)
        raise typer.Exit(1)
    # If the output directory doesnt exist create it
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)


    # If the user does not want to use the phylogeny based approach then use the reference-based approoach for mutation detection
    if (not phylogeny_based_mutation_detection) and (ref_based_mutation_detection): 

        typer.echo(f"Analyzing mutations in {input_msa_file_path} and saving results to {output_dir}...")
        mutation_dict = detect_mutations_reference_based(input_msa_file_path, output_dir)
    
    else:
        if not tree_file_path.exists():
            typer.echo("Error: Input tree file does not exist!", err=True)
            raise typer.Exit(1)
        if not ancestral_state_file_path.exists():
            typer.echo("Error: Input Ancestral state file does not exist!", err=True)
            raise typer.Exit(1)
         
        #Phylogeny based mutation detection
        typer.echo(f"Detecting mutations in {input_msa_file_path} and saving results to {output_dir}...")
        denovo_mutations, ancestral_mutations = detect_mutations_phylogeny_based(input_msa_file_path, output_dir, tree_file_path, ancestral_state_file_path, outgroup_id, include_ancestral_nodes)


def main():
    app()

if __name__ == "__main__":
    main()
