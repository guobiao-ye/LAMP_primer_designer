import subprocess
import pandas as pd
import ast

def create_blast_db(fasta_file, db_name):
    """
    Build a BLAST database using makeblastdb
    :param fasta_file: Path to the input FASTA file
    :param db_name: Name of the output database (without extension)
    """
    try:
        # Build BLAST database
        subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl', '-out', db_name], check=True)
        print(f"BLAST database {db_name} successfully created")
    except subprocess.CalledProcessError as e:
        print("*******")

def run_blast(db_name, output_file):
    """
    Run BLAST alignment using blastn
    :param query_sequence: Input query sequence (single)
    :param db_name: BLAST database name
    :param output_file: File for output results
    """
    try:
        # Run BLAST alignment
        subprocess.run(['blastn', '-query', './data/output/Intermediate_file/query.fasta', '-db', db_name, '-out', output_file, '-outfmt', '6', "-task", "blastn-short", "-word_size", "7", "-evalue", "1","-num_threads", "8"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST alignment failed: {e}")

def parse_blast_output(output_file):
    """
    Parse BLAST alignment results
    :param output_file: BLAST alignment result file
    :return: Parsed alignment results (list)
    """
    row_record=[]
    try:
        with open(output_file, "r") as f:
            for line in f:
                columns = line.split()
                row_now=int(columns[0].split(",")[0])
                primer_len=int(columns[0].split(",")[1])
                if row_now not in row_record and int(float(columns[2]))==100 and int(columns[3])==primer_len:
                    row_record.append(row_now)
        return row_record
    except Exception as e:
        print(f"Failed to parse: {e}")

def convert_to_dict(x):
    return ast.literal_eval(x)

def specificity_screening():
    fasta_file = "./data/resource/hg38.fa"  # Path to the input FASTA file
    db_name = "./data/resource/blast_db"  # BLAST database name
    output_file = "./data/output/Intermediate_file/blast_output.txt"  # Path to the BLAST output file

    df=pd.read_csv("./data/output/Intermediate_file/candidate_result.csv")
    columns = df.columns
    info_columns = [col for col in columns if 'info' in col]
    row,col=df.shape
    df = pd.read_csv("./data/output/Intermediate_file/candidate_result.csv", converters={col: convert_to_dict for col in info_columns})

    # 1. Build BLAST database
    create_blast_db(fasta_file, db_name)

    # 2. Execute BLAST alignment
    with open("./data/output/Intermediate_file/query.fasta", "w") as query_file:
        for i in range(row):
            for j in range(1,9):
                dict=df.iat[i,j]
                seq=dict["sequence"]
                length=dict["length"]
                query_file.write(f">{i},{length},{j}\n{seq}\n")
    run_blast(db_name,output_file)
    row_sub = parse_blast_output(output_file)
    result=df.drop(row_sub) 
    result.to_csv("./data/output/Intermediate_file/specific_primer.csv",index=False)