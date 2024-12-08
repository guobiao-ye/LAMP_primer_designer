import subprocess
import pandas as pd
import ast
from Bio import SeqIO

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

def parse_blast_output_withinput(output_file,input_record):
    """
    Parse BLAST alignment results with input sequence information
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
                    if not(columns[1]==input_record[0] and columns[8]>=input_record[1] and columns[9]<=input_record[2]):
                        row_record.append(row_now)
        return row_record
    except Exception as e:
        print(f"Failed to parse: {e}")

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

def align_input_seq(input_seq,db_name,input_align_output):
    input_record=[]
    with open("./data/output/Intermediate_file/input_seq.fasta","w") as f:
        f.write(f">input\n{input_seq}")
    try:
        subprocess.run(['blastn', '-query', './data/output/Intermediate_file/input_seq.fasta', '-db', db_name, '-out', input_align_output, '-outfmt', '6',"-num_threads", "8"], check=True)
        #subprocess.run(['blastn', '-query', 'query.fasta', '-db', db_name, '-out', output_file, '-outfmt', '6',"-task", "blastn-short", "-word_size", "7", "-evalue", "1","-num_threads", "8"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST alignment failed: {e}")
    try:
        with open(input_align_output, "r") as f:
            for line in f:
                columns = line.split()
                if int(float(columns[2]))==100 and int(columns[3])==len(input_seq):
                    input_record.append(columns[1])
                    input_record.append(columns[8])
                    input_record.append(columns[9])
        return input_record
    except Exception as e:
        print(f"Failed to parse: {e}")

def specificity_screening(input_file = './data/input/example.fasta', fasta_file = "./data/resource/hg38.fa"):
    db_name = "./data/resource/blast_db"  # BLAST database name
    output_file = "./data/output/Intermediate_file/blast_output.txt"  # Path to the BLAST output file
    input_align_output="./data/output/Intermediate_file/input_align_output.txt"

    input_seq = None
    quasispecies_sequences = []

    # Open the FASTA file and parse the records
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if input_seq is None:
                # Set the first sequence as the main sequence
                input_seq = str(record.seq)
            else:
                # Add subsequent sequences to the quasispecies list
                quasispecies_sequences.append(str(record.seq))

    df=pd.read_csv("./data/output/Intermediate_file/candidate_result.csv")
    columns = df.columns
    info_columns = [col for col in columns if 'info' in col]
    row,col=df.shape
    df = pd.read_csv("./data/output/Intermediate_file/candidate_result.csv", converters={col: convert_to_dict for col in info_columns})

    # 1. Build BLAST database
    create_blast_db(fasta_file, db_name)

    input_record=align_input_seq(input_seq,db_name,input_align_output)

    # 2. Execute BLAST alignment
    with open("./data/output/Intermediate_file/query.fasta", "w") as query_file:
        for i in range(row):
            for j in range(1,9):
                dict=df.iat[i,j]
                seq=dict["sequence"]
                length=dict["length"]
                query_file.write(f">{i},{length},{j}\n{seq}\n")
    run_blast(db_name,output_file)
    if len(input_record)>0:
        row_sub = parse_blast_output_withinput(output_file,input_record)
    else:
        row_sub = parse_blast_output(output_file)
    result=df.drop(row_sub) 
    result.to_csv("./data/output/Intermediate_file/specific_primer.csv",index=False)