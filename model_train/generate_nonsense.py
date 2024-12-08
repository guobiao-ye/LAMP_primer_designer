import random
import pandas as pd

# Function to generate a random DNA sequence of a given length
def generate_random_sequence(length):
    return ''.join(random.choices('ACTG', k=length))

# Function to generate random positions that correspond to the sequence length
def generate_random_position(sequence_length):
    start = random.randint(1, 100000)  # random starting position (1 to 100000)
    end = start + sequence_length - 1
    return f'{start}-{end}'

# Generate the list of rows
def generate_nonsense_primers(num_primers):
    data = []
    for i in range(1, num_primers + 1):
        base_lamp_id = f'NS{i:05d}'  # Base Nonsense primer ID (NS00001, NS00002, ...)
        genbank = base_lamp_id  # Genbank ID is the same as lamp_id
        
        # Iterate over both F3-F1 and B1-B3
        for name in ['F3', 'F2', 'F1', 'B1', 'B2', 'B3']:
            lamp_id = f'{base_lamp_id}'  # E.g., NS00001_F1, NS00001_B2, etc.
            sequence_length = random.randint(18, 25)  # Random length between 18 and 25
            sequence = generate_random_sequence(sequence_length)
            position = generate_random_position(sequence_length)
            
            row = [i, genbank, name, sequence, position, lamp_id]  # Adding id as the first column
            data.append(row)

    # Create a DataFrame and save to CSV
    df = pd.DataFrame(data, columns=['id', 'genbank', 'name', 'sequence', 'position', 'lamp_id'])
    return df

# Generate 10 base primers, each having 6 variants (F1-F3 and B1-B3) and save to a CSV file
nonsense_primers_df = generate_nonsense_primers(415)
nonsense_primers_df.to_csv('data/nonsense_primers_10w.csv', index=False)

print(nonsense_primers_df)
