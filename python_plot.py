import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file_path):
    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=['Nb reads', 'Nb unique k-mers', 'growth', 'acceleration'])

    # Read the file and parse the lines
    # one line every 100 lines
    nb_lines = 0
    with open(input_file_path, 'r') as file:
        for line in file:
            nb_lines += 1
            if line.startswith('Nb reads:') and nb_lines % 100 == 0:
                parts = line.split()
                nb_reads = int(parts[2])
                nb_unique_kmers = int(parts[6].rstrip(','))
                growth = float(parts[8].rstrip(','))
                acceleration = float(parts[10])
                df.loc[len(df)] = [nb_reads, nb_unique_kmers, growth, acceleration]

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(df['Nb reads'], df['Nb unique k-mers'], marker='o')
    plt.title('Number of Unique k-mers vs Number of Reads')
    plt.xlabel('Number of Reads')
    plt.ylabel('Number of Unique k-mers')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Plot unique k-mers vs reads from a pukc log file.')
    parser.add_argument('input_file', type=str, help='Path to the input log file')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the main function with the input file path
    main(args.input_file)
