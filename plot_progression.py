import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file_path):
    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=['Nb kmers seen', 'Nb unique k-mers', 'growth', 'acceleration'])

    # Read the file and parse the lines
    # one line every 100 lines
    nb_lines = 0
    with open(input_file_path, 'r') as file:
        for line in file:
            nb_lines += 1
            if line.startswith('Nb kmers seen:') and nb_lines % 100 == 0:
                parts = line.split()
                Nb_kmers_seen = int(parts[3])
                nb_unique_kmers = int(parts[7].rstrip(','))
                growth = float(parts[9].rstrip(','))
                acceleration = float(parts[11])
                df.loc[len(df)] = [Nb_kmers_seen, nb_unique_kmers, growth, acceleration]

    # Get the file name without extension for the plot title
    file_name = input_file_path.split('/')[-1].split('.')[0]
    # Plotting to f'{file_name}_kmers_plot.png'
    plt.figure(figsize=(10, 6))
    plt.plot(df['Nb kmers seen'], df['Nb unique k-mers'], marker='o')
    plt.title(f'{file_name} Nb distinct solid canonical k-mers vs Nb total of kmer')
    plt.xlabel('Number of kmers seen')
    plt.ylabel('Number of distinct solid canonical k-mers')
    plt.grid(True)
    plt.savefig(f'{file_name}_kmers_plot.png')
    print(f"Plot saved as {file_name}_kmers_plot.png")
    
    

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Plot unique k-mers vs reads from a pukc log file.')
    parser.add_argument('input_file', type=str, help='Path to the input log file')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the main function with the input file path
    main(args.input_file)
