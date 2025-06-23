import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file_path):
    # Initialize an empty DataFrame
    df = pd.DataFrame(columns=['Nb kmers seen', 'Nb unique k-mers'])
    # A  line is expected to be like:
    # Processed 190000 reads, 13994639 kmers seen, 3008814 unique solid kmers.
    # Read the file and parse the lines
    # one line every 10 lines
    nb_lines = 0
    last_nb_kmers_seen = 0
    last_nb_unique_kmers = 0
    total_estimated_kmers = 0
    growth = 0.0
    with open(input_file_path, 'r') as file:
        for line in file:
            nb_lines += 1
            if line.startswith('Processed'):
                nb_lines += 1
                if nb_lines % 10 == 0:
                    parts = line.split()
                    Nb_kmers_seen = int(parts[3])
                    nb_unique_kmers = int(parts[6].rstrip(','))
                    df.loc[len(df)] = [Nb_kmers_seen, nb_unique_kmers]
                    last_nb_kmers_seen = Nb_kmers_seen
                    last_nb_unique_kmers = nb_unique_kmers
            if line.startswith('Estimated total kmers'):
                parts = line.strip().split()
                if parts[-1].isdigit():
                    # line: Estimated total kmers: 13994639
                    # get the total estimated kmers
                    total_estimated_kmers = int(parts[-1])
            # line: Stable curve found 10/10: f(kmers) = ? + 0.10786222*kmers + -0.000055000186*kmers^2
            # get the growth from the "*kmer" factor, here 0.10786222
            if line.startswith('Stable curve found'):
                part = line.strip().split()[-3]
                growth = float(part.split('*')[0])
             

    # Get the file name without extension for the plot title
    file_name = input_file_path.split('/')[-1].split('.')[0]
    # Plotting to f'{file_name}_kmers_plot.png'
    plt.figure(figsize=(10, 6))
    # plot log x values: 
    # plt.xscale('log')
    # plt.yscale('log')
    plt.plot(df['Nb kmers seen'], df['Nb unique k-mers'], marker='o')
    plt.title(f'{file_name} Nb distinct solid canonical k-mers vs Nb total of kmer')
    plt.xlabel('Number of kmers seen')
    plt.ylabel('Number of distinct solid canonical k-mers')
    plt.grid(True)
    
    # Plotting on the same graph the tendency line: with x from last_nb_kmers_seen to total_estimated_kmers
    # with the plot: f(x)= last_nb_kmers_seen + growth * (x - last_nb_kmers_seen)
    range_estimation = last_nb_kmers_seen // 3 # 33% of the last nb kmers seen
    total_estimated_kmers = last_nb_kmers_seen + range_estimation
    x_values = range(last_nb_kmers_seen, total_estimated_kmers + 1, 100000)
    y_values = [last_nb_unique_kmers + growth * (x - last_nb_kmers_seen) for x in x_values]
    plt.plot(x_values, y_values, color='red', linestyle='--', label='Estimated growth trend')
    
    # also plot the estimated f(x) from last_nb_kmers_seen - range_estimation
    x_values_estimated = range(last_nb_kmers_seen - range_estimation, last_nb_kmers_seen + 1, 100000)
    y_values_estimated = [last_nb_unique_kmers + growth * (x - last_nb_kmers_seen) for x in x_values_estimated]
    plt.plot(x_values_estimated, y_values_estimated, color='green', linestyle='--', label='Estimated growth trend (before)')
    
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
