import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def extract_string(input_string, start_char, stop_char):
    start_index = input_string.find(start_char)
    stop_index = input_string.find(stop_char, start_index + 1)
    if start_index != -1 and stop_index != -1:
        return input_string[start_index:stop_index]
    return None

def select_last_energy(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('!'):
                number = extract_string(line, '-', ' Ry')
                if number is not None:
                    return float(number)
    return None

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 script.py path_to_files plot_saved_path [xlabel]")
        sys.exit(1)
    
    path_to_files = sys.argv[1]
    plot_saved_path = sys.argv[2]
    xlabel = sys.argv[3] if len(sys.argv) > 3 else "wf cutoff [Ry]"
    
    files = sorted(f for f in os.listdir(path_to_files) if f.endswith(".out"))
    
    grids = []
    energies = []
    
    for file_name in files:
        try:
            ecut_value = int(file_name.split('_')[1].split('.')[0])
            energy = select_last_energy(os.path.join(path_to_files, file_name))
            if energy is not None:
                grids.append(ecut_value)
                energies.append(energy)
        except ValueError:
            continue
    
    if not grids:
        print("No valid data found.")
        sys.exit(1)
    
    grids, energies = zip(*sorted(zip(grids, energies)))
    grids = np.array(grids)
    energies = np.array(energies)
    
    plt.figure(figsize=(16, 9))
    plt.plot(grids, energies, marker='o')
    ax = plt.gca()
    ax.yaxis.offsetText.set_fontsize(24)
    plt.ylabel('Energies [Ry]', fontsize=28)
    plt.xlabel(xlabel, fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid()
    plt.tight_layout()

    plt.savefig(plot_saved_path)
    plt.show()
    

