af_pdb_path = "./af.pdb"
rmsf_xvg_path = "./run2/rmsf.xvg"
res_limit = 161

# Read in .pdb file and parse secondary structure elements
helix_lines = []
sheet_lines = []
ca_lines = []

# Open the file and read line by line
with open(af_pdb_path, 'r') as file:
    for line in file:
        # Check if the line starts with "HELIX"
        if line.startswith("HELIX"):
            helix_lines.append(line.strip())
        elif line.startswith("SHEET"):
            sheet_lines.append(line.strip())
        elif line.startswith("ATOM"):
            if "CA" in line:
                if len(ca_lines) < res_limit:
                    ca_lines.append(line.strip())

helix_residues = []
sheet_residues = []

for line in helix_lines:
    helix_residues.append([int(line[21:25].strip()), int(line[33:37].strip())])

for line in sheet_lines:
    sheet_residues.append([int(line[22:26].strip()), int(line[34:38].strip())])

pdb_cas = []

for line in ca_lines:
    pdb_cas.append([line[7:12].strip(), line[17:20].strip()])

# Read in rmsf.xvg file and parse rmsf for c-alphas
rmsf_values = []
rmsf_values_2 = []

with open(rmsf_xvg_path, 'r') as file:
    atoms_now = 0
    for line in file:
        # Check when the atoms start
        if atoms_now == 1:
            rmsf_values.append(float(line.split()[1]))
        elif line.startswith("@TYPE xy"):
            atoms_now = 1

# Check that there are the same number of c-alphas as residues

print(len(rmsf_values), len(pdb_cas))

# Prepare final data for plotting 

residue_ss = [0 for i in range(len(pdb_cas))]
# 0: no ss, 1: alpha helix, 2: beta strand 
for element in helix_residues:
    length = element[1] - element[0] + 1
    for i in range(length):
        res = i + element[0] - 1
        if res < res_limit:
            residue_ss[res] = 1

for element in sheet_residues:
    length = element[1] - element[0] + 1
    for i in range(length):
        res = i + element[0] - 1
        if res < res_limit:
            residue_ss[res] = 2

final_data = [[], [], [], []] # [ res_num, res_name, ss_element, rmsf ]

for i in range(len(pdb_cas)):

   # final_data.append([i+1, pdb_cas[i][1], residue_ss[i], rmsf_values[i]])
    final_data[0].append(i+1)
    final_data[1].append(pdb_cas[i][1])
    final_data[2].append(residue_ss[i])
    final_data[3].append(rmsf_values[i])

# Plot rmsf.xvg for c-alphas (as residue numbers)


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


# Define the window size for a centered rolling average of +-5 residues
half_window = 5
window_size = 2 * half_window + 1  # This results in a window of size 11

# Calculate the centered rolling average
rolling_avg = np.convolve(final_data[3], np.ones(window_size) / window_size, mode='same')

rolling_avg[:window_size//2] = np.nan  # Set first 5 (or half the window size) values to NaN
rolling_avg[-window_size//2:] = np.nan # Set last 5 (or half the window size) values to NaN


# Note that using `mode='same'` keeps the output length the same as the original data,
# so rolling_avg_x will be the same as final_data[0].

fig, ax = plt.subplots(figsize=(10, 6))

# Plot RMSF values
ax.plot(final_data[0], rmsf_values, label="RMSF for C-alpha", color="black", linewidth=1)


# Plot the centered rolling average line
#ax.plot(final_data[0], rolling_avg, label="Rolling Average (+-5)", color="orange", linewidth=1.5)

# Define the y position for secondary structure visualization
ss_y_position = -0.05  # Position below RMSF plot
ss_height = 0.05       # Height of the rectangles/arrows for visualization

ax.plot(final_data[0], [ss_y_position + (ss_height * 0.5)] * len(final_data[0]), color="black", linestyle="-", linewidth=1, zorder=-1)

# Add rectangles for alpha helices
for element in helix_residues:
    start = element[0]  # Starting residue number
    end = element[1]    # Ending residue number
    if start > res_limit:
        break
    elif end > res_limit:
        end = res_limit
    width = end - start + 1
    rect = patches.Rectangle(
        (start, ss_y_position), width, ss_height,
        edgecolor="black", facecolor='lightblue', label="Alpha Helix (AlphaFold 3)" if 'Alpha Helix (AlphaFold 3)' not in ax.get_legend_handles_labels()[1] else ""
    )
    ax.add_patch(rect)

# Add arrows for beta sheets
for element in sheet_residues:
    start = element[0]
    end = element[1]
    if start > res_limit:
        break
    elif end > res_limit:
        end = res_limit
    width = end - start + 1
    # Create an arrow for the beta sheet
    arrow = patches.FancyArrow(
        start, ss_y_position + ss_height / 2, width - 0.5, 0,  # Arrow length and position
        edgecolor='black' ,width=ss_height / 2, facecolor="red", head_width=ss_height, head_length=2.0, length_includes_head=True, label="Beta Sheet (AlphaFold 3)" if 'Beta Sheet (AlphaFold 3)' not in ax.get_legend_handles_labels()[1] else ""
    )
    ax.add_patch(arrow)

# Add labels and legend
ax.set_xlabel("Cα Residue Number")
ax.set_ylabel("RMSF (nm)")
ax.legend()
plt.title("TrmD Monomer Cα RMSF (Run 2)")

plt.savefig("rmsf_plot.png", format="png", dpi=300)

plt.show()
