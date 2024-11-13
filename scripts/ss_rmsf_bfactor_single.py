xtal_pdb_path = "./1p9p.pdb"
af_pdb_path = "./af.pdb"
rmsf_xvg_path = "./run1/rmsf.xvg"
residue_limit = 161

helix_lines = []
sheet_lines = []
ca_lines = []
rmsf_values = []

# Populate ss
with open(af_pdb_path, 'r') as file:
    prev_residue_num = 0
    for line in file:
        # Check if the line starts with "HELIX"
        if line.startswith("HELIX"):
            helix_lines.append(line.strip())
        elif line.startswith("SHEET"):
            sheet_lines.append(line.strip())
        elif line.startswith("ATOM"):
            if "CA" in line:
                if int(line[22:26].strip()) == (prev_residue_num + 1):
                    if len(ca_lines) < residue_limit:
                        ca_lines.append(line.strip())
                else:
                    if len(ca_lines) < residue_limit:
                        ca_lines.append("missing")

helix_residues = []
sheet_residues = []

for line in helix_lines:
    helix_residues.append([int(line[21:25].strip()), int(line[33:37].strip())])

for line in sheet_lines:
    sheet_residues.append([int(line[22:26].strip()), int(line[34:38].strip())])

pdb_cas = []

for line in ca_lines:
    if line == "missing":
        pdb_cas.append(0)
    else:
        pdb_cas.append([line[7:12].strip(), line[17:20].strip()])