with open('.jobs_fullyMerge_v2.txt', 'r') as f:
    lines = f.readlines()

unique_lines = set()
duplicates = []
for line in lines:
    if line.strip() not in unique_lines:
        unique_lines.add(line.strip())
    else:
        duplicates.append(line.strip())

if len(duplicates) == 0:
    print("No duplicate lines found.")
else:
    print("Duplicate lines found:")
    for line in duplicates:
        print(line)
