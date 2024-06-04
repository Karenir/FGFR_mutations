import pandas as pd

file_path = 'FGFR2_mutations_cosmic.csv' #You might have to change the file name to the path where you store FGFR2_mutations_cosmic.csv

# Read the CSV file into a DataFrame 
df = pd.read_csv(file_path)

# Filter out silent mutations] and unknown mutations
non_silent_mutations = df[(df['Type'] != 'Substitution - coding silent') & (df['Type'] != 'Unknown')]

# Count of all mutations (excluding silent mutations and unknown)
print("all mutations that are not silent", non_silent_mutations ['Count'].sum()) #1583

#Mutations influencing S791
stop_codon_before_791 = df[(df['Position'] <= 791) & (df ['Type'] == 'Substitution - Nonsense')]
frameshift_before_791 = df[(df ['Position'] <= 791) & ((df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift'))]
mis_mut_S780 = df[(df ['Position'] == 780) & (df['Type'] == 'Substitution - Missense')]
mis_mut_S791 = df[(df ['Position'] == 791) & (df['Type'] == 'Substitution - Missense')]

#Add count.sum because there could be more counts per mutation
print("Stop codons before 791: ", stop_codon_before_791 ['Count'].sum()) #42
print("Frameshift before 791: ", frameshift_before_791 ['Count'].sum()) #23
print("Point mutations on S780:", mis_mut_S780 ['Count'].sum()) #1
print("Point mutations on S791:", mis_mut_S791 ['Count'].sum()) #4
print("Insertions/Deletions in frame influencing S780 and/or S791: 0") #manual count

#All mutations
all_stop_codons = df [df ['Type'] == 'Substitution - Nonsense']
all_frameshifts = df[(df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift')]
all_deletionsInsertion_inframe = df[(df['Type'] == 'Insertion - In frame') | (df['Type'] == 'Deletion - In frame') | (df['Type'] == 'Complex - deletion inframe')]
all_missense = df [df ['Type'] == 'Substitution - Missense']

print("total stop codons: ", all_stop_codons ['Count'].sum()) #44
print("total frameshifts: ", all_frameshifts ['Count'].sum()) #25
print("total deletions/insertion in frame: ", all_deletionsInsertion_inframe ['Count'].sum()) #16
print("total missense: ", all_missense ['Count'].sum()) #1498

#All mutations with catalytic site, defined to 770 aa
all_stop_codons_cat = df[(df['Position'] >= 770) & (df['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat = df[(df['Position'] >= 770) & ((df ['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site: ", all_stop_codons_cat ['Count'].sum()) #4
print("total frameshifts w/cat site: ", all_frameshifts_cat ['Count'].sum()) #3
print("total deletions/insertion in frame w/cat site: ", all_deletionsInsertion_inframe ['Count'].sum()) #all 16 (checked manually)

#Mutations influencing S791 with cat site:
all_stop_codons_cat_tail = df[(df ['Position'] >= 770) & (df['Position'] <= 791) & (df ['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat_tail = df[(df ['Position'] >= 770) & (df['Position'] <= 791) & ((df ['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site influencing S791: ", all_stop_codons_cat_tail['Count'].sum()) #2
print("total frameshifts w/cat site influencing S791: ", all_frameshifts_cat_tail['Count'].sum()) #1
print("Insertions/Deletions in frame with catalytic site influencing S780 and/or S791: 0") #manual count
