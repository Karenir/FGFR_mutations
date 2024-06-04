import pandas as pd

file_path = 'FGFR3_mutations_cosmic.csv' #You might have to change the file name to the path where you store FGFR3_mutations_cosmic.csv

# Read the CSV file into a DataFrame 
df = pd.read_csv(file_path)

# Filter out silent mutations] and unknown mutations
non_silent_mutations = df[(df['Type'] != 'Substitution - coding silent') & (df['Type'] != 'Unknown') & (df['Type'] != 'Nonstop extension')]

# Count of all mutations (excluding silent mutations and unknown)
print("all mutations that are not silent", non_silent_mutations ['Count'].sum()) #871

#Mutations influencing S776
stop_codon_before_776 = df[(df['Position'] <= 776) & (df ['Type'] == 'Substitution - Nonsense')]
frameshift_before_776 = df[(df ['Position'] <= 776) & ((df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift'))]
mis_mut_S765 = df[(df ['Position'] == 765) & (df['Type'] == 'Substitution - Missense')]
mis_mut_S776 = df[(df ['Position'] == 776) & (df['Type'] == 'Substitution - Missense')]

#Add count.sum because there could be more counts per mutation
print("Stop codons before 776: ", stop_codon_before_776 ['Count'].sum()) #40
print("Frameshift before 776: ", frameshift_before_776 ['Count'].sum()) #41
print("Pointmutations on S765:", mis_mut_S765 ['Count'].sum()) #4
print("Pointmutations on S776:", mis_mut_S776 ['Count'].sum()) #1
print("Insertions/Deletions in frame influencing S765 and/or S776: 0") #manual count

#All mutations
all_stop_codons = df [df ['Type'] == 'Substitution - Nonsense']
all_frameshifts = df[(df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Complex - frameshift')]
all_deletionsInsertion_inframe = df[(df['Type'] == 'Insertion - In frame') | (df['Type'] == 'Deletion - In frame') | (df['Type'] == 'Complex - deletion inframe')]
all_missense = df [df ['Type'] == 'Substitution - Missense']
all_nonstop_extension =  df [df ['Type'] == 'Nonstop extension']

print("total stop codons: ", all_stop_codons ['Count'].sum()) #40
print("total frameshifts: ", all_frameshifts ['Count'].sum()) #44
print("total deletions/insertion in frame: ", all_deletionsInsertion_inframe ['Count'].sum()) #0
print("total missense: ", all_missense ['Count'].sum()) #787
print("all_nonstop_extension: ", all_nonstop_extension  ['Count'].sum()) #1

#All mutations with catalytic site, defined to 755 aa
all_stop_codons_cat = df[(df['Position'] >= 755) & (df['Type'] == 'Substitution - Nonsense')] 
all_frameshifts_cat = df[(df['Position'] >= 755) & ((df ['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site: ", all_stop_codons_cat ['Count'].sum()) #0
print("total frameshifts w/cat site: ", all_frameshifts_cat ['Count'].sum()) #0
print("total deletions/insertion in frame w/cat site: ", all_deletionsInsertion_inframe ['Count'].sum()) #0

#Mutations influencing S776 with cat site:
all_stop_codons_cat_tail = df[(df ['Position'] >= 755) & (df['Position'] <= 776) & (df ['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat_tail = df[(df ['Position'] >= 755) & (df['Position'] <= 776) & ((df ['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site influencing S776: ", all_stop_codons_cat_tail['Count'].sum()) #0
print("total frameshifts w/cat site influencing S776: ", all_frameshifts_cat_tail['Count'].sum()) #0
print("Insertions/Deletions in frame with catalytic site influencing S765 and/or S776: 0") #manual count
