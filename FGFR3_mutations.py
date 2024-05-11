import pandas as pd

file_path = '/Volumes/KarenDisk/V24/Mutationsfigures/Gene_mutationsWed Mar 20 14_02_07 2024.csv'

# Read the CSV file into a DataFrame 
df = pd.read_csv(file_path)

# Filter out silent mutations] and unknown mutations
non_silent_mutations = df[(df['Type'] != 'Substitution - coding silent') & (df['Type'] != 'Unknown') & (df['Type'] != 'Nonstop extension')]

# Count of all mutations (excluding silent mutations and unknown)
print("all mutations that are not silent", non_silent_mutations ['Count'].sum()) #5241

#Mutations influencing S791
stop_codon_before_782 = df[(df['Position'] <= 782) & (df ['Type'] == 'Substitution - Nonsense')]
frameshift_before_782 = df[(df ['Position'] <= 782) & ((df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift'))]
mis_mut_S771 = df[(df ['Position'] == 771) & (df['Type'] == 'Substitution - Missense')]
mis_mut_S782 = df[(df ['Position'] == 782) & (df['Type'] == 'Substitution - Missense')]

#Add count.sum because there could be more counts per mutation
print("Stop codons before 791: ", stop_codon_before_782 ['Count'].sum()) #10
print("Frameshift before 791: ", frameshift_before_782 ['Count'].sum()) #59
print("Pointmutations on S771:", mis_mut_S771 ['Count'].sum()) #0
print("Pointmutations on S782:", mis_mut_S782 ['Count'].sum()) #0
print("Insertions/Deletions in frame influencing S771 and/or S782: 0") #manual count

#All mutations
all_stop_codons = df [df ['Type'] == 'Substitution - Nonsense']
all_frameshifts = df[(df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Complex - frameshift')]
all_deletionsInsertion_inframe = df[(df['Type'] == 'Insertion - In frame') | (df['Type'] == 'Deletion - In frame') | (df['Type'] == 'Complex - deletion inframe')]
all_missense = df [df ['Type'] == 'Substitution - Missense']
all_nonstop_extension =  df [df ['Type'] == 'Nonstop extension']

print("total stop codons: ", all_stop_codons ['Count'].sum()) #10
print("total frameshifts: ", all_frameshifts ['Count'].sum()) #71
print("total deletions/insertion in frame: ", all_deletionsInsertion_inframe ['Count'].sum()) #43
print("total missense: ", all_missense ['Count'].sum()) #5117
print("all_nonstop_extension: ", all_nonstop_extension  ['Count'].sum()) #5

#All mutations with catalytic site, defined to 761 aa
all_stop_codons_cat = df[(df['Position'] >= 761) & (df['Type'] == 'Substitution - Nonsense')] 
all_frameshifts_cat = df[(df['Position'] >= 761) & ((df ['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site: ", all_stop_codons_cat ['Count'].sum()) #0
print("total frameshifts w/cat site: ", all_frameshifts_cat ['Count'].sum()) #14
print("total deletions/insertion in frame w/cat site: ", all_deletionsInsertion_inframe ['Count'].sum()) #all (checked manually)

#Mutations influencing $789 with cat site:
all_stop_codons_cat_tail = df[(df ['Position'] >= 761) & (df['Position'] <= 782) & (df ['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat_tail = df[(df ['Position'] >= 761) & (df['Position'] <= 782) & ((df ['Type'] == 'Deletion - Frameshift') | (df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Complex - deletion inframe'))]

print("total stop codons w/cat site influencing S791: ", all_stop_codons_cat_tail['Count'].sum()) #0
print("total frameshifts w/cat site influencing S791: ", all_frameshifts_cat_tail['Count'].sum()) #4
print("Insertions/Deletions in frame with catalytic site influencing S777 and/or S789: 0") #manual count