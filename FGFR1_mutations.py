import pandas as pd

file_path = 'FGFR1_mutations_cosmic.csv'

# Read the csv file into a data frame
df = pd.read_csv(file_path)

# Filter out silent mutations] and unknown mutations
non_silent_mutations = df[(df['Type'] != 'Substitution - coding silent') & (df['Type'] != 'Unknown')]

print("All mutations that are not silent or unknown:", non_silent_mutations['Count'].sum()) #821

#Mutations influencing S789
stop_codon_before_789 = df[(df['Position'] <= 789) & (df['Type'] == 'Substitution - Nonsense')]
frameshift_before_789 = df[(df['Position'] <= 789) & ((df['Type'] == 'Insertion - Frameshift')| (df['Type'] == 'Deletion - Frameshift'))]
mis_mut_S777 = df[(df ['Position'] == 777) & (df['Type'] == 'Substitution - Missense')]
mis_mut_S789 = df[(df ['Position'] == 789) & (df['Type'] == 'Substitution - Missense')]

#Add count.sum because there could be more counts per mutation
print("Stop codons before 789: ", stop_codon_before_789['Count'].sum()) #29
print("Frameshift before 789: ", frameshift_before_789['Count'].sum())  #14
print("Pointmutations on S777:", mis_mut_S777 ['Count'].sum()) #0
print("Pointmutations on S789:", mis_mut_S789 ['Count'].sum()) #1
print("Insertions/Deletions in frame influencing S777 and/or S789: 0") #manual count

#All mutations
all_stop_codons = df[df['Type'] == 'Substitution - Nonsense']
all_frameshifts = df[(df['Type'] == 'Insertion - Frameshift') | (df['Type'] == 'Deletion - Frameshift')]
all_deletionsInsertion_inframe = df[(df['Type'] == 'Deletion - In frame') | (df['Type'] == "Insertion - In frame")]
all_missense = df[df['Type'] == 'Substitution - Missense']

print("total stop codons: ", all_stop_codons['Count'].sum()) #33
print("total frameshifts: ", all_frameshifts['Count'].sum()) #15
print("total deletions/insertion in frame: ", all_deletionsInsertion_inframe['Count'].sum()) #13
print("total missense: ", all_missense['Count'].sum()) #760

#All mutations with catalytic site, defined to 767 aa
all_stop_codons_cat = df[(df['Position'] >= 767) & (df['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat = df[(df['Position'] >= 767) & ((df['Type'] == 'Insertion - Frameshift')| (df['Type'] == 'Deletion - Frameshift'))] 

print("total stop codons w/cat site: ", all_stop_codons_cat['Count'].sum()) #4
print("total frameshifts w/cat site: ", all_frameshifts_cat['Count'].sum()) #1
print("total deletions/insertion in frame w/cat site: ", all_deletionsInsertion_inframe['Count'].sum()) #all 13 (checked manually)

#Mutations influencing S789 with cat site: 
all_stop_codons_cat_tail = df[(df['Position'] >= 767) & (df['Position'] <= 789) & (df['Type'] == 'Substitution - Nonsense')]
all_frameshifts_cat_tail = df[(df['Position'] >= 767) & (df['Position'] <= 789) & ((df['Type'] == 'Insertion - Frameshift')| (df['Type'] == 'Deletion - Frameshift'))] 

print("total stop codons w/cat site influencing S789: ", all_stop_codons_cat_tail['Count'].sum()) #0
print("total frameshifts w/cat site influencing S789:", all_frameshifts_cat_tail['Count'].sum()) #0
print("Insertions/Deletions in frame with catalytic site influencing S777 and/or S789: 0") #manual count
