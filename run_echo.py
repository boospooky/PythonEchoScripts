from echo_functions import *
dna1_final = range(0,6,1) # in nM
dna2_final = range(0,6,1) # in nM

rxn_vol = 5
# Integrase plasmid
dna1_conc = 300 # 1704 ng/ul stock
dna1_len = 4524 # bp

# Reporter plasmid
dna2_conc = 300 #564 ng/ul stock
dna2_len = 4656 # bp

# bxb1_srcwells = ['dna1 source well','dna2 source well','water source well']
bxb1_srcwells = ['C1','C2','C3']
firstwell = ['D', '2']

echo_csv_maker(rxn_vol, dna1_conc, dna1_len, dna2_conc, dna2_len, dna1_final, dna2_final, 'bxb1',bxb1_srcwells, firstwell)

# TP901
tp901_conc = 300 #1017 ng/ul stock
tp901_len = 4521 # bp

tp901r_conc = 300 #492 ng/ul stock
tp901r_len = 4538 # bp

tp901_srcwells = ['D1','D2','D3']

echo_csv_maker(rxn_vol, tp901_conc, tp901_len, tp901r_conc, tp901r_len, dna1_final, dna2_final, 'tp901',tp901_srcwells,firstwell)

## Todo
# optimize well placement based on what's been used before