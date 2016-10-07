# Written by Victoria Hsiao, August 2016
# echo_csv_maker creates a .csv file that can be directly read by the Echo liquid handler
# it also creates a .txt file that lists how much volume you need of each thing

# Inputs needed
# Integrase plasmid
# dna1_conc = 300 # ng/ul
# dna1_len = 4524 # bp

# Reporter plasmid
# dna2_conc = 300 # ng/ul
# dna2_len = 4656 # bp

# dna1_final = range(0,11,2) # in nM
# dna2_final = range(0,11,2) # in nM
# outputname = 'Bxb1' # string for naming output files

# Currently does not support matrices that are nonsymmetrical


def dna2nM_convert(dnaconc, dnalength):
    # Converts DNA ng/ul to nM
    return (dnaconc*1e6)/(660*dnalength)


def echo_csv_maker(rxn_vol, dna1_conc, dna1_len, dna2_conc, dna2_len, dna1_final, dna2_final, outputname, src_wells, firstwell):

    # This is code to calculate TXTL volumes and automatically generate an Echo picklist
    # Import numpy but use the nickname np.function instead of numpy.function
    import numpy as np
    import string

    # Define a custom function for rounding numbers to intervals of 25nL
    def myround(x, base = 25):
        return int(base *round(float(x)/base))

    final_volume = rxn_vol * 1e3 # Final reaction volume (nL)

    # Final volume of mastermix should be 75% of total volume
    # Buffer and extract are calculated separately so that it can do a final calculation on how many tubes you need.
    txtlBuffer = (5.0/9.0 * 0.75)*final_volume
    txtlExtract = (4.0/9.0 * 0.75)*final_volume
    txtlMM = txtlBuffer + txtlExtract

    # 28.50 needed, 30ul per extract tube
    # 35.63 needed, 37ul per buffer tube

    # matrix size
    mat_len = len(dna1_final)

    # DNA_conc(ng/ul) * 1e6 (uL/L) * (bp*mol)/660g * 1/dna_length(bp = DNA (nM)
    # Double stranded DNA is 660g/(bp*mol), single stranded DNA is 330g/(bp*mol)

    dna1_nM = dna2nM_convert(dna1_conc,dna1_len)
    dna2_nM = dna2nM_convert(dna2_conc,dna2_len)

    dna1_vol = np.zeros(mat_len)
    dna2_vol = np.zeros(mat_len)

    # print 'dna2_vol', dna2_vol

    # This calculates exact uL volumes for desired final conc (nM)
    # Then converts the volume to nL and rounds to the nearest interval of 25nL

    for i in range(0, mat_len):
        dna1_vol[i] = myround(dna1_final[i]*(final_volume/dna1_nM))
        dna2_vol[i] = myround(dna2_final[i]*(final_volume/dna2_nM))

    # print 'dna1_vol:', dna1_vol
    # print 'dna2_vol:', dna2_vol
    dna1_finalconc = []
    dna2_finalconc = []

    for each in dna1_vol:
        dna1_actualconc = ((each/final_volume) * dna1_nM)
        dna1_finalconc.append('%.2f' %dna1_actualconc)

    for each in dna2_vol:
        dna2_actualconc = (each/final_volume) * dna2_nM
        dna2_finalconc.append('%.2f' %dna2_actualconc)

    # print 'dna1_finalconc', dna1_finalconc

    # Create crossmatrices for dna1 and dna2 to calculate how much water to add to each well
    dna1_finalvol = np.tile(dna1_vol,[mat_len,1]).transpose()
    dna2_finalvol = np.tile(dna2_vol,[mat_len,1])

    # print 'dna1_finalvol:', dna1_finalvol
    # print 'dna2_finalvol:', dna2_finalvol

    txtlMM_finalvol = np.ones((mat_len,mat_len))*txtlMM
    dnatot_finalvol = dna1_finalvol+dna2_finalvol

    water_finalvol = np.ones((mat_len,mat_len))*final_volume - txtlMM_finalvol - dnatot_finalvol

    # 384 well plate naming
    # Rows are B - O to avoid edges
    wellalpha = string.ascii_uppercase[1:15]

    # Columns are 2 - 23 max to avoid edges
    wellnum = [str(x) for x in list(range(2, 24, 1))]

    offset_row = wellalpha.find(firstwell[0])
    offset_column = wellnum.index(firstwell[1])

    wellID = []
    for i in range(0, mat_len):
        for j in range(0, mat_len):
            wellID.append(wellalpha[i + offset_row] + wellnum[j + offset_column])

    # print wellID

    dna1_exp = np.reshape(dna1_finalvol, (1,len(wellID)))[0].tolist()
    dna2_exp = np.reshape(dna2_finalvol, (1,len(wellID)))[0].tolist()
    water_exp = np.reshape(water_finalvol, (1,len(wellID)))[0].tolist()
    txtl_exp = [txtlMM] * len(wellID)

    # Sample concentration info for the csv file
    concID1 = []
    concID2 = []

    for each in dna1_finalconc:
        concID1.append(str(each) + 'nM')
    for each in dna2_finalconc:
        concID2.append(str(each) + 'nM')

    concID1 = np.reshape(np.tile(concID1,[mat_len,1]).transpose(),[1,len(wellID)])[0].tolist()
    concID2 = np.reshape(np.tile(concID2,[mat_len,1]),[1,len(wellID)])[0].tolist()

    # print concID1
    # print concID2

    import csv

    # Identifiers for each column
    SP = ["Source[1]"]*len(wellID) # source plate
    SPtype = ["384PP_AQ_BP"]*len(wellID) # source plate type
    DPtype = ["Nunc_384_black_glassbottom"]*len(wellID)

    sampID1 = ['IntDNA']*len(wellID) # Source ID for integrase DNA
    sampID2 = ['ReporterDNA']*len(wellID) # Source ID for reporter DNA
    sampID3 = ['Water']*len(wellID) # Source ID for water
    # sampID4 = ['TXTL MM']*len(wellID) # Source ID for TXTL

    swID1 = [src_wells[0]]*len(wellID) # Source well for Integrase DNA
    swID2 = [src_wells[1]]*len(wellID) # Source well for Reporter DNA
    swID3 = [src_wells[2]]*len(wellID) # Source well for water
    # swID4 = ['A4']*len(wellID) # Source well for TXTL

    with open((outputname + '_EchoInput.csv'), 'w') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["Source Plate","Source Plate Type","Source Well","Sample ID","Destination Plate Name","Destination Well","Transfer Volume","Sample Comment"])

        export1 = zip(SP, SPtype, swID1, sampID1, DPtype, wellID, dna1_exp, concID1)
        for row in export1:
            writer.writerow(row)

        export2 = zip(SP, SPtype, swID2, sampID2, DPtype, wellID, dna2_exp, concID2)
        for row in export2:
            writer.writerow(row)

        export3 = zip(SP, SPtype, swID3, sampID3, DPtype, wellID, water_exp)
        for row in export3:
            writer.writerow(row)

        # export4 = zip(SP, SPtype, swID4, sampID4, DPtype, wellID, txtl_exp)
        # for row in export4:
        #     writer.writerow(row)

    # Total amounts needed for everything
    rxn_num = len(dna1_exp)

    dna1_tot = sum(dna1_exp)
    dna2_tot = sum(dna2_exp)
    water_tot = sum(water_exp) # nL

    txtlMM_tot = sum(txtl_exp)/1000 # needed rxn volume

    # 28.50 needed, 30ul per extract tube
    # 35.63 needed, 37ul per buffer tube

    # Calculating how much txtlMM is needed if 15ul dead volume is needed per well
    # 70ul max per well
    # 15ul min per well
    # 55ul usable volume per well

    import math
    txtlMM_wellsneeded = math.ceil(txtlMM_tot/55) # this needs to rounded up
    txtlMM_extra = txtlMM_wellsneeded*15
    txtlMM_finaltot = txtlMM_tot + txtlMM_extra

    txtl_extract_tot = (4.0/9.0)*txtlMM_finaltot
    txtl_buffer_tot = (5.0/9.0)*txtlMM_finaltot

    tubes_extract_tot = math.ceil((txtl_extract_tot)/28.50)
    tubes_buffer_tot = math.ceil(txtl_buffer_tot/35.63)

    with open((outputname+'_stocks.txt'), 'w') as text_file:
        text_file.write("Reactions: %s \n" % rxn_num)
        text_file.write("DNA1 total volume: %s uL + 15uL dead volume\n" % str(dna1_tot/1000))
        text_file.write("DNA1 conc: %s ng/ul\n" % str(dna1_conc))

        text_file.write("DNA2 total volume: %s uL + 15uL dead volume\n" % str(dna2_tot / 1000))
        text_file.write("DNA2 conc: %s ng/ul\n" % str(dna2_conc))

        text_file.write("Water total volume: %s uL + 15uL dead volume\n" % str(water_tot/1000))

        text_file.write("Extract total volume: %s uL\n" % str(txtl_extract_tot))
        text_file.write("Buffer total volume: %s uL\n" % str(txtl_buffer_tot))
        text_file.write("Tubes of extract needed: %s \n" % str(tubes_extract_tot))
        text_file.write("Tubes of buffer needed: %s \n" % str(tubes_buffer_tot))

        text_file.write("txtlMM source wells: %s \n" % str(txtlMM_wellsneeded))
        text_file.write("TXTL Echo volume: %s nL\n" % str(txtl_exp[0]))
    ## Todo
    # Change TXTL source type to SP

    return