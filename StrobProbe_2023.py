# Authored by Cullen McComb, Amanda C. Macke, Maria S. Kelly, Ashan Dayananda, Ruxandra I Dima, & Pietro Strobbia
# Last Edited - Amanda 04/14/23

import os
import sys
import time
from collections import namedtuple
import math 
import random

import numpy as np
import pandas as pd
from seqfold import dg, dg_cache, fold, Struct   # Assuming dg is calculated in kcal/mol
from Bio.Seq import Seq 

# Please find the Variable Descriptions in the READ.ME file
# Input File Location + Name
INPUT_FILE_NAME = 'Sensor_Parameters'
INPUT_FILE_LOC = './'

print('##################################################################################')
print('Welcome to the StrobProbe DNA Sensor Design Algorithm \nDeveloped in Collaboration with the Strobbia, Dima & Stan Research Groups \nUniversity of Cincinnati')
print('VERSION 10') 
print('##################################################################################\n')
# Equations and Values taken from SantaLucia&Hicks2004

###### Defining Required Functions #####
# Change in Gibbs Free Energy as a Function of Time
def GIBBS_CALC (T, H, S):
    g_t = H-(T+273.15)*S    # H & S are in kcal ; Convert temp to Kelvin
    return g_t

# Fixing the calculated Gibbs Free Energy to account for the length of the sequence and the salt
def GIBBS_FIXER (length):
    g_fix = 0.114*(length/2)*math.log(SALT_CORR)
    return(g_fix)

# Melting Temperature of a given pair of strands
def MELT_TEMP(H, S):
    t_m = (H*1000)/((S*1000)+(-0.0108)+(0.00199*math.log(0.25e-9)))+16.6*math.log(SALT_CORR)   # Converts the H to cal/mol and the S to cal/mol
    return(t_m)

# Calculate the DH, DS, & DG values for a given strand couple 
# Based on the thermodynamic nearest neighbor parameters for Watson-Crick base pairs in 1M NaCl
def STRAND_THERMO(STRAND1, STRAND2, H, S, G):  
    # Create base variables to use in the nearest neighbors approximation 
    # Set our Starting Point for these variables - Initial Energy
    i = 2
    h = 0.2   # enthalpy (h) in kcal/mol
    s = - 5.7    # entropy (s) in eu (cal/Kmol)
    g = 1.96   # Gibbs Free Energy (g) in kcal/mol
    
    #begin to adjust the variables based on the bases present
    if STRAND1[0] == 't' or STRAND1[0] == 'a':
        h = h + 2.2
        s = s + 6.9
        g = g + 0.05
    if STRAND1[len(STRAND1)-1] == 't' or STRAND1[len(STRAND1)-1] == 'a':
        h = h + 2.2
        s = s + 6.9
        g = g + 0.05
    while i <= len(STRAND1):  
        coup = STRAND2 [i-2:i]   
        #print(coup)
        if coup == 'aa' or coup == 'tt':
            h = h - 7.6
            s = s - 21.3
            g = g - 1
        if coup == 'at':
            h = h - 7.2
            s = s - 20.4
            g = g - 0.88
        if coup == 'ta':
            h = h - 7.2
            s = s - 21.3
            g = g - 0.58     
        if coup == 'ca' or coup == 'tg':
            h = h - 8.5
            s = s - 22.7
            g = g - 1.45    
        if coup == 'gt' or coup == 'ac':
            h = h - 8.4
            s = s - 22.4
            g = g - 1.44   
        if coup == 'ct' or coup == 'ag':
            h = h - 7.8
            s = s - 21.0
            g = g - 1.28
        if coup == 'ga' or coup == 'tc':
            h = h - 8.2
            s = s - 22.2
            g = g - 1.30
        if coup == 'cg':
            h = h - 10.6
            s = s - 27.2
            g = g - 2.17 
        if coup == 'gc':
            h = h - 9.8
            s = s - 24.4
            g = g - 2.24
        if coup == 'gg' or coup == 'cc':
            h = h - 8.0
            s = s - 19.9
            g = g - 1.84
        i=i+1
    s = s/1000   # Convert Entropy to (kcal/Kmol)
    thermo_list=namedtuple("thermo_list", [H, S, G])
    return thermo_list(h, s, g)

################################### User Input Information ###################################
### Input informtation is saved to .csv
PARAMETERS = pd.read_csv('{0}.csv'.format(INPUT_FILE_LOC+INPUT_FILE_NAME), header=None)

### All Information Needs to be written to a .txt file 
if os.path.exists('ProbeDesign_{0}.txt'.format(PARAMETERS.iloc[0,1])):
    os.remove('ProbeDesign_{0}.txt'.format(PARAMETERS.iloc[0,1]))
PROBE_INFO = open(r'ProbeDesign_{0}.txt'.format(PARAMETERS.iloc[0,1]), 'a')

PROBE_INFO.write('##################################################################################\n')
PROBE_INFO.write('Welcome to the StrobProbe DNA Sensor Design Algorithm \nDeveloped in Collaboration with the Strobbia, Dima & Stan Research Groups \nUniversity of Cincinnati\n')
PROBE_INFO.write('VERSION 10\n') 
PROBE_INFO.write('##################################################################################\n')

PROBE_INFO.write('------------------------------------- Probe Design {0} --------------------------------------\n\n'.format(PARAMETERS.iloc[0,1]))

# Import Target (T)
PAR_TARGET = '{0}'.format(PARAMETERS.iloc[1,1])
# All sequences must be given in lowercase letters
PAR_TARGET.lower()
TARGET = Seq(PAR_TARGET) # Sets target
print('----- Target Identified -----')
print('T= {}'.format(PAR_TARGET))

# Calculate No of Phosphates
PHOSPHATES_T = int(len(TARGET))
print('Length of Target = {}'.format(PHOSPHATES_T))

# Design Paramters from Input File
PROBE_INFO.write('Design Parameters: \n'.format(PARAMETERS.iloc[0,1]))
PROBE_INFO.write('Target Sequence: {} \n'.format(PAR_TARGET))
PROBE_INFO.write('Length of Target: {0} \n'.format(PHOSPHATES_T))
# Import Salt Correction - for Monovalent Salt
SALT_CORR = float(PARAMETERS.iloc[2,1])
PROBE_INFO.write('Salt Correction Used (M): {0} \n'.format(SALT_CORR))
# Import Termperature 
TEMPERATURE = float(PARAMETERS.iloc[3,1])
PROBE_INFO.write('Temperature Used (C): {0} \n'.format(TEMPERATURE))
# Import GIBBS Limits
GIBBS_PHT_MAX = float(PARAMETERS.iloc[4,1])
PROBE_INFO.write('Max Gibbs Limit Used for PH-T (kcal/mol): {0} \n'.format(GIBBS_PHT_MAX))
GIBBS_PHT_MIN = float(PARAMETERS.iloc[5,1])
PROBE_INFO.write('Min Gibbs Limit Used for PH-T (kcal/mol): {0} \n'.format(GIBBS_PHT_MIN))

PROBE_INFO.write('\n------------------------- PLACE HOLDER GENERATOR -------------------------\n') 
print('\n--- Analyzing Target with Place Holder 1 (including Toe Hold 1) ---')
PH1 = TARGET.reverse_complement()   # PH1 = Complement of Target (includes TH1)
print('PH= {}'.format(PH1))
# Determining Target - Place Holder 1 Thermodynamics 
h_TPH1, s_TPH1, g_TPH1=STRAND_THERMO(TARGET, PH1, 'h_TPH1', 's_TPH1', 'g_TPH1')
print('Target-PlaceHolder1 DG: {} (kcal/mol)'.format(round(g_TPH1,4)))

time.sleep(2)
# Save Target - PH1 to .txt
PROBE_INFO.write(' \n')
PROBE_INFO.write('--- Target - Place Holder 1 Thermodynamics ---\n')
PROBE_INFO.write('- Theoretical Thermodynamics based on Santa Lucia 2004 -\n')
PROBE_INFO.write('DH - Enthalpy (kcal/mol): {0} \n'.format(round(h_TPH1, 4)))
PROBE_INFO.write('DS - Entropy (kcal/Kmol): {0} \n'.format(round(s_TPH1, 4)))
PROBE_INFO.write('DG - Gibbs Free Energy (kcal/mol): {0} \n'.format(round(g_TPH1, 4)))
time.sleep(2)

# Calculating the Gibbs Free Energy to account for the Temperature, Length, & Salt Correction
Gt_TPH1 = GIBBS_CALC (TEMPERATURE, h_TPH1, s_TPH1)
Gfixer_TPH1 = GIBBS_FIXER (PHOSPHATES_T)
Gcorr_TPH1 = Gt_TPH1 - Gfixer_TPH1
time.sleep(2)
PROBE_INFO.write(' \n')
PROBE_INFO.write('- Corrected Values -\n')
PROBE_INFO.write('Temp & Salt corrected Gibbs Free Energy (kcal/mol): {0} \n'.format(round(Gcorr_TPH1, 4)))
time.sleep(2)

# Calculate the Melting Temperature
TM_TPH1 = (MELT_TEMP(h_TPH1, s_TPH1))-273.15
PROBE_INFO.write('Melting Temperature (C): {0} \n'.format(round(TM_TPH1, 4)))
time.sleep(2)

# CHECK POINT
print('\n--- CHECK POINT ---')
# Target - Placeholder
# Corrected Gibbs has to be between GIBBS_PHT_MIN & GIBBS_PHT according to current information
if Gcorr_TPH1 < GIBBS_PHT_MAX: 
    PROBE_INFO.write('ERROR: The resulting Gibbs is too large/negative - you can try using a shorter target sequence')
    sys.exit(' ERROR: The resulting Gibbs is too large/negative - you can try using a shorter target sequence \n Gibbs PH-P Limit = {0} \n Corrected Gibbs T-PH = {1}'.format(GIBBS_PHT_MAX, Gcorr_TPH1))
elif Gcorr_TPH1 > GIBBS_PHT_MIN: 
    PROBE_INFO.write('ERROR: The resulting Gibbs is too small/positive - you can try using a longer target sequence')
    sys.exit('ERROR: The resulting Gibbs is too small/positive - you can try using a longer target sequence \n Gibbs PH-P min Limit = {0} \n Corrected Gibbs T-PH = {1}'.format(GIBBS_PHT_MIN, Gcorr_TPH1))
else:
    PROBE_INFO.write('{0} kcal/mol > Corrected Gibbs Free Energy > {1} kcal/mol \n Calculated Target-Placeholder Hybridization = {2} kcal/mol \n'.format(GIBBS_PHT_MIN, GIBBS_PHT_MAX, round((Gcorr_TPH1),1)))
    print('{0} kcal/mol > Corrected Gibbs Free Energy > {1} kcal/mol \n Calculated Target-Placeholder Hybridization = {2} kcal/mol'.format(GIBBS_PHT_MIN, GIBBS_PHT_MAX, round((Gcorr_TPH1),1)))

#####################################   Probe 1   #####################################
PROBE_INFO.write('\n------------------------- PROBE 1 GENERATOR -------------------------\n')
print('\n--------------- Finding Toe Hold 1 ---------------')
TH1_min = int(PARAMETERS.iloc[6,1])     # Minimum Toe Hold 1 length
DDG_PHT_PPH = int(PARAMETERS.iloc[7,1])     # DDG Value
CHECK = True
while CHECK==True:
    print('\n----- No. of elements removed:', TH1_min, '-----')   # First iteration is 6
    PH2 = PH1[TH1_min:len(PH1)]     # Placeholder without Toe Holder 1 Sequence
    TH1 = PH1[0:TH1_min]     # Toe Hold 1 Sequence
    print("Toe Hold 1 CHECK: ", TH1)
    PROBE_CHECK=PH2.reverse_complement()
    PHOSPHATES_P=int(len(PROBE_CHECK))
    # Determining Place Holder 2 - Probe Check Thermodynamics
    h_PH2P, s_PH2P, g_PH2P=STRAND_THERMO(PH2, PROBE_CHECK, 'h_PH2P', 's_PH2P', 'g_PH2P')
    time.sleep(2)
    # Calculating the Gibbs Free Energy to account for the Temperature, Length, & Salt Correction
    Gt_PH2P = GIBBS_CALC (TEMPERATURE, h_PH2P, s_PH2P)
    Gfixer_PH2P = GIBBS_FIXER (PHOSPHATES_P)
    Gcorr_PH2P = Gt_PH2P - Gfixer_PH2P
    time.sleep(2)
    # CHECK POINT Probe - Placeholder
    # The PH2P has to be more positive than Gcorr_TPH1 + DDG_PHT_PPH
    # DDG_PHT_PPH = Gcorr_PH2P(-) - Gcorr_TPH1(-)
    if (Gcorr_PH2P - Gcorr_TPH1 > DDG_PHT_PPH):  
        PROBE_INFO.write('Corrected PH2P > Corrected TPH1 + DDG PHT-PPH\n DDG_PHT_PPH = {1} kcal/mol\nCalculated Probe 1 - Place Holder 2 = {0} kcal/mol \n'.format(round((Gcorr_PH2P),3), DDG_PHT_PPH))

        print('Calculated Probe 1 - Place Holder 2 = {0} kcal/mol'.format(round((Gcorr_TPH1 + DDG_PHT_PPH),3)))
        CHECK = False
    else:
        if TH1_min == (len(PH1)-1):
            PROBE_INFO.write('ERROR: A Probe for the indicated Target cannot be found')
            sys.exit('ERROR: A Probe for the indicated Target cannot be found')
        TH1_min = TH1_min+1

P1=PROBE_CHECK
print('\n----- Identified Toe Hold 1 -----')

print('\n--- Analyzing Place Holder 2 with Probe 1 (NOT including Toe Hold 1) ---\n')
time.sleep(2)
# Save PH2 -P1 to .txt
PROBE_INFO.write(' \n')
PROBE_INFO.write('--- Place Holder 2 - Probe 1 Thermodynamics ---\n')
PROBE_INFO.write('- Toe Hold 1 ({0} elements): {1}'.format(TH1_min, TH1))
PROBE_INFO.write('- Theoretical Thermodynamics based on Santa Lucia 2004 -\n')
PROBE_INFO.write('DH - Enthalpy (kcal/mol): {0} \n'.format(round(h_PH2P, 4)))
PROBE_INFO.write('DS - Entropy (kcal/Kmol): {0} \n'.format(round(s_PH2P, 4)))
PROBE_INFO.write('DG - Gibbs Free Energy (kcal/mol): {0} \n'.format(round(g_PH2P, 4)))
PROBE_INFO.write(' \n')
PROBE_INFO.write('- Calculated Values -\n')
PROBE_INFO.write('Temp & Salt corrected Gibbs Free Energy (kcal/mol): {0} \n'.format(round(Gcorr_PH2P, 4)))
# Calculate the Melting Temperature
TM_PH2P = (MELT_TEMP(h_PH2P, s_PH2P))-273.15
PROBE_INFO.write('Melting Temperature (C): {0} \n'.format(round(TM_PH2P, 4)))
time.sleep(2)

#####################################   FUEL GENERATOR   #####################################
# Set initial Fuel strand as the Target minus TH1
PROBE_INFO.write('\n ------------------------- FUEL GENERATOR -------------------------\n')
F1 = TARGET[:-len(TH1)]
F1 = F1.lower()
TH2_min = len(TH1)     # Toe Hold 2 has to be larger than Toe Hold 1
# DDG has to be within a specified threshold
DDG_FPH_TPH_max = int(PARAMETERS.iloc[12,1])     # DDG maximum
DDG_FPH_TPH_min = int(PARAMETERS.iloc[13,1])     # DDG minimum
# Set conditional variables
gen_fuel_attempt = 1
fuel_found = False

DUMMY_COUNTER=TARGET
max_attempts=15
# Gives 5 chances before failing
while gen_fuel_attempt <= max_attempts:
    time.sleep(2)
    print('\n--- Attempt {} at identifying a fuel strand ---'.format(gen_fuel_attempt))

    # Start with empty TH2
    TH2 = str()
    if fuel_found==False:
    
        # Add nucleotides until len(TARGET)-2 and check
        while len(TH2) < len(DUMMY_COUNTER)-2:
            # Add random nucleotide to TH2 and find complement
            TH2_ADD = ' '.join([random.choice('tgca') for x in range(1)])
            TH2 = TH2_ADD + TH2
            TH2_COMP = Seq(TH2).reverse_complement()
            
            if len(TH2) <= TH2_min:   # Check that TH2 is at least TH_minimum elements long
                continue
            
            print("Toe Hold 2 CHECK ({0} elements): ". format(len(TH2)), TH2)
            # Construct F2 and PH3 using generated TH2
            F2 = TH2 + F1
            PH3 = PH2 + TH2_COMP   # TH1 is not included in the Place Holder-Fuel hybridization - PH3 is PH2 + TH2 complement

            # CHECK POINT
            # FUEL STRUCTURE CHECK: Looking for potential secondary structures in build strand
            Ghp_Ftest = dg(F2, temp=TEMPERATURE)     # Gibbs Free Energy of the Probe without hairpin
            Ghpcorr_Ftest = Ghp_Ftest + .0000015*(1000-(SALT_CORR*1000))**2    
            if (Ghpcorr_Ftest <= -4): 
                print('\n~~~ WARNING: Unwanted Secondary Structure in Fuel ~~~ \n')
        
                PROBE_INFO.write(' \n')
                PROBE_INFO.write('\n~~~~~ WARNING CHECK SECONDARY STRUCTURE IN FUEL ~~~~~\n')
                PROBE_INFO.write(' \n')
            else:
                print('\n~~~ No Unwanted Secondary Structure in Fuel Detected ~~~\n') 

            # Add Fuel Loop
            FUEL_LOOP = Seq('{0}'.format(PARAMETERS.iloc[10,1]))
            
            F3 = F2 + FUEL_LOOP
            
            # CHECK POINT
            # FUEL HAIRPIN CHECK: Find Neck of Fuel and run through check
            print('Fuel Hairpin Check') 
            time.sleep(2)
            k = int(PARAMETERS.iloc[8,1])     # Minimum No. of elements in the Hairpin
            CHECK = True
            while CHECK==True:
                NECK_F = F2[len(F2)-(k):len(F2)]    # 3' Segment of Fuel Identified for Neck
                NECK_F_COMP = NECK_F.reverse_complement()

                # Fuel including complement of neck
                FUEL_CHECK = F3 + NECK_F_COMP      

                # Use Seq Fold to determine if the hairpin will fold at TEMPERATURE
                Ghp_F = dg(FUEL_CHECK, temp=TEMPERATURE)     # Gibbs Free Energy of the Folded Probe
                Ghpcorr_F = Ghp_F + .0000015*(1000-(SALT_CORR*1000))**2

                if (Ghpcorr_F >= -6) == True and (Ghpcorr_F <= -2) == True: 
                    print('Fuel Hairpin Loop length: {0}'.format(k))
                    print('Fuel Hairpin Neck length: {0}'.format(NECK_F))
                    print('Hairpin passes Gibbs check = ', Ghpcorr_F, '(kcal/mol)')

                    # Check if the hairpin folds & unfolds
                    # Calculate the TM of the hairpin by determining the temperature when DG from Seq Fold is 0
                    # Starting from the calculated Ghpcorr_F
                    TEMP_CHECK = TEMPERATURE + 2     # Starting with our current temp
                    while CHECK == True:
                        Gcheck_HP = dg(FUEL_CHECK, temp=TEMP_CHECK)
                        Gcheckcorr_HP = Gcheck_HP + .0000015*(1000-(SALT_CORR*1000))**2

                        if (Gcheckcorr_HP < 0) == True:
                            # Check to see if the hairpin will open at 50 C (AT MELTING TEMP?)
                            G50_HP = dg(FUEL_CHECK, temp=TEMP_CHECK-1)
                            G50corr_HP = G50_HP + .0000015*(1000-(SALT_CORR*1000))**2
                            if G50corr_HP <= 1 and G50corr_HP >= -1:
                                print('The fuel hairpin will open at {0} C in {1} M monovalent salt.\n'.format(TEMP_CHECK-1, SALT_CORR))
                                CHECK=False
                            else:
                                TEMP_CHECK = TEMP_CHECK + 2
                                continue  
                        else:
                            TEMP_CHECK = TEMP_CHECK + 2  
                            if TEMP_CHECK > 50:
                                PROBE_INFO.write('This hairpin will not unfold - try a larger hairpin loop')
                                time.sleep(2)
                                break
                            else:
                                continue
                else:
                    k = k + 1 #Add to k
                    if k == ((len(FUEL_CHECK)/2)-2):
                        break
                    else:
                        continue
                
            # FUEL-PLACEHOLDER HYBRIDIZATION CHECK: Checking if TH2 is acceptable
            #Determining Fuel - Place Holder 3 Thermodynamics
            #h_FPH3, s_FPH3, g_FPH3=STRAND_THERMO(FUEL_CHECK, PH3, 'h_FPH3', 's_FPH3', 'g_FPH3')
            h_FPH3, s_FPH3, g_FPH3=STRAND_THERMO(F2, PH3, 'h_FPH3', 's_FPH3', 'g_FPH3')
            time.sleep(2)

            # Calculating the Gibbs Free Energy to account for the Temperature, Length, & Salt Correction
            Gt_FPH3 = GIBBS_CALC(TEMPERATURE, h_FPH3, s_FPH3)
            Gfixer_FPH3 = GIBBS_FIXER(PHOSPHATES_T)
            Gcorr_FPH3 = Gt_FPH3 - Gfixer_FPH3
            time.sleep(2)
            # Calculate the Melting Temperature
            TM_FPH3 = (MELT_TEMP(h_FPH3, s_FPH3))-273.15
            time.sleep(2)
            
            # CHECK POINT
            print('--- Check Point ---')          
            CHECK = True  
            ddg_FPH_TPH = Gcorr_FPH3 - Gcorr_TPH1
            print(FUEL_CHECK)
            print('DD Gibbs FPH - TPH = {0} kcal/mol'.format(round(ddg_FPH_TPH,1)))
            if DDG_FPH_TPH_min < ddg_FPH_TPH and  ddg_FPH_TPH < DDG_FPH_TPH_max:
                print('\n----- Fuel Identified -----')
                print('{0} kcal/mol < {1} kcal/mol < {2} kcal/mol'.format(DDG_FPH_TPH_min, round(ddg_FPH_TPH,3), DDG_FPH_TPH_max))
                PROBE_INFO.write('{1} kcal/mol < DDG_FPH3_TPH1 < {2} kcal/mol \nCalculated DDG Fuel - Place Holder 3 & Target - Place Holder 1 = {0} kcal/mol \n'.format(round(ddg_FPH_TPH,3), DDG_FPH_TPH_min, DDG_FPH_TPH_max))
                
                time.sleep(2)
                
                # Calculate the Melting Temperature
                TM_FPH3 = (MELT_TEMP(h_FPH3, s_FPH3))-273.15
                time.sleep(2)
                
                PROBE_INFO.write('\nHairpin Check:\n')
                #PROBE_INFO.write('The fuel hairpin will open at {0} C in {1} M monovalent salt\n'.format(TEMP_CHECK-1, SALT_CORR))
                PROBE_INFO.write('Fuel Hairpin LOOP of {0} elements\n'.format(len(FUEL_LOOP)))
                PROBE_INFO.write('Fuel Hairpin NECK of {0} elements\n'.format(len(NECK_F_COMP)))

                FUEL_final = FUEL_CHECK
                print(FUEL_final)
                PH_final = TH1 + PH3   # Add TH1 back onto PH sequence that also includes the TH2
                
                # Need the following lines to exit out of loops
                fuel_found = True
                DUMMY_COUNTER = [range(0, len(DUMMY_COUNTER)+1)]
                gen_fuel_attempt = max_attempts+1
                
                # Save Fuel to .txt
                PROBE_INFO.write(' \n')
                PROBE_INFO.write('--- Fuel 2 - Place Holder 3 Thermodynamics ---\n')
                PROBE_INFO.write('Toe Hold 2 ({0} elements): {1} \n'.format(len(TH2), TH2))
                PROBE_INFO.write(' \n')
                PROBE_INFO.write('- Theoretical Thermodynamics based on Santa Lucia 2004 -\n')
                PROBE_INFO.write('DH - Enthalpy (kcal/mol): {0} \n'.format(round(h_FPH3, 4)))
                PROBE_INFO.write('DS - Entropy (kcal/Kmol): {0} \n'.format(round(s_FPH3, 4)))
                PROBE_INFO.write('DG - Gibbs Free Energy (kcal/mol): {0} \n'.format(round(g_FPH3, 4)))
                PROBE_INFO.write(' \n')
                PROBE_INFO.write('- Calculated Values -\n')
                PROBE_INFO.write('Temp & Salt corrected Gibbs Free Energy (kcal/mol): {0} \n'.format(round(Gcorr_FPH3, 4)))
                PROBE_INFO.write('Melting Temperature (C): {0} \n'.format(round(TM_FPH3, 4)))
                time.sleep(2)
                break
            
            else:
                if len(TH2) == (len(F2)-2):
                    print(' \n')
                    print('--- ERROR --- \n'.format(gen_fuel_attempt))
                    print('ERROR: Toe Hold 2 is too long \n')
                    print('Trying Next Attempt')
                    time.sleep(2)
                    gen_fuel_attempt += 1
                    break
                if ddg_FPH_TPH < DDG_FPH_TPH_min:
                    print(' \n')
                    print('--- ERROR --- \n'.format(gen_fuel_attempt))
                    print('ERROR: DD Gibbs FPH - TPH is too negative \n')
                    print('Trying Next Attempt')
                    time.sleep(2)
                    gen_fuel_attempt += 1
                    if gen_fuel_attempt == max_attempts+1:   #Ends while loop
                        PROBE_INFO.write('ERROR: Could Not identify a Toe Hold 2 for Fuel Strand after 5 attempts')
                        time.sleep(2)
                        sys.exit('\nERROR: Could Not identify a Toe Hold 2 for Fuel Strand after 5 attempts')
                    break
            if gen_fuel_attempt == max_attempts+1:   #Ends while loop
                PROBE_INFO.write('ERROR: Could Not identify a Toe Hold 2 for Fuel Strand after {0} attempts'.format(max_attempts))
                time.sleep(2)
                sys.exit('\nERROR: Could Not identify a Toe Hold 2 for Fuel Strand after {0} attempts'.format(max_attempts))
        time.sleep(2)
    
    # When fuel_found=True, you passed both checks
time.sleep(2)
print('\n----- Identified Toe Hold 2 -----')
print("Final FUEL Strand Identified:")
print(FUEL_final)
time.sleep(2)
print("Final PLACE HOLDER Strand Identified:")
print(PH_final)
time.sleep(2)

#####################################   PROBE GENERATOR   #####################################
# Spacer - TH2* - PH* - Hairpin Loop - Neck Complement
PROBE_INFO.write('\n------------------------- PROBE HAIRPIN GENERATOR -------------------------\n')
# Add TH2 complement onto beginning of P1 - Complement of PH3 
P2 = PH3.reverse_complement()   # PH3 = PH2 + TH2 complement (No TH1)

# CHECK POINT
# FUEL HAIRPIN CHECK: Find Neck of Fuel and run through check
PROBE_INFO.write('Hairpin Check:')
# Generating the Neck to add to the 5' end to form the Probe Hairpin 
k=int(PARAMETERS.iloc[8,1])     # Minimum No. of elements in the Hairpin

CHECK = True
print('\nProbe Hairpin Check:') 


Ghp_Ptest = dg(P2, temp=TEMPERATURE)     # Gibbs Free Energy of the Probe without hairpin
Ghpcorr_Ptest = Ghp_Ptest + .0000015*(1000-(SALT_CORR*1000))**2
if (Ghpcorr_Ptest <= -3): 
    print('\n~~~ WARNING: Unwanted Secondary Structure in Fuel ~~~ \n')
    # 
    PROBE_INFO.write(' \n')
    PROBE_INFO.write('\n~~~~~ WARNING CHECK SECONDARY STRUCTURE IN FUEL ~~~~~\n')
    PROBE_INFO.write(' \n')            
else:
    print('\n~~~ No Unwanted Secondary Structure in Fuel Detected ~~\n') 


while CHECK==True:
    NECK_P = P2[len(P2)-(k):len(P2)]    # Extra portion of probe to make the hairpin loop
    NECK_P_comp = NECK_P.reverse_complement()     # Add the portion to end of probe for the neck of the hairpin
    PROBE_CHECK= NECK_P_comp + P2     # Complete Probe Strand to check (NECK COMP + P1 including NECK)
    print('Checking Hairpin: {}'.format(PROBE_CHECK))
    # CHECK POINT
    # Use Seq Fold to determine if the hairpin will fold at TEMPERATURE
    Ghp_P = dg(PROBE_CHECK, temp=TEMPERATURE)     # Gibbs Free Energy of the Folded Probe
    Ghpcorr_P = Ghp_P + .0000015*(1000-(SALT_CORR*1000))**2
    print('\nChecking Corrected Hairpin DG: {} (kcal/mol)'.format(Ghpcorr_P))
    if (Ghpcorr_P >= -7) == CHECK and (Ghpcorr_P <= -5) == CHECK:
        print('\nProbe Hairpin Loop length: {0}'.format(k)) 
        print('Probe Hairpin Neck: {0}'.format(NECK_P)) 
        print('Hairpin passes Gibbs check: ', Ghpcorr_P, '(kcal/mol)')
        
        SPACER = '{0}'.format(PARAMETERS.iloc[9,1])
        PROBE_test = SPACER + PROBE_CHECK
        # Check if the hairpin folds & unfolds
        # Calculate the TM of the hairpin by determining the temperature when DG from Seq Fold is 0
        # Starting from the calculated Ghpcorr_P
        TEMP_CHECK = TEMPERATURE + 2     # Starting with our current temp
        while CHECK == True:
            #print('Temp Iteration', TEMP_CHECK)
            Gcheck_HP = dg(PROBE_test, temp=TEMP_CHECK)
            Gcheckcorr_HP = Gcheck_HP + .0000015*(1000-(SALT_CORR*1000))**2
            if (Gcheckcorr_HP < 0) == CHECK:
                # Check to see if the hairpin will open at 50 C (AT MELTING TEMP?)
                G50_HP = dg(PROBE_test, temp=TEMP_CHECK-1)   # -1 for the temperature before it melts
                G50corr_HP = G50_HP + .0000015*(1000-(SALT_CORR*1000))**2
                if G50corr_HP <= 1 and G50corr_HP >= -1:
                    print('The probe hairpin will open at', TEMP_CHECK-1, 'C in', SALT_CORR, 'M monovalent salt.')
                    time.sleep(2)
                    CHECK=False
                else:
                    TEMP_CHECK = TEMP_CHECK + 2
                    time.sleep(2)
                    continue  
            else:
                TEMP_CHECK = TEMP_CHECK + 2  
                if TEMP_CHECK > 50:
                    print('This hairpin will not unfold - try a larger haripin loop')
                    break
        #CHECK=False
    if k == ((len(PROBE_CHECK)/2)-2):
        PROBE_INFO.write('ERROR: The Hairpin does not work - you need a different Target')
        sys.exit('ERROR: The Hairpin does not work - you need a different Target')
    k = k + 1
print('\n --- Hairpin Generated ---')
PROBE_final = PROBE_test
# Save Hairpin to .txt
time.sleep(2)
PROBE_INFO.write(' \n')
PROBE_INFO.write('-- Hairpin Thermodynamics with Seq Fold --\n')
PROBE_INFO.write('Probe Hairpin LOOP of {0} elements\n'.format(k))
PROBE_INFO.write('Probe Hairpin NECK of {0} elements\n'.format(len(NECK_P_comp)))
PROBE_INFO.write('The probe hairpin will open at {0} C in {1} M monovalent salt\n'.format(TEMP_CHECK-1,SALT_CORR))
PROBE_INFO.write('- Calculated Values -\n')
PROBE_INFO.write('Temp & Salt corrected Gibbs Free Energy (kcal/mol): {0} \n'.format(round(Ghpcorr_P, 4)))
PROBE_INFO.write('Melting Temp - Hairpin will open (C): {0} \n'.format(TEMP_CHECK-1))
PROBE_INFO.write(' \n')
time.sleep(2)
print("Final PROBE Strand Identified:")
print(PROBE_final)

### Save Final Sequences to .txt ###
PROBE_INFO.write('##################################################################################\n')
PROBE_INFO.write('------------------------------ FINAL SEQUENCES ------------------------------\n')
PROBE_INFO.write('Identified Target:\n')
PROBE_INFO.write('{0} \n'.format(PAR_TARGET))
PROBE_INFO.write('Length of Target: {0}\n'.format(len(PAR_TARGET)))
PROBE_INFO.write(' \n')
time.sleep(2)
PROBE_INFO.write('Identified Final Place Holder: \n')
PROBE_INFO.write('{0} \n'.format(PH_final))
PROBE_INFO.write('Length of Final Place Holder: {0}\n'.format(len(PH_final)))
PROBE_INFO.write(' \n')
time.sleep(2)
PROBE_INFO.write('Identified Final Fuel : \n')
PROBE_INFO.write('{0} \n'.format(FUEL_final))
PROBE_INFO.write('Length of Final Fuel: {0}\n'.format(len(FUEL_final)))
PROBE_INFO.write(' \n')
time.sleep(2)
PROBE_INFO.write('Identified Final Probe: \n')
PROBE_INFO.write('{0} \n'.format(PROBE_final))
PROBE_INFO.write('Length of Final Probe: {0}\n'.format(len(PROBE_final)))
PROBE_INFO.write(' \n')
time.sleep(2)
PROBE_INFO.write('Identified Final Toe Hold 1 ({0} elements): \n'.format(len(TH1)))
PROBE_INFO.write('{0} \n'.format(TH1))
PROBE_INFO.write(' \n')
time.sleep(2)
PROBE_INFO.write('Identified Final Toe Hold 2 ({0} elements): \n'.format(len(TH2)))
PROBE_INFO.write('{0} \n'.format(TH2))
PROBE_INFO.write('\nThanks for using StrobeProbe.py for your DNA Sensor Design Needs!\n')
PROBE_INFO.close() 
time.sleep(2)
print('\n--------------------------------------------------------------------------------------------')
print('All pertinent information about the Sensor Design has been exported to the ProbeDesign_{Target_Name}.txt file in the current directory')
print('--------------------------------------------------------------------------------------------')

print('\nThanks for using StrobProbe for your DNA Sensor Design Needs!')
print('~~~ Done :-) ~~~')
