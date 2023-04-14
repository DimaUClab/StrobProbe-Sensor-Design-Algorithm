# StrobProbe-Sensor-Design-Algorithm
StrobProbe Catalytic Homogenous DNA SERS Sensor Design Algorithm

Published by: Amanda C. Macke - Dima Group @ University of Cincinnati
Authored by: Steven Quarin, Amanda C. Macke, Lyndsay Kissell, Maria S. Kelly, Ashan Dayananda, Joseph Ungvary, George Stan, Ruxandra I. Dima and Pietro Strobbia
In Collaboration with the Strobbia, Dima and Stan Groups

Details of the applied analysis can be found in the original publication: (Link Coming)

########################################################################### 

# Welcome to the StrobProbe DNA Sensor Design Algorithm  
# Developed in Collaboration with the Strobbia, Dima & Stan Research Groups  
# University of Cincinnati  
## VERSION 10 

############################################################################ 

This python code was delevoped on a linux system but can be ran in any python environment. In linux, the code is most easily executed from a terminal with the following command: 

> python3 StrobProbe_2023.py 

###########################################################################

Programmed with Python 3.8 & access to python libraries included in Anaconda (numpy, pandas, seqfold & Biopython) 
From seqfold (https://pypi.org/project/seqfold/) calculates the ∆G for a hairpin at a certain T. 
From biopython (https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html) 
    you need to install the biopython library to allow parsing of DNA sequences and accept user input and complement the sequence 

Dependencies Used: 
numpy == 1.22.4 
pandas == 1.4.4 
seqfold == 0.7.15 
biopython == 1.81 

###########################################################################

Abbreviations Used: 
P: Probe - PH: PlaceHolder - T: Target – F: Fuel 

General Description Algorithm Logic (StrobeProbe_2023.py): 
This program is used to determine the sequences and thermodynamics for designing a catalytic homogeneous DNA sensor using SERS based on a provided target (T). There is an input file (Sensor_Parameters.csv) that contains the needed user input – all input information must be included. The placeholder (PH) is the complement of the input T strand. Toe Hold 1 (TH1) for the T to hybridize with PH is determined by truncating the beginning of the PH strand. The beginner fuel (F) is determined from the T sequence as it needs to be the complement of the PH. The final F strand has to complete a haripin secondary structure, so an input series of “T” bases are added to the end of the beginner F sequence. A series of complementary bases are added to form the neck of the hairpin when the hairpin is closed. Finally, Toe Hold 2 (TH2) for PH to hybridize with F is determined by truncating the opposite end of the beginner F sequence that complements the opposite end of PH from TH1. This is the initiating site for PH-F hybridization. Because this piece is already a part of the PH and F, it does not need to be extracted. The probe (P) is determined from the complement of the PH. The P also has to form a hairpin, so the complementary sequence of bases from the opposite end of the P is added to the end. The sensor molecule has to have enough space to interact with the nanostar without sterically hindering the formation of the hairpin, so a spacer of “A” bases are added to the end that attaches to the nanostar to complete the probe sequence. All sequence hybridization thermodynamic properties are calculated based on data from [Santa Lucia & Hicks Annual Review of Biophysics and Biomolecular Structure 2004] and are used as strand hybridization checks for formation from experimental values. SeqFold is used to check energetics associated with hairpin formation. 

###########################################################################
Input file (Sensor_Parameters.csv): 
- A file including the necessary information regarding the system and sensor trying to be built. The name and location is typically ran from the same location, but variables are provided in lines 16 (INPUT_FILE_NAME) & 17 (INPUT_FILE_LOC) to be changed for your convenience.  
Note: Current folder = './' - 1 folder back = '../' - Into a new location from the current folder = ‘./new_folder/’ 

- This should be a .csv file formatted in the same way as provided in the example (with default settings) that includes the following information in the following order: 
    
    Target_Name: This is required so that the output of a given design is saved. If you do not change the name of this, you will overwrite a previous file. Make sure to always use systematic and unique file names to prevent losing files.  Always backup your data! 

    Target: This is the sequence input you are trying to target. The input may be provided in upper OR lower case letters corresponding to the nucleic acids one would find in a DNA sequence. These letters will then be converted to lowercase letters for the analysis. 

    Salt_Correction: This is the concentration of your monovalent salt. It is dependent on the experimental buffer used (and or is the concentration of the salt used) 

    Temperature: This is the experimental temperature used. Keep in mind, room temperature is usually considered as 25 *C. The temperature used is extremely important for these experiments.  

    ΔGibbs_PHT_Maximum: The ΔGibbs Maximum for your range in kcal/mol for the PlaceHolder-Target hybrid 

    ΔGibbs_PHT_Minimum: The ΔGibbs Minimum for your range in kcal/mol for the PlaceHolder-Target hybrid 

    ToeHold_Minimum: This is used for the minimum length of Toe Hold 1 which is on the final PlaceHolder for the Target to grab onto. Toe Hold 2 has to be longer than Toe Hold 1  

    Δ Δ Gibbs_PHT_PPH: The change in Gibbs between the PlaceHolder-Target hybrid and the Probe-PlaceHolder hybrid. This is used as a thermodynamic check to evaluate if the placeholder will preferentially bind as expected 

    P_Hairpin_Minimum: This is the minimum length of the probe hairpin. For now, it is 3 which is the minimum number of nucleic acids required to form a hairpin structure. 

    Probe_Spacer: This is the sequence of A added to the nano-star end of the probe sequence 

    Fuel_Loop: Sequence given to create loop between neck and complement for proper hairpin formation  

    F_Hairpin_Minimum: Min length the fuel hairpin can be 

    Δ Δ Gibbs_FPH_TPH_Maximum: You want the Δ Δ Gibbs to be more negative than Δ Δ Gibbs_PHT_PPH from a previous variable input 

    Δ Δ Gibbs_FPH_TPH_Minimum: Needed to provide a targeted window. Experiments show that this value should be more negative to drive the reaction, but if it is too negative the PlaceHolder will not behave as needed. 
###########################################################################
Output File (ProbeDesign_{TARGET_NAME}.txt): 
- All of the pertinent sequence strand, hybrid, and hairpin information is stored in the output file. The sections are separated by which sequences are being built and evaluated. The appropriate sequences associated with the design parameters are displayed at the bottom. 

###########################################################################
Variable Descriptors: 
  Place Holder 1 (PH1) : Reverse Complement of Target Input  
  Place Holder 2 (PH2) : Reverse Complement of Target Input (NO Toe Hold 1)  
  Place Holder 3 (PH3) : Reverse Complement of Target Input (NO Toe Hold 1) + Toe Hold 2 Complement
  Place Holder FINAL (PH_final) : Reverse Complement of Target Input including Toe Hold 1 & Toe Hold 2 Complement  
  
  Probe 1 (P1) : Reverse Complement of Place Holder 2 (NO Toe Hold 1 & NO Toe Hold 2)  
  Probe 2 (P2) : Reverse Complement of Place Holder 3 (NO Toe Hold 1 WITH Toe Hold 2)  
  Probe FINAL (PROBE_final) : SPACER + Neck Complement + P2 (NO Toe Hold 1)  

  Fuel 1 (F1) : Target Sequence without last X elements corresponding to TH1  
  Fuel 2 (F2) : TH2 + Target Sequence (NO TH1)  
  Fuel 3 (F3) : F2 + FUEL LOOP input  
  Fuel FINAL (FUEL_final) : TH2 + F1 + FUEL_LOOP input + Neck Complement   

  Toe Hold 2 (TH2) : Includes the last X elements (randomly generated) of the Place Holder NOT included in the Target  
  Toe Hold 1 (TH1) : Includes the first X elements of the Place Holder NOT included in the Probe 
