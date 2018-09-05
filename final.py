# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 09:43:17 2017

@author: HP
"""

#program to create dataset(protein_features_refiningresult.csv) for  protein collected from databank of pdb files(protein_repository) 
#complexity O(n.m) 



#import library 
from pydpi.protein import getpdb,AAComposition #for ammino acid sequence and others
from pydpi.pypro import CTD #compostion transition distribution descriptors
from pydpi.pypro import PyPro  #for pseudo amino acid composition descriptors
#import prody as pdy
#from pylab import *
import time
import os
import glob
import shutil
import sys
import multiprocessing
#arguments
args=sys.argv
if len(args)!=2:
    print "arguments typed:",args
    print "\nError !!! Wrong number of parameters"
    print "\nUsages: $python final.py <sr_no>" #sno should start from 0
    exit()  
    
#parameters
filename="protein_dataset"
path1 =os.getcwd()+"\\"+"protein_repository"
path2 =os.getcwd()
pdb_dic={}
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X'}
sno=int(args[1])

startTime = time.time()
#writing result in csv file
fp = open(str(sno)+filename+".csv",'w')
#fp.write("protien_name,rmsd,tmscore,gdt_ts,amino_acid_sequence_length,A,C,E,D,G,F,I,H,K,M,L,N,Q,P,S,R,T,W,V,Y,_NormalizedVDWVD1075,_PolarityD1075,_SecondaryStrD3025,_PolarityD3100,_ChargeD1100,_SecondaryStrT23,_PolarityD3025,_NormalizedVDWVC1,_NormalizedVDWVC3,_HydrophobicityT23,_SolventAccessibilityD1025,_PolarityD1100,_NormalizedVDWVD2100,_ChargeD3050,_PolarityD2001,_SolventAccessibilityD2025,_SecondaryStrD2025,_PolarityT12,_PolarityT13,_ChargeT23,_HydrophobicityD2100,_SolventAccessibilityD3001,_ChargeD1075,_SecondaryStrD1075,_PolarizabilityD2050,_SolventAccessibilityD3100,_PolarizabilityD1075,_HydrophobicityD3025,_PolarizabilityD3025,_SecondaryStrD3050,_SolventAccessibilityD2050,_HydrophobicityD2001,_SecondaryStrD2075,_PolarityD3050,_PolarityD1001,_ChargeD2075,_NormalizedVDWVD3025,_SolventAccessibilityD1100,_SecondaryStrD1100,_NormalizedVDWVD1001,_NormalizedVDWVD2050,_NormalizedVDWVC2,_PolarizabilityC2,_PolarizabilityC3,_PolarizabilityC1,_SecondaryStrC2,_SecondaryStrC3,_ChargeT12,_ChargeT13,_ChargeD1001,_NormalizedVDWVD3001,_NormalizedVDWVD2025,_NormalizedVDWVD3050,_SecondaryStrD1025,_SolventAccessibilityD1050,_PolarizabilityD3001,_PolarityD2075,_SecondaryStrD3001,_ChargeD2025,_HydrophobicityD2050,_PolarizabilityD1100,_HydrophobicityC1,_HydrophobicityC2,_SolventAccessibilityD2100,_PolarizabilityD2001,_PolarizabilityD1025,_NormalizedVDWVT13,_PolarizabilityD3050,_SolventAccessibilityT13,_SolventAccessibilityT12,_PolarityD1050,_ChargeD3075,_SolventAccessibilityD1001,_HydrophobicityC3,_HydrophobicityD3100,_SecondaryStrD1050,_ChargeD3001,_SecondaryStrD2050,_PolarityD2025,_NormalizedVDWVD2001,_HydrophobicityD3075,_ChargeD1025,_SolventAccessibilityD2001,_HydrophobicityD1025,_ChargeD1050,_SolventAccessibilityD3025,_PolarizabilityD3075,_PolarityD2100,_SecondaryStrD3100,_PolarizabilityT13,_PolarizabilityT12,_SecondaryStrD2001,_NormalizedVDWVD1050,_PolarizabilityD1050,_SolventAccessibilityT23,_ChargeC3,_ChargeC2,_ChargeC1,_NormalizedVDWVT23,_NormalizedVDWVD3100,_PolarizabilityD2025,_SolventAccessibilityD3075,_HydrophobicityD2075,_HydrophobicityD1075,_ChargeD2050,_ChargeD3100,_SecondaryStrC1,_ChargeD3025,_PolarityD3001,_PolarityD3075,_PolarityD1025,_SecondaryStrT13,_SecondaryStrT12,_HydrophobicityD1001,_SolventAccessibilityD1075,_HydrophobicityD1050,_HydrophobicityT13,_HydrophobicityT12,_SecondaryStrD1001,_NormalizedVDWVD1100,_NormalizedVDWVD3075,_NormalizedVDWVD1025,_SolventAccessibilityD2075,_SecondaryStrD2100,_PolarizabilityD3100,_SecondaryStrD3075,_NormalizedVDWVT12,_PolarizabilityT23,_SolventAccessibilityD3050,_PolarityT23,_PolarityD2050,_HydrophobicityD2025,_PolarizabilityD2075,_HydrophobicityD3050,_HydrophobicityD3001,_HydrophobicityD1100,_ChargeD2001,_SolventAccessibilityC1,_SolventAccessibilityC2,_SolventAccessibilityC3,_PolarizabilityD1001,_NormalizedVDWVD2075,_ChargeD2100,_PolarizabilityD2100,_PolarityC1,_PolarityC3,_PolarityC2\n")
fp.write("protien_name,rmsd,tmscore,gdt_ts,amino_acid_sequence_length,A,C,E,D,G,F,I,H,K,M,L,N,Q,P,S,R,T,W,V,Y,_NormalizedVDWVD1075,_PolarityD1075,_SecondaryStrD3025,_PolarityD3100,_ChargeD1100,_SecondaryStrT23,_PolarityD3025,_NormalizedVDWVC1,_NormalizedVDWVC3,_HydrophobicityT23,_SolventAccessibilityD1025,_PolarityD1100,_NormalizedVDWVD2100,_ChargeD3050,_PolarityD2001,_SolventAccessibilityD2025,_SecondaryStrD2025,_PolarityT12,_PolarityT13,_ChargeT23,_HydrophobicityD2100,_SolventAccessibilityD3001,_ChargeD1075,_SecondaryStrD1075,_PolarizabilityD2050,_SolventAccessibilityD3100,_PolarizabilityD1075,_HydrophobicityD3025,_PolarizabilityD3025,_SecondaryStrD3050,_SolventAccessibilityD2050,_HydrophobicityD2001,_SecondaryStrD2075,_PolarityD3050,_PolarityD1001,_ChargeD2075,_NormalizedVDWVD3025,_SolventAccessibilityD1100,_SecondaryStrD1100,_NormalizedVDWVD1001,_NormalizedVDWVD2050,_NormalizedVDWVC2,_PolarizabilityC2,_PolarizabilityC3,_PolarizabilityC1,_SecondaryStrC2,_SecondaryStrC3,_ChargeT12,_ChargeT13,_ChargeD1001,_NormalizedVDWVD3001,_NormalizedVDWVD2025,_NormalizedVDWVD3050,_SecondaryStrD1025,_SolventAccessibilityD1050,_PolarizabilityD3001,_PolarityD2075,_SecondaryStrD3001,_ChargeD2025,_HydrophobicityD2050,_PolarizabilityD1100,_HydrophobicityC1,_HydrophobicityC2,_SolventAccessibilityD2100,_PolarizabilityD2001,_PolarizabilityD1025,_NormalizedVDWVT13,_PolarizabilityD3050,_SolventAccessibilityT13,_SolventAccessibilityT12,_PolarityD1050,_ChargeD3075,_SolventAccessibilityD1001,_HydrophobicityC3,_HydrophobicityD3100,_SecondaryStrD1050,_ChargeD3001,_SecondaryStrD2050,_PolarityD2025,_NormalizedVDWVD2001,_HydrophobicityD3075,_ChargeD1025,_SolventAccessibilityD2001,_HydrophobicityD1025,_ChargeD1050,_SolventAccessibilityD3025,_PolarizabilityD3075,_PolarityD2100,_SecondaryStrD3100,_PolarizabilityT13,_PolarizabilityT12,_SecondaryStrD2001,_NormalizedVDWVD1050,_PolarizabilityD1050,_SolventAccessibilityT23,_ChargeC3,_ChargeC2,_ChargeC1,_NormalizedVDWVT23,_NormalizedVDWVD3100,_PolarizabilityD2025,_SolventAccessibilityD3075,_HydrophobicityD2075,_HydrophobicityD1075,_ChargeD2050,_ChargeD3100,_SecondaryStrC1,_ChargeD3025,_PolarityD3001,_PolarityD3075,_PolarityD1025,_SecondaryStrT13,_SecondaryStrT12,_HydrophobicityD1001,_SolventAccessibilityD1075,_HydrophobicityD1050,_HydrophobicityT13,_HydrophobicityT12,_SecondaryStrD1001,_NormalizedVDWVD1100,_NormalizedVDWVD3075,_NormalizedVDWVD1025,_SolventAccessibilityD2075,_SecondaryStrD2100,_PolarizabilityD3100,_SecondaryStrD3075,_NormalizedVDWVT12,_PolarizabilityT23,_SolventAccessibilityD3050,_PolarityT23,_PolarityD2050,_HydrophobicityD2025,_PolarizabilityD2075,_HydrophobicityD3050,_HydrophobicityD3001,_HydrophobicityD1100,_ChargeD2001,_SolventAccessibilityC1,_SolventAccessibilityC2,_SolventAccessibilityC3,_PolarizabilityD1001,_NormalizedVDWVD2075,_ChargeD2100,_PolarizabilityD2100,_PolarityC1,_PolarityC3,_PolarityC2,PAAC29,PAAC24,PAAC28,PAAC25,PAAC23,PAAC8,PAAC9,PAAC2,PAAC3,PAAC1,PAAC6,PAAC7,PAAC4,PAAC5,PAAC21,PAAC20,PAAC30,PAAC22,PAAC18,PAAC19,PAAC27,PAAC26,PAAC14,PAAC15,PAAC16,PAAC17,PAAC10,PAAC11,PAAC12,PAAC13,'_NormalizedVDWVD1075,MoranAuto_Mutability28,164,MoreauBrotoAuto_AvFlexibility19,MoreauBrotoAuto_AvFlexibility18,165,MoreauBrotoAuto_AvFlexibility15,MoreauBrotoAuto_AvFlexibility14,MoreauBrotoAuto_AvFlexibility17,MoreauBrotoAuto_AvFlexibility16,MoreauBrotoAuto_AvFlexibility11,MoreauBrotoAuto_AvFlexibility10,MoreauBrotoAuto_AvFlexibility13,MoreauBrotoAuto_AvFlexibility12,167,160,GearyAuto_AvFlexibility11,QSOSW35,MoreauBrotoAuto_Hydrophobicity30,344,345,346,347,340,341,342,343,APAAC20,163,GW,GV,GT,GS,GR,GQ,GP,MoreauBrotoAuto_Hydrophobicity1,MoreauBrotoAuto_Hydrophobicity3,MoreauBrotoAuto_Hydrophobicity2,MoreauBrotoAuto_Hydrophobicity5,MoreauBrotoAuto_Hydrophobicity4,GY,MoreauBrotoAuto_Hydrophobicity6,GG,GF,GE,GD,GC,GA,GN,GM,GL,GK,GI,GH,MoranAuto_Mutability30,011,_SolventAccessibilityD2100,012,MoranAuto_AvFlexibility26,MoranAuto_AvFlexibility27,MoranAuto_AvFlexibility24,MoranAuto_AvFlexibility25,MoranAuto_AvFlexibility22,MoranAuto_AvFlexibility23,MoranAuto_AvFlexibility20,MoranAuto_AvFlexibility21,GearyAuto_AvFlexibility1,GearyAuto_AvFlexibility2,GearyAuto_AvFlexibility3,GearyAuto_AvFlexibility4,GearyAuto_AvFlexibility5,MoranAuto_AvFlexibility28,GearyAuto_AvFlexibility7,016,017,270,271,272,273,274,_NormalizedVDWVD3025,_SecondaryStrD1050,277,tausw5,MoranAuto_Polarizability11,tausw42,374,642,K,taugrant42,taugrant7,MoreauBrotoAuto_FreeEnergy30,taugrant43,102,103,100,101,106,taugrant40,104,105,245,taugrant41,244,QSOSW34,247,246,GearyAuto_Mutability13,taugrant45,_ChargeD1025,620,612,W,077,242,074,054,055,056,057,050,051,052,053,_PolarizabilityD2001,tausw8,073,GearyAuto_AvFlexibility28,GearyAuto_Steric22,GearyAuto_Steric23,GearyAuto_Steric20,GearyAuto_Steric21,GearyAuto_Steric26,GearyAuto_Steric27,GearyAuto_Steric24,GearyAuto_Steric25,GearyAuto_Steric28,GearyAuto_Mutability18,tausw4,610,GearyAuto_AvFlexibility24,211,tausw2,tausw3,tausw36,tausw37,tausw34,tausw35,tausw32,tausw33,tausw30,V,tausw1,tausw38,tausw39,SY,615,SS,SR,SQ,SP,SW,SV,ST,SK,SI,SH,SN,SM,SL,SC,SA,SG,SF,SE,SD,MoranAuto_Hydrophobicity26,MoreauBrotoAuto_Hydrophobicity9,MoranAuto_Hydrophobicity24,MoranAuto_Hydrophobicity25,MoranAuto_Hydrophobicity22,MoranAuto_Hydrophobicity23,MoranAuto_Hydrophobicity20,MoranAuto_Hydrophobicity21,555,554,557,556,551,550,MoranAuto_Hydrophobicity28,MoranAuto_Hydrophobicity29,756,LF,LG,LD,LE,GearyAuto_Steric7,LC,GearyAuto_Steric5,GearyAuto_Steric4,LN,QSOgrant23,LL,LM,LK,LH,LI,LV,LW,LT,LR,LS,LP,LQ,614,QSOgrant20,LY,A,_ChargeD3001,MoreauBrotoAuto_ResidueASA30,QSOSW48,QSOSW49,QSOSW46,QSOSW47,QSOSW44,QSOSW45,QSOSW42,QSOSW43,AS,QSOSW41,511,MoreauBrotoAuto_Hydrophobicity7,APAAC8,APAAC9,APAAC6,APAAC7,APAAC4,APAAC5,APAAC2,APAAC3,644,APAAC1,GearyAuto_ResidueVol11,GearyAuto_ResidueVol10,GearyAuto_ResidueVol13,GearyAuto_ResidueVol12,GearyAuto_ResidueVol15,GearyAuto_ResidueVol14,GearyAuto_ResidueVol17,GearyAuto_ResidueVol16,GearyAuto_ResidueVol19,GearyAuto_ResidueVol18,514,MoreauBrotoAuto_ResidueVol2,MoreauBrotoAuto_ResidueVol3,MoreauBrotoAuto_ResidueVol1,MoreauBrotoAuto_ResidueVol6,MoreauBrotoAuto_ResidueVol7,MoreauBrotoAuto_ResidueVol4,MoreauBrotoAuto_ResidueVol5,MoreauBrotoAuto_ResidueVol8,MoreauBrotoAuto_ResidueVol9,_PolarityD3025,_SecondaryStrT23,_HydrophobicityD1001,310,GearyAuto_FreeEnergy24,GearyAuto_FreeEnergy25,GearyAuto_FreeEnergy26,GearyAuto_FreeEnergy27,GearyAuto_FreeEnergy20,GearyAuto_FreeEnergy21,GearyAuto_FreeEnergy22,GearyAuto_FreeEnergy23,GearyAuto_FreeEnergy28,GearyAuto_FreeEnergy29,MoranAuto_ResidueASA30,_SolventAccessibilityD3100,MoranAuto_Steric30,MoranAuto_Hydrophobicity14,625,_SecondaryStrD3050,624,407,406,405,404,403,402,401,400,452,GearyAuto_ResidueASA7,453,GearyAuto_ResidueASA6,454,L,GearyAuto_Hydrophobicity30,455,MoreauBrotoAuto_FreeEnergy27,GearyAuto_ResidueASA4,456,GearyAuto_ResidueASA3,GearyAuto_Polarizability20,_HydrophobicityD2050,GearyAuto_ResidueASA1,MoreauBrotoAuto_FreeEnergy9,265,_PolarizabilityD1025,MoranAuto_ResidueVol16,027,647,MoranAuto_ResidueVol17,371,370,373,372,375,MoranAuto_ResidueVol14,377,376,MoreauBrotoAuto_AvFlexibility20,MoranAuto_AvFlexibility19,MoranAuto_ResidueVol15,MoranAuto_AvFlexibility13,MoranAuto_ResidueVol12,MoranAuto_AvFlexibility11,MoranAuto_AvFlexibility10,MoranAuto_AvFlexibility17,MoranAuto_AvFlexibility16,MoranAuto_AvFlexibility15,MoranAuto_ResidueVol13,_SolventAccessibilityD2001,677,GearyAuto_ResidueASA9,GearyAuto_ResidueASA8,301,MoranAuto_Mutability25,MoranAuto_Mutability24,_SolventAccessibilityD2075,MoranAuto_Mutability26,MoranAuto_Mutability21,MoranAuto_Mutability20,MoranAuto_Mutability23,MoranAuto_Mutability22,MoranAuto_Mutability29,_SolventAccessibilityD3075,MoreauBrotoAuto_AvFlexibility28,MoreauBrotoAuto_ResidueVol13,MoreauBrotoAuto_AvFlexibility29,EM,EL,EN,EI,EH,EK,EE,ED,EG,EF,EA,EC,GearyAuto_Mutability17,GearyAuto_Mutability16,GearyAuto_Mutability15,GearyAuto_Mutability14,EY,GearyAuto_Mutability12,GearyAuto_Mutability11,GearyAuto_Mutability10,ET,EW,EV,EQ,EP,ES,ER,003,MoranAuto_Polarizability29,002,001,445,000,GearyAuto_AvFlexibility6,007,MoranAuto_AvFlexibility29,412,006,005,004,Y,177,176,175,174,173,172,171,170,MoreauBrotoAuto_Steric4,MoreauBrotoAuto_Steric5,MoreauBrotoAuto_Steric6,MoreauBrotoAuto_Steric7,MoreauBrotoAuto_Steric1,MoreauBrotoAuto_Steric2,MoreauBrotoAuto_Steric3,MoreauBrotoAuto_Steric8,MoreauBrotoAuto_Steric9,QSOSW40,061,060,063,062,065,064,067,066,QSOgrant38,QSOgrant39,_SolventAccessibilityD1050,QSOgrant31,QSOgrant32,QSOgrant33,QSOgrant34,QSOgrant35,QSOgrant36,QSOgrant37,_ChargeD2025,GearyAuto_Steric30,PAAC40,_ChargeD1050,_PolarityD1075,276,tausw43,_HydrophobicityD3075,tausw41,tausw40,tausw45,tausw44,_HydrophobicityD2075,_HydrophobicityD3001,357,356,MoranAuto_ResidueVol18,355,MoranAuto_Polarizability9,MoranAuto_Polarizability8,354,_ChargeD2050,GearyAuto_AvFlexibility13,MoranAuto_Polarizability1,MoranAuto_FreeEnergy27,MoranAuto_Polarizability3,MoranAuto_Polarizability2,MoranAuto_Polarizability5,_ChargeT12,MoranAuto_Polarizability7,MoranAuto_Polarizability6,775,774,216,_ChargeT13,tausw9,MoranAuto_Hydrophobicity30,773,772,GearyAuto_AvFlexibility16,350,M,GearyAuto_AvFlexibility17,QSOgrant4,477,GearyAuto_AvFlexibility14,QSOgrant5,GearyAuto_AvFlexibility15,MoranAuto_FreeEnergy29,QQ,QP,QS,MoranAuto_FreeEnergy28,QT,QW,QV,QY,266,QSOgrant1,QA,QC,QE,QD,QG,QF,QI,QH,QK,QSOgrant3,QM,QL,QN,MoreauBrotoAuto_ResidueASA27,MoreauBrotoAuto_ResidueASA26,MoreauBrotoAuto_ResidueASA25,MoreauBrotoAuto_ResidueASA24,MoreauBrotoAuto_ResidueASA23,MoreauBrotoAuto_ResidueASA22,MoreauBrotoAuto_ResidueASA21,MoreauBrotoAuto_ResidueASA20,QSOSW33,QSOSW32,QSOSW31,QSOSW30,QSOSW37,_PolarityD1100,MoreauBrotoAuto_ResidueASA29,MoreauBrotoAuto_ResidueASA28,QSOgrant30,MoreauBrotoAuto_ResidueVol30,_HydrophobicityD2100,212,524,RN,GearyAuto_ResidueVol28,GearyAuto_ResidueVol29,525,GearyAuto_ResidueVol24,GearyAuto_ResidueVol25,GearyAuto_ResidueVol26,GearyAuto_ResidueVol27,GearyAuto_ResidueVol20,GearyAuto_ResidueVol21,GearyAuto_ResidueVol22,GearyAuto_ResidueVol23,GearyAuto_ResidueVol9,GearyAuto_ResidueVol8,540,527,546,547,544,545,GearyAuto_ResidueVol1,520,GearyAuto_ResidueVol3,GearyAuto_ResidueVol2,GearyAuto_ResidueVol5,GearyAuto_ResidueVol4,GearyAuto_ResidueVol7,GearyAuto_ResidueVol6,761,GearyAuto_ResidueASA5,523,_SecondaryStrD3001,767,GearyAuto_FreeEnergy30,MoranAuto_ResidueASA29,MoranAuto_ResidueASA28,MoranAuto_ResidueASA25,MoranAuto_ResidueASA24,MoranAuto_ResidueASA27,MoranAuto_ResidueASA26,MoranAuto_ResidueASA21,MoranAuto_ResidueASA20,MoranAuto_ResidueASA23,MoranAuto_ResidueASA22,MoranAuto_Steric21,MoranAuto_Steric20,MoranAuto_Steric23,MoranAuto_Steric22,MoranAuto_Steric25,MoranAuto_Steric24,MoranAuto_Steric27,MoranAuto_Steric26,MoranAuto_Steric29,MoranAuto_Steric28,GearyAuto_Polarizability16,GearyAuto_Polarizability17,GearyAuto_Polarizability14,GearyAuto_Polarizability15,GearyAuto_Polarizability12,GearyAuto_Polarizability13,GearyAuto_Polarizability10,GearyAuto_Polarizability11,GearyAuto_Polarizability18,GearyAuto_Polarizability19,_HydrophobicityD1100,GearyAuto_Hydrophobicity29,GearyAuto_Hydrophobicity28,GearyAuto_Hydrophobicity25,GearyAuto_Hydrophobicity24,GearyAuto_Hydrophobicity27,GearyAuto_Hydrophobicity26,GearyAuto_Hydrophobicity21,GearyAuto_Hydrophobicity20,GearyAuto_Hydrophobicity23,GearyAuto_Hydrophobicity22,_PolarityD2100,166,GearyAuto_ResidueASA2,GearyAuto_Mutability9,GearyAuto_Mutability8,GearyAuto_Mutability7,GearyAuto_Mutability6,GearyAuto_Mutability5,GearyAuto_Mutability4,GearyAuto_Mutability3,GearyAuto_Mutability2,GearyAuto_Mutability1,443,366,RK,364,365,362,363,360,361,C,_HydrophobicityT13,_HydrophobicityT12,MoreauBrotoAuto_Polarizability11,MoreauBrotoAuto_Hydrophobicity15,MoreauBrotoAuto_Hydrophobicity14,MoreauBrotoAuto_Hydrophobicity17,MoreauBrotoAuto_Hydrophobicity16,MoreauBrotoAuto_Hydrophobicity11,MoreauBrotoAuto_Hydrophobicity10,MoreauBrotoAuto_Hydrophobicity13,MoreauBrotoAuto_Hydrophobicity12,_PolarityT23,367,MoreauBrotoAuto_Hydrophobicity19,MoreauBrotoAuto_Hydrophobicity18,MoranAuto_Mutability10,MoranAuto_Mutability11,MoranAuto_Mutability12,MoranAuto_Mutability13,MoranAuto_Mutability14,MoranAuto_Mutability15,MoranAuto_Mutability16,MoranAuto_Mutability17,MoranAuto_Mutability18,MoranAuto_Mutability19,663,114,QSOgrant2,MoreauBrotoAuto_FreeEnergy29,MoreauBrotoAuto_FreeEnergy28,MoreauBrotoAuto_Polarizability18,_SecondaryStrD2075,MoranAuto_Polarizability10,MoreauBrotoAuto_FreeEnergy23,MoranAuto_Polarizability12,MoreauBrotoAuto_FreeEnergy25,MoreauBrotoAuto_FreeEnergy24,MoranAuto_Polarizability17,MoreauBrotoAuto_FreeEnergy26,DM,252,253,250,251,256,_SecondaryStrD1075,254,_NormalizedVDWVD3001,107,670,CK,CI,CH,CN,CM,MoreauBrotoAuto_Polarizability8,CC,MoreauBrotoAuto_Polarizability6,MoreauBrotoAuto_Polarizability5,MoreauBrotoAuto_Polarizability4,CG,MoreauBrotoAuto_Polarizability2,CE,N,CY,353,CS,CR,_PolarizabilityD3001,CP,CW,CV,162,CT,PAAC14,GearyAuto_ResidueASA26,GearyAuto_ResidueASA27,GearyAuto_ResidueASA24,GearyAuto_ResidueASA25,GearyAuto_ResidueASA22,GearyAuto_ResidueASA23,GearyAuto_ResidueASA20,GearyAuto_ResidueASA21,GearyAuto_ResidueASA28,GearyAuto_ResidueASA29,076,QSOSW28,HF,075,072,PAAC15,070,QSOSW29,662,MoreauBrotoAuto_Steric12,QSOgrant28,QSOgrant27,QSOgrant26,QSOgrant25,QSOgrant24,_NormalizedVDWVT23,QSOgrant22,QSOgrant21,MoreauBrotoAuto_ResidueASA17,QSOSW15,MoreauBrotoAuto_AvFlexibility8,_PolarizabilityD2025,112,QSOSW17,PAAC32,PAAC33,PAAC30,PAAC31,PAAC36,PAAC37,PAAC34,PAAC35,MoreauBrotoAuto_Steric30,QSOSW11,PAAC38,PAAC39,QSOSW22,HD,QSOSW10,VA,VC,VD,VE,VF,VG,VH,VI,VK,VL,VM,VN,671,VP,VQ,VR,VS,VT,_PolarizabilityT23,VV,VW,_HydrophobicityD3100,VY,673,MoreauBrotoAuto_Mutability30,434,435,MoreauBrotoAuto_AvFlexibility5,QSOgrant50,641,T,QSOSW18,422,430,431,PAAC11,_NormalizedVDWVD3050,MoranAuto_ResidueVol8,MoranAuto_ResidueVol9,MoranAuto_ResidueVol4,MoranAuto_ResidueVol5,MoranAuto_ResidueVol6,MoranAuto_ResidueVol7,766,MoranAuto_ResidueVol1,MoranAuto_ResidueVol2,MoranAuto_ResidueVol3,MoranAuto_Steric2,MoranAuto_Steric3,MoranAuto_Steric1,MoranAuto_Steric6,MoranAuto_Steric7,MoranAuto_Steric4,MoranAuto_Steric5,MoranAuto_Steric8,MoranAuto_Steric9,MoreauBrotoAuto_Polarizability13,MoreauBrotoAuto_Polarizability12,PAAC8,PAAC9,MoreauBrotoAuto_Polarizability17,MoreauBrotoAuto_Polarizability16,MoreauBrotoAuto_Polarizability15,MoreauBrotoAuto_Polarizability14,PAAC2,PAAC3,MoreauBrotoAuto_Polarizability19,PAAC1,PAAC6,PAAC7,PAAC4,PAAC5,334,PAAC13,337,MoreauBrotoAuto_ResidueASA12,MoreauBrotoAuto_ResidueASA13,MoreauBrotoAuto_ResidueASA10,MoreauBrotoAuto_ResidueASA11,MoreauBrotoAuto_ResidueASA16,D,MoreauBrotoAuto_ResidueASA14,MoreauBrotoAuto_ResidueASA15,QSOSW20,QSOSW21,MoreauBrotoAuto_ResidueASA18,MoreauBrotoAuto_ResidueASA19,QSOSW24,QSOSW25,QSOSW26,QSOSW27,752,WN,645,GearyAuto_ResidueVol30,MoranAuto_Polarizability4,537,536,535,534,533,532,531,530,_SolventAccessibilityT13,_SolventAccessibilityT12,HY,_ChargeT23,HR,HS,HP,HQ,HV,HW,HT,MoranAuto_ResidueASA18,HK,HH,HI,HN,HL,HM,MoranAuto_ResidueASA10,MoranAuto_ResidueASA11,MoranAuto_ResidueASA12,HA,_SecondaryStrD3025,MoranAuto_ResidueASA15,MoranAuto_ResidueASA16,MoranAuto_ResidueASA17,MoranAuto_Steric14,MoranAuto_Steric15,MoranAuto_Steric16,MoranAuto_Steric17,MoreauBrotoAuto_ResidueVol29,MoreauBrotoAuto_ResidueVol28,MoranAuto_Steric12,MoranAuto_Steric13,MoreauBrotoAuto_ResidueVol25,GearyAuto_AvFlexibility8,MoreauBrotoAuto_ResidueVol27,MoreauBrotoAuto_ResidueVol26,MoranAuto_Steric18,MoranAuto_Steric19,MoreauBrotoAuto_ResidueVol23,MoreauBrotoAuto_ResidueVol22,421,420,423,_ChargeD3025,425,424,427,646,GearyAuto_Hydrophobicity10,GearyAuto_Hydrophobicity11,GearyAuto_Hydrophobicity12,GearyAuto_Hydrophobicity13,_HydrophobicityD1050,GearyAuto_Hydrophobicity15,GearyAuto_Hydrophobicity16,GearyAuto_Hydrophobicity17,GearyAuto_Hydrophobicity18,GearyAuto_Hydrophobicity19,302,_NormalizedVDWVD1025,304,305,306,307,_PolarityD3050,GearyAuto_Mutability19,_PolarityD2050,_SolventAccessibilityC1,_SolventAccessibilityC2,_SolventAccessibilityC3,_PolarizabilityD1001,147,QSOgrant17,_HydrophobicityT23,_SolventAccessibilityD2025,_NormalizedVDWVC3,_PolarityT12,_PolarityT13,643,_SolventAccessibilityD3001,_PolarizabilityD1075,MoranAuto_FreeEnergy18,MoranAuto_FreeEnergy19,MoranAuto_FreeEnergy12,MoranAuto_FreeEnergy13,MoranAuto_FreeEnergy10,MoranAuto_FreeEnergy11,MoranAuto_FreeEnergy16,MoranAuto_FreeEnergy17,MoranAuto_FreeEnergy14,MoranAuto_FreeEnergy15,GearyAuto_AvFlexibility9,_PolarizabilityC2,_PolarizabilityC3,_PolarizabilityC1,MoranAuto_Polarizability24,MoranAuto_Polarizability25,MoranAuto_Polarizability26,MoranAuto_Polarizability27,MoranAuto_Polarizability20,MoranAuto_Polarizability21,MoranAuto_Polarizability22,_ChargeD1001,QR,_PolarityD2025,227,226,225,224,223,222,221,220,GearyAuto_Steric29,MoranAuto_ResidueASA19,771,117,_HydrophobicityD1025,151,150,153,152,155,113,157,156,045,QSOSW14,TW,HC,770,MoranAuto_ResidueASA13,141,TY,TV,MoranAuto_ResidueASA14,TT,TR,TS,TP,HG,TN,TL,TM,GearyAuto_ResidueASA30,TH,TI,TF,TG,TD,HE,TC,TA,AA,AC,E,AD,AG,AF,AI,AH,AK,AM,AL,AN,GearyAuto_Hydrophobicity6,GearyAuto_Hydrophobicity7,_PolarizabilityD3100,GearyAuto_Hydrophobicity5,GearyAuto_Hydrophobicity2,GearyAuto_Hydrophobicity3,AW,GearyAuto_Hydrophobicity1,AY,MoranAuto_Steric10,GearyAuto_Hydrophobicity8,GearyAuto_Hydrophobicity9,PAAC21,PAAC20,GearyAuto_Steric19,GearyAuto_Steric18,PAAC25,PAAC24,PAAC27,PAAC26,GearyAuto_Steric13,PAAC28,GearyAuto_Steric11,GearyAuto_Steric10,_PolarityC1,GearyAuto_Steric16,_PolarityC3,_PolarityC2,MoreauBrotoAuto_ResidueVol21,tausw31,_SolventAccessibilityD1025,QSOSW39,MoreauBrotoAuto_ResidueVol20,010,tausw26,_HydrophobicityD3050,QSOSW38,014,015,MoreauBrotoAuto_Mutability29,MoreauBrotoAuto_Mutability28,MoreauBrotoAuto_Mutability27,MoreauBrotoAuto_Mutability26,MoreauBrotoAuto_Mutability25,MoreauBrotoAuto_Mutability24,MoreauBrotoAuto_Mutability23,MoreauBrotoAuto_Mutability22,MoreauBrotoAuto_Mutability21,MoreauBrotoAuto_Mutability20,QSOgrant49,QSOgrant48,GearyAuto_Mutability28,QSOgrant41,QSOgrant40,QSOgrant43,QSOgrant42,QSOgrant45,QSOgrant44,QSOgrant47,QSOgrant46,GearyAuto_Mutability26,_ChargeD2075,WM,QSOSW12,235,GearyAuto_Mutability25,717,WK,715,714,713,tausw20,711,QSOSW36,P,GearyAuto_Mutability20,MoranAuto_Hydrophobicity19,MoranAuto_Hydrophobicity18,MoranAuto_Hydrophobicity17,MoranAuto_Hydrophobicity16,MoranAuto_Hydrophobicity15,GearyAuto_Mutability21,MoranAuto_Hydrophobicity13,MoranAuto_Hydrophobicity12,MoranAuto_Hydrophobicity11,MoranAuto_Hydrophobicity10,MoreauBrotoAuto_Polarizability26,MoreauBrotoAuto_Polarizability27,MoreauBrotoAuto_Polarizability24,MoreauBrotoAuto_Polarizability25,MoreauBrotoAuto_Polarizability22,MoreauBrotoAuto_Polarizability23,MoreauBrotoAuto_Polarizability20,MoreauBrotoAuto_Polarizability21,MoreauBrotoAuto_FreeEnergy6,MoreauBrotoAuto_FreeEnergy7,MoreauBrotoAuto_FreeEnergy4,MoreauBrotoAuto_FreeEnergy5,MoreauBrotoAuto_FreeEnergy2,MoreauBrotoAuto_FreeEnergy3,MoreauBrotoAuto_Polarizability28,MoreauBrotoAuto_Polarizability29,426,AR,_PolarityD1050,WQ,tausw22,GearyAuto_Hydrophobicity14,313,MoreauBrotoAuto_Steric28,tausw25,_ChargeD2100,taugrant15,taugrant14,taugrant17,taugrant16,taugrant11,taugrant10,taugrant13,taugrant12,312,taugrant19,taugrant18,303,tausw24,GearyAuto_Polarizability30,442,441,440,447,446,712,444,QSOSW19,ME,MD,MG,MF,MA,MC,MM,ML,MN,MI,MH,MK,MoreauBrotoAuto_AvFlexibility9,MT,MW,MV,MQ,MP,MS,MR,MoreauBrotoAuto_AvFlexibility1,437,MoreauBrotoAuto_AvFlexibility3,MoreauBrotoAuto_AvFlexibility2,MY,MoreauBrotoAuto_AvFlexibility4,MoreauBrotoAuto_AvFlexibility7,MoreauBrotoAuto_AvFlexibility6,FP,FQ,FR,FS,FT,FV,FW,335,FY,710,336,331,330,333,332,FA,317,FC,FD,FE,FF,FG,FH,FI,FK,FL,FM,FN,745,542,746,543,747,316,541,041,QSOgrant16,F,326,553,552,YY,_PolarizabilityD1100,QSOgrant11,YI,YH,YK,MoranAuto_Polarizability30,YM,YL,YN,YA,YC,YE,YD,YG,YF,_NormalizedVDWVD2001,234,_SecondaryStrD2050,YQ,YP,YS,YR,230,YT,YW,YV,237,716,231,MoranAuto_Hydrophobicity27,232,233,146,GearyAuto_Steric3,144,NV,taugrant9,taugrant8,140,GearyAuto_Steric2,_SecondaryStrD1001,taugrant4,_NormalizedVDWVD3075,taugrant6,taugrant1,GearyAuto_Steric1,taugrant3,taugrant2,Q,_ChargeC3,GearyAuto_Steric6,_ChargeC1,121,LA,GearyAuto_Steric9,125,GearyAuto_Steric8,GearyAuto_AvFlexibility30,133,MoranAuto_FreeEnergy1,MoranAuto_FreeEnergy2,MoranAuto_FreeEnergy3,MoranAuto_FreeEnergy4,MoranAuto_FreeEnergy5,MoranAuto_FreeEnergy6,MoranAuto_FreeEnergy7,MoranAuto_FreeEnergy8,MoranAuto_FreeEnergy9,_PolarizabilityD2050,RT,RV,RW,RP,RQ,RR,RS,RY,137,130,RD,RE,RF,RG,PAAC18,RA,MoreauBrotoAuto_Steric18,MoreauBrotoAuto_Steric19,MoreauBrotoAuto_Steric16,MoreauBrotoAuto_Steric17,PAAC16,MoreauBrotoAuto_Steric15,RH,MoreauBrotoAuto_Steric13,PAAC12,MoreauBrotoAuto_Steric11,_NormalizedVDWVD3100,MoranAuto_Steric11,025,024,MoreauBrotoAuto_Mutability18,MoreauBrotoAuto_Mutability19,021,020,023,022,MoreauBrotoAuto_Mutability12,MoreauBrotoAuto_Mutability13,MoreauBrotoAuto_Mutability10,MoreauBrotoAuto_Mutability11,MoreauBrotoAuto_Mutability16,MoreauBrotoAuto_Mutability17,MoreauBrotoAuto_Mutability14,MoreauBrotoAuto_Mutability15,GearyAuto_Polarizability8,GearyAuto_Polarizability9,GearyAuto_Polarizability4,GearyAuto_Polarizability5,GearyAuto_Polarizability6,GearyAuto_Polarizability7,GearyAuto_Polarizability1,GearyAuto_Polarizability2,GearyAuto_Polarizability3,510,_NormalizedVDWVD2075,704,705,706,707,700,701,702,703,145,275,142,143,414,_HydrophobicityD2025,762,415,MoreauBrotoAuto_Polarizability30,_SecondaryStrD2025,416,taugrant5,417,613,410,MoreauBrotoAuto_FreeEnergy8,QSOSW16,_PolarityD1025,tausw29,611,763,411,G,tausw28,616,413,617,MoranAuto_Mutability6,MoranAuto_Mutability7,MoranAuto_Mutability4,MoranAuto_Mutability5,MoranAuto_Mutability2,MoranAuto_Mutability3,MoranAuto_Mutability1,MoranAuto_ResidueVol19,MoranAuto_Mutability8,MoranAuto_Mutability9,taugrant28,taugrant29,760,515,taugrant20,taugrant21,taugrant22,taugrant23,taugrant24,taugrant25,taugrant26,taugrant27,MoranAuto_Hydrophobicity9,MoranAuto_Hydrophobicity8,MoranAuto_Hydrophobicity7,MoranAuto_Hydrophobicity6,MoranAuto_Hydrophobicity5,MoranAuto_Hydrophobicity4,MoranAuto_Hydrophobicity3,MoranAuto_Hydrophobicity2,MoranAuto_Hydrophobicity1,516,623,622,621,_PolarityD3100,627,626,GearyAuto_Polarizability29,GearyAuto_Polarizability28,GearyAuto_Polarizability27,GearyAuto_Polarizability26,GearyAuto_Polarizability25,GearyAuto_Polarizability24,GearyAuto_Polarizability23,GearyAuto_Polarizability22,GearyAuto_Polarizability21,_ChargeD3050,tausw27,236,NG,517,MoreauBrotoAuto_Hydrophobicity8,NY,656,657,654,655,652,653,650,651,311,MoranAuto_ResidueVol22,MoreauBrotoAuto_Mutability8,MoreauBrotoAuto_Mutability9,GearyAuto_AvFlexibility18,GearyAuto_AvFlexibility19,_HydrophobicityD1075,GearyAuto_AvFlexibility12,MoreauBrotoAuto_Mutability1,GearyAuto_AvFlexibility10,MoreauBrotoAuto_Mutability3,MoreauBrotoAuto_Mutability4,MoreauBrotoAuto_Mutability5,MoreauBrotoAuto_Mutability6,MoreauBrotoAuto_Mutability7,R,315,764,314,_PolarityD2075,KC,KA,KG,KF,KE,KD,KK,_ChargeC2,KI,KH,KN,KM,KL,KS,KR,KQ,KP,KW,KV,KT,765,KY,NW,DN,753,DL,MoreauBrotoAuto_AvFlexibility30,DK,DH,DI,DF,DG,DD,DE,MoreauBrotoAuto_Mutability2,DC,DA,026,_SolventAccessibilityD3025,DY,DV,DW,DT,DR,DS,DP,DQ,MoreauBrotoAuto_Steric29,_PolarizabilityD1050,036,352,MoranAuto_ResidueVol10,037,WL,323,320,321,MoranAuto_FreeEnergy30,327,324,325,035,562,_NormalizedVDWVD1050,201,200,203,202,205,204,207,206,217,351,MoreauBrotoAuto_FreeEnergy1,_SecondaryStrD2100,513,MoranAuto_Polarizability23,MoreauBrotoAuto_ResidueASA1,MoreauBrotoAuto_ResidueASA2,MoreauBrotoAuto_ResidueASA3,MoreauBrotoAuto_ResidueASA4,MoreauBrotoAuto_ResidueASA5,MoreauBrotoAuto_ResidueASA6,MoreauBrotoAuto_ResidueASA7,MoreauBrotoAuto_ResidueASA8,MoreauBrotoAuto_ResidueASA9,572,WG,WF,WE,WD,WC,WA,GearyAuto_Mutability29,_NormalizedVDWVC1,GearyAuto_Mutability27,GearyAuto_Mutability24,_NormalizedVDWVC2,GearyAuto_Mutability22,GearyAuto_Mutability23,WI,WH,WW,WV,WT,WS,WR,MoranAuto_ResidueVol11,WP,_ChargeD1075,WY,664,526,674,MoranAuto_FreeEnergy26,450,GearyAuto_ResidueASA19,GearyAuto_ResidueASA18,GearyAuto_ResidueASA17,GearyAuto_ResidueASA16,GearyAuto_ResidueASA15,GearyAuto_ResidueASA14,GearyAuto_ResidueASA13,GearyAuto_ResidueASA12,GearyAuto_ResidueASA11,GearyAuto_ResidueASA10,_SecondaryStrD1100,675,H,120,MoreauBrotoAuto_Polarizability10,122,123,124,MoreauBrotoAuto_Steric23,126,127,214,777,_PolarityD2001,676,132,_HydrophobicityC1,_HydrophobicityC2,_HydrophobicityC3,131,512,MoreauBrotoAuto_Steric22,_PolarizabilityD3050,776,136,QSOSW13,521,_SolventAccessibilityD1001,135,134,032,033,030,031,_PolarizabilityT13,_PolarizabilityT12,034,522,PR,PS,PP,PQ,PV,PW,PT,PY,PC,PA,PF,PG,PD,PE,PK,PH,PI,PN,PL,PM,_ChargeD1100,QSOgrant29,MoreauBrotoAuto_Steric27,S,_PolarizabilityD3025,672,IQ,tausw19,MoranAuto_ResidueVol30,tausw15,tausw16,tausw17,tausw10,tausw11,tausw12,tausw13,605,604,573,MoreauBrotoAuto_Steric25,571,570,_PolarityD3075,576,575,574,QSOSW1,QSOSW3,QSOSW2,QSOSW5,QSOSW4,QSOSW7,QSOSW6,MoreauBrotoAuto_Polarizability3,QSOSW8,154,726,603,taugrant39,MoreauBrotoAuto_Steric24,MoranAuto_AvFlexibility12,taugrant36,taugrant35,taugrant34,taugrant33,taugrant32,taugrant31,taugrant30,731,730,733,732,735,734,CF,736,506,507,504,505,502,503,500,PAAC19,630,631,632,633,634,635,636,637,_SecondaryStrC2,_SecondaryStrC3,467,RC,461,460,463,462,RL,457,577,RM,tausw18,MoreauBrotoAuto_Steric14,465,PAAC17,PAAC10,RI,tausw14,MoreauBrotoAuto_Steric10,MoranAuto_ResidueASA8,MoranAuto_ResidueASA9,MoranAuto_ResidueASA6,MoranAuto_ResidueASA7,MoranAuto_ResidueASA4,MoranAuto_ResidueASA5,MoranAuto_ResidueASA2,MoranAuto_ResidueASA3,taugrant44,MoranAuto_ResidueASA1,322,GearyAuto_AvFlexibility29,_NormalizedVDWVD1001,GearyAuto_AvFlexibility27,GearyAuto_AvFlexibility26,GearyAuto_AvFlexibility25,_PolarityD3001,GearyAuto_AvFlexibility23,GearyAuto_AvFlexibility22,GearyAuto_AvFlexibility21,GearyAuto_AvFlexibility20,_SecondaryStrD3100,_SecondaryStrD3075,MoranAuto_AvFlexibility9,MoranAuto_AvFlexibility8,MoranAuto_AvFlexibility7,MoranAuto_AvFlexibility6,MoranAuto_AvFlexibility5,MoranAuto_AvFlexibility4,MoranAuto_AvFlexibility3,MoranAuto_AvFlexibility2,MoranAuto_AvFlexibility1,MoranAuto_Polarizability28,466,IY,MoreauBrotoAuto_AvFlexibility21,MoreauBrotoAuto_AvFlexibility22,MoreauBrotoAuto_AvFlexibility23,MoreauBrotoAuto_AvFlexibility24,MoreauBrotoAuto_AvFlexibility25,MoreauBrotoAuto_AvFlexibility26,MoreauBrotoAuto_AvFlexibility27,I,IP,IS,IR,IT,IW,IV,II,IH,IK,IM,IL,_SolventAccessibilityD3050,IN,IA,IC,IE,ID,IG,IF,MoreauBrotoAuto_Hydrophobicity28,MoreauBrotoAuto_Hydrophobicity29,MoreauBrotoAuto_Hydrophobicity24,MoreauBrotoAuto_Hydrophobicity25,MoreauBrotoAuto_Hydrophobicity26,MoreauBrotoAuto_Hydrophobicity27,MoreauBrotoAuto_Hydrophobicity20,MoreauBrotoAuto_Hydrophobicity21,MoreauBrotoAuto_Hydrophobicity22,MoreauBrotoAuto_Hydrophobicity23,MoranAuto_FreeEnergy23,MoranAuto_FreeEnergy22,MoranAuto_FreeEnergy21,MoranAuto_FreeEnergy20,APAAC18,APAAC19,MoranAuto_FreeEnergy25,MoranAuto_FreeEnergy24,APAAC14,APAAC15,APAAC16,APAAC17,APAAC10,APAAC11,APAAC12,APAAC13,013,MoranAuto_Polarizability19,QSOSW23,215,071,213,210,MoranAuto_Polarizability18,241,TQ,_SolventAccessibilityD2050,MoreauBrotoAuto_FreeEnergy21,432,MoranAuto_AvFlexibility30,MoreauBrotoAuto_FreeEnergy20,MoranAuto_Polarizability13,TK,_NormalizedVDWVD2050,MoreauBrotoAuto_FreeEnergy22,_SecondaryStrD2001,MoranAuto_Polarizability15,GearyAuto_Mutability30,MoranAuto_Polarizability14,263,262,261,260,267,_NormalizedVDWVD2025,_SecondaryStrD1025,264,MoranAuto_Polarizability16,436,TE,MoreauBrotoAuto_FreeEnergy10,MoreauBrotoAuto_FreeEnergy11,MoreauBrotoAuto_FreeEnergy12,MoreauBrotoAuto_FreeEnergy13,MoreauBrotoAuto_FreeEnergy14,MoreauBrotoAuto_FreeEnergy15,MoreauBrotoAuto_FreeEnergy16,MoreauBrotoAuto_FreeEnergy17,MoreauBrotoAuto_FreeEnergy18,MoreauBrotoAuto_FreeEnergy19,GearyAuto_FreeEnergy19,GearyAuto_FreeEnergy18,QSOgrant8,115,QSOSW9,740,116,111,110,667,_PolarizabilityD3075,AE,243,_PolarizabilityD2100,QSOgrant9,257,741,666,_NormalizedVDWVD1100,255,CQ,_SolventAccessibilityD1075,433,742,GearyAuto_FreeEnergy14,665,_NormalizedVDWVT12,_NormalizedVDWVT13,047,046,tausw6,044,043,042,_HydrophobicityD3025,040,_PolarizabilityD2075,AQ,QSOgrant14,QSOgrant15,QSOgrant12,QSOgrant13,QSOgrant10,AP,_ChargeD2001,taugrant38,743,GearyAuto_Hydrophobicity4,QSOgrant18,QSOgrant19,GearyAuto_FreeEnergy1,taugrant37,GearyAuto_FreeEnergy3,GearyAuto_FreeEnergy2,GearyAuto_FreeEnergy5,GearyAuto_FreeEnergy4,GearyAuto_FreeEnergy7,GearyAuto_FreeEnergy6,GearyAuto_FreeEnergy9,GearyAuto_FreeEnergy8,AT,_SolventAccessibilityD1100,AV,_NormalizedVDWVD2100,240,MoranAuto_AvFlexibility18,161,NH,NI,NK,NL,NM,NN,tausw21,NA,tausw23,NC,ND,NE,NF,_HydrophobicityD2001,MoranAuto_ResidueVol27,MoranAuto_ResidueVol26,MoranAuto_ResidueVol25,MoranAuto_ResidueVol24,MoranAuto_ResidueVol23,_PolarityD1001,MoranAuto_ResidueVol21,MoranAuto_ResidueVol20,NP,NQ,NR,NS,NT,MoranAuto_ResidueVol29,MoranAuto_ResidueVol28,300,QSOgrant6,751,750,757,737,755,754,560,561,661,563,564,565,566,567,MoranAuto_AvFlexibility14,QSOgrant7,660,tausw7,727,724,725,722,723,720,721,501,_ChargeD3075,QSOSW50,607,606,601,600,640,602,MoreauBrotoAuto_ResidueVol14,MoreauBrotoAuto_ResidueVol15,MoreauBrotoAuto_ResidueVol16,MoreauBrotoAuto_ResidueVol17,MoreauBrotoAuto_ResidueVol10,MoreauBrotoAuto_ResidueVol11,MoreauBrotoAuto_ResidueVol12,MoreauBrotoAuto_ResidueVol24,MoreauBrotoAuto_ResidueVol18,MoreauBrotoAuto_ResidueVol19,_SolventAccessibilityT23,MoreauBrotoAuto_Polarizability9,CL,PAAC23,MoreauBrotoAuto_Polarizability7,PAAC22,464,_ChargeD3100,CA,744,_SecondaryStrT13,_SecondaryStrT12,_SecondaryStrC1,GearyAuto_FreeEnergy11,GearyAuto_FreeEnergy10,GearyAuto_FreeEnergy13,GearyAuto_FreeEnergy12,GearyAuto_FreeEnergy15,MoranAuto_Mutability27,GearyAuto_FreeEnergy17,GearyAuto_FreeEnergy16,PAAC29,MoreauBrotoAuto_Polarizability1,GearyAuto_Steric12,CD,MoreauBrotoAuto_Steric21,MoreauBrotoAuto_Steric20,GearyAuto_Steric17,451,MoreauBrotoAuto_Steric26,472,473,470,471,476,GearyAuto_Steric15,474,475,GearyAuto_Steric14\n")

#changing the current directory to protein_repository
os.chdir(path2)
#start program
#read the filenames  and save the list of files in dirs
#dirs =os.listdir(path)
#or we can use glob instead to match pdb files automatically
dirs= glob.glob('*.pdb')
print dirs

class PdbSeqError(Exception):
    """
    Error class for this module.
    """

    pass

#locating the period to remove the extension from the filename
def locPeriod(x):
    for i in range(len(x)):
        if x[i]=='.':
            return i
            break #it will return the first period
def locSlash(x):
    for i in range(len(x)):
        if x[i]=='\n':
            return i
            break #it will return the first period
def locComma(x):
    for i in range(len(x)):
        if x[i]==',':
            return i+1
            break
def locChar(x,y):
    for i in range(len(x)):
        if x[i]==y:
            return i
            break
#reading all the files from list file
fcsv = open('native_pdb.csv','r')
for line in fcsv:
    array =line.split(',')
    if(array[0].startswith('T0') and len(array)==2):
        pdb_dic[array[0]]=(array[1])[:locChar(array[1],'\t')]
    else:
        continue
print pdb_dic
list_dic=pdb_dic.keys()
list_dic.sort()

#dividing into 20


#function to get those keys in pdb_dic that have corresponding native pdb files in protein repository
def delInValidfile(name):
    pop=1
    pdbl=glob.glob("*.pdb")
    if(pdb_dic[name]+".pdb" in pdbl):
        pop=0
    if(pop==1):
        pdb_dic.pop(name)  
    return pop

#python script to get the seq from atom file
def pdbSeq(pdb,use_atoms=False):
    """
    use the ATOM 
    entries.  Return dictionary of sequences keyed to chain and type of
    sequence used.
    
    """
    seq_type = "ATOM  "
    atoms = []
    for l in pdb:
        if l[0:6] == "ATOM  " and l[13:16] == "CA ":
                
            # Check to see if this is a second conformation of the previous
            # atom
            if len(atoms) != 0:
                if atoms[-1][17:26] == l[17:26]:
                    continue

            atoms.append(l)
        elif l[0:6] == "HETATM" and l[13:16] == "CA " and l[17:20] == "MSE":

            # Check to see if this is a second conformation of the previous
            # atom
            if len(atoms) != 0:
                if atoms[-1][17:26] == l[17:26]:
                    continue

            atoms.append(l)

    chain_dict = dict([(l[21],[]) for l in atoms])
    for c in chain_dict.keys():
        chain_dict[c] = [l[17:20] for l in atoms if l[21] == c]

    return chain_dict,seq_type

def convertModifiedAA(chain_dict,pdb):
    """
    Convert modified amino acids to their normal counterparts.
    """

    # See if there are any non-standard amino acids in the pdb file.  If there
    # are not, return
    modres = [l for l in pdb if l[0:6] == "MODRES"]
    if len(modres) == 0:
        return chain_dict

    # Create list of modified residues
    mod_dict = dict([(l[12:15],l[24:27]) for l in modres])

    # Replace all entries in chain_dict with their unmodified counterparts.
    for c in chain_dict.keys():
        for i, a in enumerate(chain_dict[c]):
            if mod_dict.has_key(a):
                chain_dict[c][i] = mod_dict[a]

    return chain_dict

def pdbSeq2Fasta(pdb,pdb_id="",chain="all",use_atoms=True):
    """
    Extract sequence from pdb file and write out in FASTA format.
    """

    # Grab sequences
    chain_dict, seq_type = pdbSeq(pdb,use_atoms)

    # Convert modified amino acids to their natural counterparts
    chain_dict = convertModifiedAA(chain_dict,pdb)

    # Determine which chains are being written out
    if chain == "all":
        chains_to_write = chain_dict.keys()
        chains_to_write.sort()
    else:
        if chain in chain_dict.keys():
            chains_to_write = [chain]
        else:
            err = "Chain \"%s\" not in pdb!" % chain
            raise PdbSeqError(err)

    # Convert sequences to 1-letter format and join strings
    for c in chains_to_write:
        for aa_index, aa in enumerate(chain_dict[c]):
            try:
                chain_dict[c][aa_index] = d[aa]
            except KeyError:
                chain_dict[c][aa_index] = "X"

    out = []
    for c in chains_to_write:
        #out.append(">%s%s_%s\n" % (pdb_id,c,seq_type))

        # Write output in lines 80 characters long
        seq_length = len(chain_dict[c])
        num_lines = seq_length / 80
        
        for i in range(num_lines+1):
            out.append("".join([aa for aa in chain_dict[c][80*i:80*(i+1)]]))
            #out.append("\n")
        out.append("".join([aa for aa in chain_dict[c][80*(i+1):]]))
        #out.append("\n")


    return "".join(out)


#getting folders name according to the serial no of the program running
no_of_cores=multiprocessing.cpu_count()
size=int(len(list_dic)/no_of_cores)
s=sno*size
e=size+s

print "size   ",size
print list_dic[s:e]
#starting of creating dataset
for key in list_dic[s:e]:   
    time1=time.time()
    #will skip that key that has no native file in protein repository or have no key in protein repository
    os.chdir(path1)
    lst=glob.glob("T0*")
    lst.sort()
    if key not in lst:
        print "\nfolder of \t"+key+"  doesnt exist\n"
        continue
    if not pdb_dic[key]:   #native file of the folder is not present in list_dic
        print"native file of the folder\t"+key+"\tdoesnt exist\n"
        continue
    
    print "calculating features of the folder:\t",key
    #evaluation parameters calculation
    val=pdb_dic[key]
    #way to calculate the features of modelled structures and native and also write their scores
    try:
        source=path1+"\\"+val+".pdb"
        dest=path1+"\\"+key
        shutil.copy(source,dest) #copies the file to the destination
        shutil.copy(path1+"\\"+"ref_score.exe",dest) #copies the exe file to each folder
    except:
        print "\nnative file :\t"+val+"  doesnt exist "
        continue
    #now we have to change extension of all the files to .pdb
    os.chdir(dest)
    match=key+"*"
    to_file=glob.glob(match)
    for f in to_file:    
        os.rename(f,f[:locPeriod(f)]+".pdb")
    
    get_file=glob.glob("*.pdb")
    for f in get_file:    
        try:    
            model=f
            native=val[:4]
        
            cmd="ref_score.exe"+" "+model+" "+native+".pdb"+""+">score.txt"
            os.system("F:")
            os.system("cd "+path1)
            os.system(cmd)
            ft = open('score.txt','r')
            line = ft.readlines()
            if(len(line)<=3):
                break
            #features calculation
        
            #working on protein sequences
            
            filep=open(f,'r')
            seq= pdbSeq2Fasta(filep)
            #filep.close()
            
            #seq = getpdb.GetSeqFromPDB(path1+'\\'+val[:4]+".pdb")
            f15=len(seq)
            #calculating amino acid composition descriptors  total 20 descriptors
            res=AAComposition.CalculateAAComposition(seq)
            k_list1=','.join(c for c in res.keys())
            v_list1=','.join(str(c) for c in res.values())
            #CTD  calculates 147 descriptors
            ctd = CTD.CalculateCTD(seq)
            k_list2=','.join(c for c in ctd.keys())
            v_list2=','.join(str(c) for c in ctd.values())
            
            #pseudo amino acid composition descriptors  total 30
            protein=PyPro()   #protein object
            protein.ReadProteinSequence(seq)
            paac = protein.GetPAAC(lamda=10,weight=0.05)
            k_list3=','.join(c for c in paac.keys())
            v_list3=','.join(str(c) for c in paac.values())
            #all protein descriptors total 2049
            allp=protein.GetALL()
            k_list4=','.join(c for c in allp.keys())
            v_list4=','.join(str(c) for c in allp.values())
        except:
            print "\nerror while calculation of features in:\t",f
            continue
        #fp.write(str(key+"_"+f.replace(".pdb",""))+","+str((line[0])[:locSlash(line[0])-1])+","+str((line[1])[:locSlash(line[1])-1])+","+str((line[2])[:(locSlash(line[2])-1)])+","+str(f15)+","+str(v_list1)+","+str(v_list2)+"\n")
        fp.write(str(key+"_"+f.replace(".pdb",""))+","+str((line[0])[:locSlash(line[0])-1])+","+str((line[1])[:locSlash(line[1])-1])+","+str((line[2])[:(locSlash(line[2])-1)])+","+str(f15)+","+str(v_list1)+","+str(v_list2)+","+str(v_list3)+","+str(v_list4)+"\n")
        fileTime=time.time()-time1
        print "\ntime for file\t"+f+"is\t",fileTime
    folderTime=time.time()-time1
    print "\ncompletion time for folder\t"+key+"is\t",folderTime
fp.close()


print "done\n"
print "total time \t",round(time.time()-startTime,2)
