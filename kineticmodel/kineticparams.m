% Vmax: µmol/(min*mg)
% Km: mM
% Kcat: 1/s
% E0: mM

%Create a struct for all kinetic params
function params=kineticparams(Vmax)
        
    params=struct();

    % Constants
    params.constant_CoA = 0.5;%mM
    params.NH4_MALratio = 6.735/9.158;%NH4/malate
    params.NH4_GLYratio = 6.735/0.005;%NH4/glycine

    %% Nutrients exchange
    %Reaction: id=SourceMAL
    %MAL_cytosol => MAL_Bacteroid
    params.SourceMAL_KmMAL = 1;
    params.SourceMAL_Vmax= Vmax.MALOut;
    params.SourceMAL_Keq=1;

    %Reaction: id=SourceSUC
    %SUC_cytosol => SUC_Bacteroid
    params.SourceSUC_KmSUC = 2.5;
    params.SourceSUC_Vmax= Vmax.SUCOut;

    % Reaction: id=SourceFUM
    % FUM_cytosol => FUM_Bacteroid
    params.SourceFUM_KmFUM = 3.91;
    params.SourceFUM_Vmax= Vmax.FUMOut;

    % Reaction: id=ToNH4
    % B=>NH4
    params.ToNH4_KmB=0.5;
    params.ToNH4_Vmax=Vmax.ToNH4;

    %Reaction: id=AllantoinOut
    %Allantoin_nodule => Allantoin_root
    params.ALTNOut_Km = 0.0762;
    params.ALTNOut_Vmax = Vmax.ALTNOut;

    %% Glycolysis
    % Reaction: id = HKI, name = Hexokinase
    % GLC => G6P
    params.HKI_KmGLC=0.075;
    params.HKI_KmG6P=0.5;
    params.HKI_Keq=850;
    params.HKI_Vmax=Vmax.HKI;
    
    % Reaction: id = PGI, name = PGI
    % G6P <=> F6P
    params.PGI_KmG6P=0.27;
    params.PGI_KmF6P=0.48;
    params.PGI_Keq=0.276;
    params.PGI_KiG6P_inh6PG=0.013;
    params.PGI_Vmax= Vmax.PGI;
    
    % Reaction: id = PFK, name = PFK 
    % ATP+F6P=>ADP+FDP
    params.PFK_KmF6P=1.5;
    params.PFK_Vmax=Vmax.PFK;
    
    % Reaction: id = FBP 
    % FDP => F6P+PO4
    params.FBP_KmFDP = 0.25;
    params.FBP_Vmax=Vmax.FBP;
    
    % Reaction: id = FBA, name = Aldolase
    % FDP <=> DAP+GAP
    params.FBA_KmFDP=0.0167;
    params.FBA_KmGAP=0.0605375;
    params.FBA_KmDAP=0.104375;
    params.FBA_Keq=0.19;
    params.FBA_KiFDP_inhR5P=2.2;
    params.FBA_Vmax=Vmax.FBA;
    
    % Reaction: id = TPI, name = TPI 
    % DAP <=> GAP
    params.TPI_KmDAP=0.812;
    params.TPI_KmGAP=0.245;
    params.TPI_Ki3PG = 0.4;
    params.TPI_KiPEP = 0.661;
    Kcat_DAP=1080;
    Kcat_GAP=6170;
    params.TPI_Keq=Function_Keq1(Kcat_DAP,Kcat_GAP,params.TPI_KmDAP,params.TPI_KmGAP);
    params.TPI_Vmax = Vmax.TPI;
    
    %Reaction: id = GDH, name = GDH
    %GAP + NAD + PO4 <=> BPG + NADH 
    params.GDH_KmGAP = 0.074;
    params.GDH_KmPO4 = 9;
    params.GDH_KmBPG = 0.036;
    params.GDH_Keq=20;
    params.GDH_Vmax= Vmax.GDH;
    
    % Reaction: id = PGK 
    % ADP+BPG <=> ATP+PGA3
    params.PGK_KmBPG = 0.1;
    params.PGK_KmPGA3 = 0.146;
    params.PGK_Keq = 100;
    params.PGK_Vmax= Vmax.PGK;
    
    % Reaction: id = GPM 
    % PGA3 <=> PGA2
    params.GPM_KmPGA3 = 0.1;
    params.GPM_KmPGA2 = 0.369;
    GPM_Kcat_PGA3=3.01;
    GPM_Kcat_PGA2=4.01;
    params.GPM_Keq=Function_Keq1(GPM_Kcat_PGA3,GPM_Kcat_PGA2,params.GPM_KmPGA3,params.GPM_KmPGA2);
    params.GPM_Vmax= Vmax.GPM;
    
    % Reaction: id = ENO
    % PGA2 <=> PEP
    params.ENO_KmPGA2 = 0.19;
    params.ENO_KmPEP = 0.534;
    params.ENO_Keq = 3;
    params.ENO_Vmax= Vmax.ENO;
    
    % Reaction: id = PYK
    % ADP+PEP => ATP+PYR
    params.PYK_KmPEP = 0.15;
    params.PYK_Vmax=Vmax.PYK;
    
    % Reaction: id = PDH 
    % COA+NAD+PYR => ACCOA+NADH+HCO3
    params.PDH_KmPYR = 0.4;
    params.PDH_KmCoA = 0.00061;
    params.PDH_KiPYR=31;
    params.PDH_Vmax=Vmax.PDH;
    
    %% PP pathway
    % Reaction: id = ZWF 
    % G6P+NADP => GL6P+NADPH
    params.ZWF_KmG6P = 0.1073;
    params.ZWF_KiG6P_inhNH4=0.00226;
    params.ZWF_KmGL6P=0.329;
    params.ZWF_Keq=1000;
    params.ZWF_KdG6P=0.192;
    params.ZWF_KdGL6P=0.02;
    params.ZWF_Vmax=Vmax.ZWF;

    % Reaction: id = PGL 
    % GL6P => PGN(6PG)
    params.PGL_KmGL6P = 0.83;
    params.PGL_KmPGN=1;
    params.PGL_Keq=42.7;
    params.PGL_KiG6P=2;
    params.PGL_Vmax=Vmax.PGL;

    %Reaction: id = GND 
    %NADP+PGN(6PG)=>NADPH+RB5P+HCO3
    params.GND_KmPGN = 0.01;
    params.GND_KmRB5P=45.2;
    params.GND_KdRB5P=0.044;
    params.GND_Keq=10;
    params.GND_Vmax=Vmax.GND;

    % Reaction: id = RPI 
    % RB5P <=> R5P
    params.RPI_KmRB5P = 0.01;
    params.RPI_KmR5P = 0.88;
    params.RPI_Keq = 0.333;
    params.RPI_KiR5P_inhE4P=0.1;
    params.RPI_KiRB5P_inhE4P=0.14;
    params.RPI_KiR5P_inhGAP=0.5;
    params.RPI_KiR5P_inhPGA3=1.2;
    params.RPI_KiR5P_InhPGN=7;
    params.RPI_Vmax=Vmax.RPI;

    % Reaction: id = RPE 
    % X5P <=> RB5P 
    params.RPE_KmX5P = 0.067;
    params.RPE_KmRB5P = 2.5;
    params.RPE_Keq = 1;
    params.RPE_Vmax=Vmax.RPE;

    % Reaction: id = TAL, name = F6P_E4P_TAL
    % F6P+E4P <=> S7P+GAP  
    params.TAL_Keq = 26.6266;
    params.TAL_KmF6P = 1.2;
    params.TAL_KmE4P = 0.09;    
    params.TAL_KmS7P = 0.285;  
    params.TAL_KmGAP = 0.038;    
    params.TAL_Vmax = Vmax.TAL;

    % Reaction: id = TKT1, name = GAP_E4P_TKT
    % GAP+F6P <=> E4P+X5P
    params.TKT1_Keq = 1;
    params.TKT1_KmP5P = 0.616;
    params.TKT2_KmP5P = 0.118;
    params.TKT1_KmGAP = 0.2727;
    params.TKT1_KmF6P = 0.5443;
    TKT1_Vmax = 100;
    params.TKT1_Vmax = Vmax.TKT1; 


    % Reaction: id = TKT2, name = GAP_S7P_TKT	
    % GAP+S7P <=> R5P+X5P
    params.TKT2_Keq = 1000;
    params.TKT2_KmS7P = 0.01576;
    params.TKT2_KmGAP = 0.09078;
    TKT2_Vmax = 100;%mM/s
    params.TKT2_Vmax = Vmax.TKT2;   
    
    %% PEP to OAA    
    % Reaction: id = PEPC
    % PEP+HCO3 => OAA+Pi
    params.PEPC_KmPEP = 0.09;
    params.PEPC_KiMAL = 37.1;
    params.PEPC_KiGLU = 11.5;
    params.PEPC_KiASP = 12.6;
    params.PEPC_Vmax = Vmax.PEPC;
    
    %% TCA cycle
    % Reaction: id = GLT 
    % ACCOA+OAA <=> CIT+COA
    params.GLT_KmACCOA = 0.12;
    params.GLT_KmOAA = 0.02;
    params.GLT_KmCIT = 1.16;
    params.GLT_KmCOA = 0.0001;
    params.GLT_KiOAA_inhOAA=0.033;
    params.GLT_KiOAA_inh2kg=0.76;
    params.GLT_Keq = 8300;
    params.GLT_Vmax = Vmax.GLT;

    % Reaction: id = ACN1
    % CIT <=> ACO
    params.ACN1_KmCIT = 0.1;
    params.ACN1_KmACO = 0.3;
    ACN1_Kcat_CIT = 5.3;
    ACN1_Kcat_ACO = 5.2;
    params.ACN1_Keq = Function_Keq1(ACN1_Kcat_CIT,ACN1_Kcat_ACO,params.ACN1_KmCIT,params.ACN1_KmACO);    
    params.ACN1_Vmax = Vmax.ACN1;

    % Reaction: id = ACN2
    % ACO <=> ICIT
    params.ACN2_KmACO = 0.2;  
    params.ACN2_KmICIT = 0.58;
    ACN2_Kcat_ACO = 5.2;
    ACN2_Kcat_ICIT = 1.1;
    params.ACN2_Keq = Function_Keq1(ACN2_Kcat_ACO,ACN2_Kcat_ICIT,params.ACN2_KmACO,params.ACN2_KmICIT);
    params.ACN2_Vmax = Vmax.ACN2;

    % Reaction: id = ICD 
    % ICIT+NADP+ => AKG+CO2+NADPH+H+
    params.ICD_KmICIT = 0.01;
    params.ICD_KmAKG = 0.5;
    params.ICD_Ki2kg_inhICIT=0.012;
    ICD_Kcat_ICIT = 51.8;
    ICD_Kcat_AKG = 172.667;
    params.ICD_Keq = Function_Keq1(ICD_Kcat_ICIT,ICD_Kcat_AKG,params.ICD_KmICIT,params.ICD_KmAKG);
    params.ICD_Vmax = Vmax.ICD;

    % Reaction: id=LPD
    % AKG+COA+NAD => SUCCOA+HCO3-+NADH
    params.LPD_KmAKG = 10;
    params.LPD_Vmax = Vmax.LPD;

    % Reaction: id = SK 
    % ADP+SUCCOA+PO4 <=> ATP+SUC+COA
    params.SK_KmSUCCOA = 0.41;
    params.SK_KmPO4 = 0.72;
    params.SK_KmSUC = 1.5;
    params.SK_KmCOA = 0.1;
    SK_Kcat_SUCCOA=22;
    SK_Kcat_SUC=29;
    params.SK_Keq=1000;
    params.SK_Vmax= Vmax.SK;

    % Reaction: id=SDH (E.coli)
    % Q+SUC <=> FUM+QH2
    params.SDH_KmSUC = 0.005;
    params.SDH_KmFUM = 0.02;
    params.SDH_KiSUC_inhOAA = 0.00007;
    SDH_Kcat_SUC=15;
    SDH_Kcat_FUM=32;
    params.SDH_Keq=1000;
    params.SDH_Vmax = Vmax.SDH;

    % Reaction: id=FUMA
    % FUM <=> MAL
    params.FUMA_KmFUM = 0.1;
    params.FUMA_KmMAL = 0.3;
    params.KiFUM_inhPEP = 0.7;
    params.KiFUM_inhPYR = 1.6;
    params.KiFUM_inhCIT=0.5;
    params.KiFUM_inhAkg=0.8;
    params.KiFUM_inhOAA=1.2;
    FUMA_Kcat_FUM = 21;
    FUMA_Kcat_MAL = 15.3;
    params.FUMA_Keq = 310;
    params.FUMA_Vmax = Vmax.FUMA;

    % Reaction: id=MDH 
    % NAD+OAA <=> MAL+NADH+H
    params.MDH_KmOAA = 0.018;
    params.MDH_KmMAL = 2.6;
    MDH_KcatOAA = 960;
    MDH_KcatMAL = 250;    
    params.MDH_Keq = 0.5;
    params.MDH_Vmax = Vmax.MDH;
     
    %% Glyoxylate shunt
    % Reaction: id=ICL 
    %ICIT => GLX+SUC
    params.ICL_KmICIT = 8;
    params.ICL_Vmax= Vmax.ICL;

    % Reaction: id=MALS 
    % ACCOA+GLX <=> MAL+COA
    params.MALS_KmMAL=8;
    params.MALS_KmACCOA=0.8;
    params.MALS_KmGLX=30;
    params.MALS_KiGLX_inhPYR=0.54;
    params.MALS_KiGLX_inhOAA=1.5;
    params.MALS_Keq = 100;
    params.MALS_Vmax= Vmax.MALS;


    %% GS/GOGAT
    % Reaction: id=GD
    % GLU+H2O+NAD+ <=> AKG+NH4+NADH+H+
    params.GD_KmAKG = 2.1;
    params.GD_KmNH4 = 15.8;
    params.GD_KmGLU = 1.67;
    GD_Vmax = 0.96;
    GD_Kcat_AKG=165.3;
    GD_Kcat_GLU=121.2;
    params.GD_Ki_2kg_inhGLU = 20;
    params.GD_Keq = Function_Keq12(GD_Kcat_GLU,GD_Kcat_AKG, params.GD_KmGLU,params.GD_KmAKG,params.GD_KmNH4);
    params.GD_Vmax = Vmax.GD;

    % Reaction: id=GS
    % GLU+NH4+ATP => GLN+PO4+ADP
    params.GS_KmGLU = 9;
    params.GS_KmNH4 = 5.2;
    params.GS_KmGLN = 1.3;
    params.GS_Ki_GLU_inhGLN = 6.6;
    params.GS_Ki_NH4_inhGLN = 7.4;
    params.GS_Ki_GLN_InhALN = 0.03;
    params.GS_Ki_GLN_InhGLY = 0.1;    
    params.GS_Keq = 1000;
    GS_Kcat_GLU = 1520;
    GS_Kcat_NH4 = 13.5;
    params.GS_Vmax = Vmax.GS;

    % Reaction: id=GOGAT
    % GLN+AKG+NADH+H+ =>2GLU+NAD+
    params.GOGAT_KmAKG = 0.039;
    params.GOGAT_KmGLN = 4;
    params.GOGAT_KiAKG_inhGLU = 0.7;
    params.GOGAT_KiAKG_inhOAA = 5;
    params.GOGAT_KiAKG_inhASP = 2.7;
    GOGAT_Kcat = 14000;
    params.GOGAT_Vmax = Vmax.GOGAT;

    %% Alanine metabolism
    % Reaction: id=AGT
    % alanine + glyoxylate <=> pyruvate + glycine
    params.AGT_KmALN = 0.31;
    params.AGT_KmGLX = 0.23;
    params.AGT_KmPYR = 0.21;
    params.AGT_KmGLY = 22;
    params.AGT_KiGLX_inhPYR=2.3;
    params.AGT_KiALN_inhPYR_Kis=15.8;
    params.AGT_KiALN_inhPYR_Kii=22.8;
    AGT_Kcat_ALN = 45;
    AGT_Kcat_GLY = 10.33;
    params.AGT_Keq = Function_Keq2(AGT_Kcat_ALN, AGT_Kcat_GLY, params.AGT_KmALN, params.AGT_KmGLX,params.AGT_KmPYR, params.AGT_KmGLY);
    params.AGT_Vmax = Vmax.AGT;

    % Reaction: id=ALTSS 
    % PYR+GLU <=> ALN+AKG
    params.ALTSS_KmPYR = 0.11;
    params.ALTSS_KmGLU = 0.68;  
    params.ALTSS_KmALN = 0.4;   
    params.ALTSS_KmAKG = 0.25; 
    ALTSS_KcatALN = 4.98;
    ALTSS_KcatPYR = 6.86;
    params.ALTSS_Keq = Function_Keq2(ALTSS_KcatPYR, ALTSS_KcatALN, params.ALTSS_KmPYR, params.ALTSS_KmGLU, params.ALTSS_KmALN, params.ALTSS_KmAKG);
    params.ALTSS_Vmax = Vmax.ALTSS;

    %% Glycine and serine
    % Reaction: id=PGDH 
    % 3PG+NAD+ <=> 3PHP+NADH+H+
    params.PGDH_Km3PG = 0.29;
    params.PGDH_Km3PHP = 0.1;
    params.PGDH_Keq = 10;   
    PGDH_Vmax = 0.365158077;
    params.PGDH_Vmax = Vmax.PGDH;

    % Reaction: id=PSTS 
    % 3PHP+GLU <=> serine+AKG
    params.PSTS_Km3PHP = 0.08;
    params.PSTS_KmGLU = 5;
    params.PSTS_KmSerine = 0.225;
    params.PSTS_KmAKG = 0.0651;
    PSTS_Kcat_3PHP = 1.75;
    PSTS_Kcat_AKG = 0.63;
    params.PSTS_Keq = Function_Keq2(PSTS_Kcat_3PHP, PSTS_Kcat_AKG, params.PSTS_Km3PHP, params.PSTS_KmGLU, params.PSTS_KmSerine, params.PSTS_KmAKG);
    params.PSTS_Vmax = Vmax.newPSTS;

    % Reaction: id=SHMT
    % Serine+THF <=> GLY+CH2THF+H2O
    params.SHMT_KmSerine = 2.5;
    params.SHMT_KmTHF = 0.25;
    params.SHMT_KmGLY = 0.66;
    params.SHMT_KmCH2THF = 0.98;
    params.SHMT_KiSRN_inhGLY=3;
    params.SHMT_KiTHF_inhCH2THF=2.9;
    SHMT_Kcat_THF = 14.17;
    SHMT_Kcat_GLY = 0.188;
    params.SHMT_Keq = 10;
    params.SHMT_Vmax = Vmax.SHMT;

    % Reaction: id=Mthfr
    % CH2THF+NAD+ => CHTHF+NADH+H+
    params.MTHFR_KmCH2THF = 0.0004;
    MTHFR_Kcat=2.2;
    params.MTHFR_KiCH2THF=0.061;
    params.MTHFR_Vmax=Vmax.MTHFR;

    % Reaction: id=Mthfd
    % CHTHF+H2O => CHOTHF
    params.MTHFD_KmCHTHF=0.026;
    MTHFD_Kcat=40.9;
    params.MTHFD_Vmax=Vmax.MTHFD;

    % Reaction: id=GLYDC
    % glycine + THF + NAD => CH2THF + NADH + CO2 + NH3
    params.GLYDC_KmGLY = 6;
    params.GLYDC_KmTHF = 7;
    params.GLYDC_KiGLY_inhSRN = 4;
    GLYDC_Kcat = 100; 
    params.GLYDC_Vmax = Vmax.GLYDC;

    %% de novo Purine synthesis
    % Reaction: id=PRS
    % R5P+ATP => PRPP+AMP
    params.PRS_KmR5P = 0.11;
    params.PRS_KmPRPP = 0.5;
    params.PRS_KisR5P_inhPRPP=0.82;
    params.PRS_KiiR5P_inhPRPP=1.7; 
    params.PRS_Keq = 0.1;
    params.PRS_Vmax= Vmax.PRS;

    % Reaction: id=PRAT 
    %GLN+PRPP+H2O => PRA+Ppi+GLU
    params.PRAT_KmGLN=18;
    params.PRAT_KmPRPP=0.4;
    params.PRAT_KmPRA=0.8;
    params.PRAT_KmGLU=20;
    params.PRAT_KiGLN_inhNH4 = 16;
    params.PRAT_KiGLN_inhGLU = 30;
    params.PRAT_KiPRPP_inhIMP = 1.7;
    params.PRAT_KiPRPP_inhXMP = 1.2;
    params.PRAT_Keq=1e3;
    params.PRAT_Vmax= Vmax.PRAT;

    % Reaction: id=GARS 
    %PRA + ATP + glycine <=> ADP + PO4 + beta-GAR
    params.GARS_KmPRA = 7;
    params.GARS_KmGLY = 0.27;
    params.GARS_KmGAR = 0.3;
    params.GARS_KmPO4 = 0.54;
    GARS_Kcat_GLY = 7; 
    GARS_Kcat_GAR = 2.3;
    params.GARS_Keq = Function_Keq2(GARS_Kcat_GLY, GARS_Kcat_GAR, params.GARS_KmPRA, params.GARS_KmGLY, params.GARS_KmGAR, params.GARS_KmPO4);
    params.GARS_Vmax = Vmax.GARS;

    % Reaction: id=GARTF
    % CHOTHF10 + beta-GAR => Tthdf + FGAR
    params.GARTF_KmCHOTHF = 0.2;
    params.GARTF_KmGAR = 19.2;
    params.GARTF_CHOTHF_inhTHF=0.012;
    GARTF_Kcat = 16.1;
    params.GARTF_Vmax = Vmax.GARTF;

    % Reaction: id=FGAMS
    % FGAR+GLN+ATP+H2O => GLU+FGAM+PO4+ADP
    params.FGAMS_KmFGAR = 4.68;
    params.FGAMS_KmGLN = 0.03;
    params.FGAMS_KiGLN_inhGLU = 1.6;
    params.FGAMS_KisFGAR_inhGLU=13;
    params.FGAMS_KiiFGAR_inhGLU=52;
    FGAMS_Kcat = 5;
    params.FGAMS_Vmax = Vmax.FGAMS;

    % Reaction: id=AIRS
    % FGAM + ATP => AIR+PO4
    params.AIRS_KmFGAM = 2.7;
    params.AIRS_KiiFGAM_inhAIR = 0.0647;
    params.AIRS_KisFGAM_inhAIR = 0.0654;
    params.AIRS_Vmax = Vmax.AIRS;
     
    % Reaction: id=CAIRS
    % AIR+CO2<=>CAIR
    params.CAIRS_KmAIR = 0.76;
    params.CAIRS_KmCAIR = 1;
    CAIRS_Kcat_AIR=40;
    CAIR_Kcat_CAIR=20;
    params.CAIRS_Keq = Function_Keq1(CAIRS_Kcat_AIR,CAIR_Kcat_CAIR,params.CAIRS_KmAIR,params.CAIRS_KmCAIR);
    params.CAIRS_Vmax = Vmax.CAIRS;

    % Reaction: id=SAICARS
    % ATP+CAIR+ASP <=> ADP+PO4+SAICAR
    params.SS_KmCAIR = 0.5;
    params.SS_KmASP = 1.4;
    params.SS_KmSAICAR = 1.5;
    params.SS_KmPO4 = 0.1;
    params.SS_Keq = 10;
    params.SS_KiCAIR_inhIMP=9.1;
    SS_Kcat_CAIR=240;
    SS_MW = 76000;
    params.SS_Vmax = Vmax.SS;

    % Reaction: id=ADSL 
    % SAICAR <=> FUM+AICAR
    params.ADSL_KmSAICAR = 1;
    params.ADSL_KmFUM = 0.52;  
    params.ADSL_KmAICAR = 0.009;
    ADSL_Kcat_SAICAR=337;
    ADSL_Kcat_FUM=2.9;
    params.ADSL_Keq = Function_Keq12(ADSL_Kcat_SAICAR,ADSL_Kcat_FUM,params.ADSL_KmSAICAR,params.ADSL_KmFUM,params.ADSL_KmAICAR);
    ADSL_Vmax = 1.29;
    ADSL_MW = 191000;
    params.ADSL_Vmax = Vmax.ADSL;

    % Reaction: id=AICARTF 
    % CHOTHF+AICAR => FAICAR+TTHDF
    params.AICARTF_KmCHOTHF = 0.2;
    params.AICARTF_KmAICAR = 0.31;
    AICARTF_Kcat = 4.97;
    params.AICARTF_Vmax= Vmax.AICARTF;

    % Reaction: id=IMPCH 
    % FAICAR <=> IMP+H2O
    params.IMPCH_KmFAICAR = 0.78;
    params.IMPCH_KmIMP = 0.5;
    IMPCH_Kcat_FAICAR=1.32;
    IMPCH_Kcat_IMP=0.5;
    params.IMPCH_Keq = Function_Keq1(IMPCH_Kcat_FAICAR,IMPCH_Kcat_IMP,params.IMPCH_KmFAICAR,params.IMPCH_KmIMP);
    params.IMPCH_Vmax = Vmax.IMPCH;

    % Reaction: id=IMPDH 
    % IMP+NAD+H2O => XMP+NADH+H+
    params.IMPDH_KmIMP = 0.5;
    IMPDH_Kcat= 6.1;
    params.IMPDH_Vmax = Vmax.IMPDH;

    % Reaction: id=NT 
    % XMP + H2O => Xanthosine + PO4
    params.NT_KmXMP = 0.77;
    NT_Vmax = 5;    
    params.NT_Vmax = Vmax.NT;

    % Reaction: id=PNP 
    % Xanthosine+PO4 => Xanthine+R1P
    params.PNP_KmXao = 0.51;
    params.PNP_KmPO4 = 0.76;
    PNP_Vmax=2.6;%U/mg
    params.PNP_Vmax = Vmax.PNP;

    % Reaction: id=XOR 
    % Xanthine+NAD+H2O => Urate+NADH+H+
    params.XOR_KmXan = 0.005;
    params.XOR_KiXAN_inhUrate = 0.18;
    XOR_Vmax = 0.02;
    params.XOR_Vmax=Vmax.XOR;

    % Reaction: id=UOD 
    % Urate+O2+H2O => HDS5+H2O2
    params.UOD_KmURATE = 0.01;
    params.UOD_KmO2 = 0.031;
    params.UOD_Ki_Urate_inhXAN = 0.01;
    UOD_Vmax = 9.58;      
    params.UOD_Vmax = Vmax.UOD;

    % Reaction: id=HIUHS 
    % HDS5+H2O = OHCU
    params.HIUHS_KmHIUH = 0.015;
    HIUHS_Kcat = 11;
    params.HIUHS_Vmax = Vmax.HIUHS;

    % Reaction: id=OHCUD 
    % OHCU => Allantoin+CO2
    params.OHCUD_KmOHCU = 0.151;
    OHCUD_Kcat = 122;
    params.OHCUD_Vmax = Vmax.OHCUD;

    %% Aspartate metabolism
    % Reaction: id=ATS
    % OAA+GLU <=> ASP+AKG
    params.ATS_KmOAA = 0.02;
    params.ATS_KmGLU = 12;   
    params.ATS_KmASP = 0.05;   
    params.ATS_KmAKG = 0.2;
    params.ATS_KiAKG_inhGLU=10;
    params.ATS_KiASP_inhGLU=12;
    params.ATS_KiASP_inhOAA=0.027;
    params.ATS_KiAKG_inhOAA=0.056;
    params.ATS_KiOAA_inhOAA=0.2;
    ATS_VmaxOAA = 320;
    ATS_VmaxAKG = 160;
    ATS_MW=96000;
    params.ATS_Keq = Function_Keq2(ATS_VmaxOAA, ATS_VmaxAKG, params.ATS_KmOAA, params.ATS_KmGLU, params.ATS_KmASP, params.ATS_KmAKG);
    params.ATS_Vmax = Vmax.SS;

  
end

% A<=>B
function z=Function_Keq1(Kcatf,Kcatb,KmA,KmB)
z = (Kcatf*KmB)/(Kcatb*KmA);
end
% A<=>C+D
function z=Function_Keq12(Kcatf,Kcatb,KmA,KmC,KmD)
z = (Kcatf*KmC*KmD)/((Kcatb)^2*KmA);
end
% A+B<=>C
function z=Function_Keq21(Kcatf,Kcatb,KmA,KmB,KmC)
z = ((Kcatf^2)*KmC)/(Kcatb*KmA*KmB);
end
% A+B<=>C+D
function z=Function_Keq2(Kcatf, Kcatb, KmA, KmB, KmC, KmD)
z = Kcatf^2/Kcatb^2 * (KmC*KmD)/(KmA*KmB);
end

function z=Enzyme_Concentration(E, MW)
%assumed E:mmol/L, MW:g/mol
z = (E/1000)*(MW*1000);%mg/L
end

function z=Vmax_UnitConvertion(Vmax, proteinconc)
%assume Vmax:U/mg, proteinconc:mg/L
z = Vmax*proteinconc/1000/60;%mM/s
end

function z=Vmax_Calculation(Kcat, E)
%assume Kcat:1/s, E:mM
z = Kcat*E;
end


