%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script is to calculate flux values.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Flux = JRRModel_Flux(t,x,Vmax)

    %load params
    params = kineticparams(Vmax);    

    %metabolites
    GLC = x(1); 
    G6P = x(2);
    F6P = x(3);
    FDP = x(4);
    GAP = x(5);
    DAP = x(6);
    BPG = x(7);
    PGA3 = x(8);
    PGA2 = x(9);
    PEP = x(10);
    PYR = x(11);
    %PP pathway
    GL6P = x(12);
    PGN = x(13);
    RB5P = x(14);
    R5P = x(15);
    X5P = x(16);
    E4P = x(17);
    S7P = x(18);
    % % TCA cycle
    OAA = x(19);
    MAL = x(20); 
    ACCOA = x(21);
    CIT = x(22);
    ACO = x(23);
    ICIT = x(24);
    AKG = x(25);
    GLX = x(26);

    SUCCOA = x(27);
    SUC = x(28);
    FUM = x(29);

    B=x(30);

    % GS/GOGAT
    NH4 = x(31);
    GLU = x(32);
    GLN = x(33);

    % Purine
    PRPP = x(34);

    %glycine and serine
    PHP3 = x(35);
    SRN = x(36);
    GLY = x(37);
    PRA = x(38);
    GAR = x(39);
    FGAR = x(40);
  
    ALN = x(41);
    THF = x(42);
    CH2THF = x(43);
    CHTHF = x(44);
    CHOTHF = x(45);
     
    FGAM = x(46);
    AIR = x(47);
    CAIR = x(48);
    ASP = x(49);
    SAICAR = x(50);
    AICAR = x(51);
    FAICAR = x(52);
    IMP = x(53);
    XMP = x(54);
    XAO = x(55);
    XAN = x(56);
    URATE = x(57);
    HIUH = x(58);
    OHCU = x(59);
    ALTN = x(60);

    %% Metabolite transfer
    MALOut = rxns.MALOut(params,MAL);
    SUCOut = rxns.SUCOut(params,SUC);
    FUMOut = rxns.FUMOut(params,FUM);
    ToNH4 = params.ToNH4_Vmax * B / (params.ToNH4_KmB + B);

    % Glycolysis
    HKI = rxns.reaction_HKI(params, GLC,G6P);
    PGI = rxns.reaction_PGI(params, G6P, F6P);
    PFK = rxns.reaction_PFK(params, F6P);
    FBP = rxns.reaction_FBP(params, FDP);
    FBA = rxns.reaction_FBA(params,FDP,GAP,DAP);
    TPI = rxns.reaction_TPI(params,DAP,GAP,PGA3,PEP);
    GDH = rxns.reaction_GDH(params,GAP,BPG);
    PGK = rxns.reaction_PGK(params,BPG,PGA3);
    GPM = rxns.reaction_GPM(params,PGA3,PGA2);
    ENO = rxns.reaction_ENO(params,PGA2,PEP);
    PYK = rxns.reaction_PYK(params, PEP);
    PDH = rxns.reaction_PDH(params, PYR);

    % % PP pathway
    ZWF = rxns.reaction_ZWF(params, G6P,GL6P);
    PGL = rxns.PGL(params,GL6P, PGN, G6P);
    GND = rxns.reaction_GND(params, PGN,RB5P);
    RPE = rxns.reaction_RPE(params,X5P,RB5P);  
    RPI = rxns. reaction_RPI(params,RB5P,R5P,E4P,GAP,PGA3,PGN);
    TKT1 = rxns.TKT1(params,GAP, F6P, S7P, X5P, E4P, R5P);
    TKT2 = rxns.TKT2(params,GAP, F6P, S7P, X5P, E4P, R5P);
    TAL = rxns.TAL(params, F6P, E4P, S7P, GAP);
    % 
    %% Anaplerotic reactions
    PEPC = rxns.PEPC(params,PEP,MAL);
    
    %% TCA cycle
    MDH = rxns.reaction_MDH(params,OAA,MAL);
    GLT = rxns.reaction_GLT(params, ACCOA, OAA,CIT,AKG);
    ACN1 = rxns.reaction_ACN1(params, CIT,ACO);
    ACN2 = rxns.reaction_ACN2(params, ACO,ICIT);
    ICDH = rxns.reaction_ICD(params, ICIT, AKG);
    LPD = rxns.reaction_LPD(params, AKG);
     
    SK = rxns.reaction_SK(params, SUCCOA);
    FUMA = rxns.reaction_FUMA(params, FUM,MAL,PEP,PYR,CIT,AKG,OAA);
    SDH = rxns.reaction_SDH(params, SUC, FUM, OAA);
     
    ICL = rxns.reaction_ICL(params, ICIT);
    MALS = rxns.reaction_MALS(params,ACCOA,GLX,MAL,PYR,OAA);
    
    %% GS/GOGAT
    GS = rxns.reaction_GS(params, GLU, NH4,GLN,ALN,GLY);
    GOGAT = rxns.reaction_GOGAT(params, GLN, AKG, GLU, OAA, ASP);
    GD = rxns.reaction_GD(params,GLU,AKG,NH4);
    
    %% Alanine metbaolism
    ALTSS = rxns.reaction_ALTSS(params,PYR,GLU,ALN,AKG);
    AGT = rxns.reaction_AGT(params,ALN,GLX,PYR,GLY);
 
    %% glycine and serine
    PGDH = rxns.reaction_PGDH(params,PGA3, PHP3);
    newPSTS = rxns.reaction_newPSTS(params,PHP3,GLU,SRN,AKG);
    SHMT = rxns.reaction_SHMT(params,SRN,THF,GLY,CH2THF);
    MTHFR = rxns.reaction_MTHFR(params,CH2THF);
    MTHFD = rxns.reaction_MTHFD(params,CHTHF);
    GLYDC = rxns.reaction_GLYDC(params,GLY,THF);
    GLYsink = MALOut*(0.005/9.158);
     
    %% de novo purine synthesis
    PRS = rxns.reaction_PRS(params, R5P, PRPP);
    PRAT = rxns.reaction_PRAT(params, PRPP, GLN, PRA, GLU, NH4,IMP);
    GARS = rxns.reaction_GARS(params, PRA, GLY, GAR);
    GARTF = rxns.reaction_GARTF(params, GAR,CHOTHF,THF);
    FGAMS = rxns.reaction_FGAMS(params, FGAR, GLN, GLU);
    AIRS = rxns.reaction_AIRS(params, FGAM);
    CAIRS = rxns.reaction_CAIRS(params, AIR);
    SS = rxns.reaction_SS(params, CAIR,ASP);
    ADSL = rxns.reaction_ADSL(params, SAICAR,FUM,AICAR);
    AICARTF = rxns.reaction_AICARTF(params,AICAR,CHOTHF);
    IMPCH = rxns.reaction_IMPCH(params, FAICAR,IMP);
    IMPDH = rxns.reaction_IMPDH(params, IMP);
    NT = rxns.reaction_NT(params, XMP);
    PNP = rxns.reaction_PNP(params, XAO);
    XOR = rxns.reaction_XOR(params, XAN, URATE);
    UOD = rxns.reaction_UOD(params, URATE, XAN);
    HIUHS = rxns.reaction_HIUHS(params, HIUH);
    OHCUD = rxns.reaction_OHCUD(params, OHCU);
    ALTNOut = rxns.ALTNOut(params,ALTN);
 
    %% Aspartate metabolism
    ATS = rxns.reaction_ATS(params,OAA,GLU,ASP,AKG);


    Flux = [HKI ...
        PGI ...
        PFK ...
        FBP ...
        FBA ...
        TPI ...
        GDH ...
        PGK ...
        GPM ...
        ENO ...
        PYK ...
        PDH ...
        ZWF ...
        PGL ...
        GND ...
        RPE ...
        RPI ...
        TKT1 ...
        TKT2 ...
        TAL ...
        PRS ...
        PEPC ...
        MDH ...
        MALOut ...
        GLT ... 
        ACN1 ...
        ACN2 ... 
        ICL ...
        MALS ...
        ICDH ...
        LPD ...
        SK ...
        SUCOut ...
        SDH ...
        FUMOut ...
        FUMA ...
        ToNH4 ...
        GS ...
        GOGAT ...
        GD ...
        PRAT ...
        PGDH ...
        newPSTS ...
        SHMT ...
        GARS ...
        GARTF ...
        FGAMS ...
        ALTSS ...
        AGT ...
        MTHFR ...
        MTHFD ...
        GLYDC ...
        AIRS ...
        CAIRS ...
        SS ...
        ATS ...
        ADSL ...
        AICARTF ...
        IMPCH ...
        IMPDH ...
        GLYsink ...
        NT ...
        PNP ...
        XOR ...
        UOD ...
        HIUHS ...
        OHCUD ...
        ALTNOut 
        ];


end

