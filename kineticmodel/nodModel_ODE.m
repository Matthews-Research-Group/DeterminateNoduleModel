%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script includes the ODE system. 
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = nodModel_ODE(t,x,Vmax)
    
    Flux = nodModel_Flux(t,x,Vmax);

    v_HKI = Flux(1);
    v_PGI = Flux(2);
    v_PFK = Flux(3);
    v_FBP = Flux(4);
    v_FBA = Flux(5);
    v_TPI = Flux(6);
    v_GDH = Flux(7);
    v_PGK = Flux(8);
    v_GPM = Flux(9);
    v_ENO = Flux(10);
    v_PYK = Flux(11);
    v_PDH = Flux(12);
    v_ZWF = Flux(13);
    v_PGL = Flux(14);
    v_GND = Flux(15);
    v_RPE = Flux(16);
    v_RPI = Flux(17);
    v_TKT1 = Flux(18);
    v_TKT2 = Flux(19);
    v_TAL = Flux(20);
    
    v_PRS = Flux(21);

    v_PEPC = Flux(22);

    v_MDH = Flux(23);
    v_MALOut = Flux(24);
    v_GLT = Flux(25);
    v_ACN1 = Flux(26);
    v_ACN2 = Flux(27);
    v_ICL = Flux(28);
    v_MALS = Flux(29);
    v_ICDH = Flux(30);
    v_LPD = Flux(31);
    v_SK = Flux(32);
    v_SUCOut = Flux(33);
    v_SDH = Flux(34);
    v_FUMOut = Flux(35);
    v_FUMA = Flux(36);
    
    v_ToNH4 = Flux(37);
     
    v_GS = Flux(38);
    v_GOGAT = Flux(39);
    v_GD = Flux(40);

    v_PRAT = Flux(41);
    v_PGDH = Flux(42);
    v_PSTS = Flux(43);
    v_SHMT = Flux(44);

    v_GARS = Flux(45);
    v_GARTF = Flux(46);
    v_FGAMS = Flux(47);
   
    v_ALTSS = Flux(48);
    v_AGT = Flux(49);
    v_MTHFR = Flux(50);
    v_MTHFD = Flux(51);
    v_GLYDC = Flux(52);
 
    v_AIRS = Flux(53);

    v_CAIRS = Flux(54);
    v_SS = Flux(55);
    v_ATS = Flux(56);
    v_ADSL = Flux(57);
    v_AICARTF = Flux(58);
    v_IMPCH = Flux(59);
    v_IMPDH = Flux(60);

    v_GLYsink = Flux(61);

    v_NT = Flux(62);
    v_PNP = Flux(63);
    v_XOR = Flux(64);
    v_UOD = Flux(65);
    v_HIUHS = Flux(66);
    v_OHCUD = Flux(67);
    v_ALTNOut = Flux(68);

    %differential equations
    GLCdot = 0;
	G6Pdot = v_HKI - v_PGI - v_ZWF;
	F6Pdot = v_PGI - v_PFK + v_FBP - v_TKT1 - v_TAL;
	FDPdot = v_PFK - v_FBA - v_FBP;
	GAPdot = v_FBA + v_TPI - v_GDH + v_TAL - v_TKT1 - v_TKT2;
	DAPdot = v_FBA - v_TPI;
	BPGdot = v_GDH - v_PGK;
	PGA3dot = v_PGK - v_GPM - v_PGDH;
	PGA2dot = v_GPM - v_ENO;
	PEPdot = v_ENO - v_PYK - v_PEPC;
	PYRdot = v_PYK - v_PDH - v_ALTSS + v_AGT;

    GL6Pdot = v_ZWF - v_PGL;
    PGNdot = v_PGL - v_GND;
    RB5Pdot = v_GND + v_RPE - v_RPI;
    R5Pdot = v_RPI + v_TKT2 - v_PRS;
    X5Pdot = v_TKT1 + v_TKT2 - v_RPE;
    E4Pdot = v_TKT1 - v_TAL;
    S7Pdot = v_TAL - v_TKT2;

    OAAdot = v_PEPC - v_MDH - v_GLT - v_ATS;
    MALdot = v_MDH - v_MALOut + v_FUMA + v_MALS ;

    ACCOAdot = v_PDH - v_GLT - v_MALS;
    CITdot = v_GLT - v_ACN1;
    ACOdot = v_ACN1 - v_ACN2;
    GLXdot = v_ICL - v_MALS - v_AGT;
    ICITdot = v_ACN2 - v_ICDH - v_ICL;
    AKGdot = v_ICDH - v_LPD - v_GOGAT + v_GD + v_PSTS + v_ALTSS + v_ATS;

    SUCCOAdot = v_LPD - v_SK;
    SUCdot = v_SK - v_SUCOut + v_ICL - v_SDH;
    FUMdot = v_SDH - v_FUMOut - v_FUMA + v_ADSL;

    Bdot = v_MALOut + v_GLYsink - ((9.158+0.005)/6.735) * v_ToNH4;

    NH4dot = v_ToNH4 - v_GS + v_GD + v_GLYDC;
    
    GLUdot = - v_GS + 2 * v_GOGAT - v_GD - v_PSTS + v_PRAT - v_ALTSS + v_FGAMS - v_ATS;
    GLNdot = v_GS - v_GOGAT - v_PRAT - v_FGAMS;
     
    PRPPdot = v_PRS - v_PRAT;

    PHP3dot = v_PGDH - v_PSTS;
    SRNdot = v_PSTS - v_SHMT;
    GLYdot = v_SHMT - v_GARS + v_AGT - v_GLYDC - v_GLYsink;

    THFdot = - v_SHMT + v_GARTF - v_GLYDC + v_AICARTF;
    CH2THFdot = v_SHMT - v_MTHFR + v_GLYDC;
    CHTHFdot = v_MTHFR - v_MTHFD;
    CHOTHFdot = v_MTHFD - v_GARTF - v_AICARTF;
    
    PRAdot = v_PRAT - v_GARS;
    GARdot = v_GARS - v_GARTF;
    FGARdot = v_GARTF - v_FGAMS;

    ALNdot = v_ALTSS - v_AGT;

    FGAMdot = v_FGAMS - v_AIRS;
    AIRdot = v_AIRS - v_CAIRS;
    CAIRdot = v_CAIRS - v_SS;
    SAICARdot = v_SS - v_ADSL;
    AICARdot = v_ADSL - v_AICARTF;
    FAICARdot = v_AICARTF - v_IMPCH;
    IMPdot = v_IMPCH - v_IMPDH;
    XMPdot = v_IMPDH - v_NT;
    XAOdot = v_NT - v_PNP;
    XANdot = v_PNP - v_XOR;
    URATEdot = v_XOR - v_UOD;
    HIUHdot = v_UOD - v_HIUHS;
    OHCUdot = v_HIUHS - v_OHCUD;
    ALTNdot = v_OHCUD - v_ALTNOut;

    ASPdot = v_ATS - v_SS;

    %All together
    dxdt = [GLCdot; G6Pdot; F6Pdot; FDPdot; GAPdot; %5
        DAPdot; BPGdot; PGA3dot; PGA2dot; PEPdot; %10
        PYRdot; GL6Pdot; PGNdot; RB5Pdot; R5Pdot; %15
        X5Pdot; E4Pdot; S7Pdot; OAAdot; MALdot; %20
        ACCOAdot; CITdot; ACOdot; ICITdot; AKGdot; %25
        GLXdot; SUCCOAdot; SUCdot; FUMdot; Bdot; %30
        NH4dot; GLUdot; GLNdot; PRPPdot; PHP3dot; %35
        SRNdot;GLYdot;PRAdot; GARdot; FGARdot; %40
        ALNdot; THFdot; CH2THFdot; CHTHFdot;CHOTHFdot; %45
        FGAMdot; AIRdot; CAIRdot; ASPdot; SAICARdot; %50
        AICARdot; FAICARdot; IMPdot; XMPdot; XAOdot;%55
        XANdot; URATEdot; HIUHdot; OHCUdot;ALTNdot %60
        ];        


end


