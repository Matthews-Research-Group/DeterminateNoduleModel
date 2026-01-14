%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script includes the rate laws for all reactions included 
% in the model.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef rxns
    methods(Static)
        %% Glycolysis
        Vmax=Vmaxvalue();
        params=kineticparams(Vmax);
        
        function r = reaction_HKI(params, GLC, G6P)
            r=funs12.Function_for_Revers1(GLC,G6P,params.HKI_Vmax,params.HKI_KmGLC,params.HKI_KmG6P,params.HKI_Keq);end

        function r = reaction_PGI(params, G6P, F6P)
            r=funs12.Function_for_Revers1(G6P,F6P,params.PGI_Vmax,params.PGI_KmG6P,params.PGI_KmF6P,params.PGI_Keq);end

        function r = reaction_PFK(params, F6P)
            r=funs12.Function_for_1sub_ternary(F6P,params.PFK_KmF6P,params.PFK_Vmax);end
                    
        function r = reaction_FBP(params, FDP)
            r=funs12.Function_for_1sub_ternary(FDP,params.FBP_KmFDP,params.FBP_Vmax);end
        
        function r = reaction_FBA(params,FDP,GAP,DAP) 
            r=funs12.Function_for_Revers12(FDP,GAP,DAP,params.FBA_Vmax,params.FBA_KmFDP,params.FBA_KmGAP,params.FBA_KmDAP,params.FBA_Keq);end
        
       
        function r = reaction_TPI(params,DAP,GAP,PGA3,PEP)
            r=(params.TPI_Vmax*DAP)*(DAP-GAP/params.TPI_Keq)/(params.TPI_KmDAP)/(1+GAP/params.TPI_KmGAP+DAP/params.TPI_KmDAP+PGA3/params.TPI_Ki3PG+PEP/params.TPI_KiPEP);
        end

        function r = reaction_GDH(params,GAP,BPG)
            r=funs12.Function_for_Revers1(GAP,BPG,params.GDH_Vmax,params.GDH_KmGAP,params.GDH_KmBPG,params.GDH_Keq);end

        function r = reaction_PGK(params,BPG,PGA3)
            r=funs12.Function_for_Revers1(BPG,PGA3,params.PGK_Vmax,params.PGK_KmBPG,params.PGK_KmPGA3,params.PGK_Keq);end
        
        function r = reaction_GPM(params,PGA3,PGA2)
            r=funs12.Function_for_Revers1(PGA3,PGA2,params.GPM_Vmax,params.GPM_KmPGA3,params.GPM_KmPGA2,params.GPM_Keq);end
        
        function r = reaction_ENO(params,PGA2,PEP)
            r=funs12.Function_for_Revers1(PGA2,PEP,params.ENO_Vmax,params.ENO_KmPGA2,params.ENO_KmPEP,params.ENO_Keq);end
        
        function r = reaction_PYK(params, PEP)
            r=funs12.Function_for_1sub_ternary(PEP,params.PYK_KmPEP,params.PYK_Vmax);end

        function r = reaction_PDH(params, PYR)
            PDH_KmPYR_new = params.PDH_KmPYR*(1+PYR/params.PDH_KiPYR);
            r=funs12.Function_for_1sub_ternary(PYR,PDH_KmPYR_new, params.PDH_Vmax);end
  

        %% PP Pathway
        function r = reaction_ZWF(params, G6P,GL6P)
            r=(params.ZWF_Vmax)*(G6P-GL6P/params.ZWF_Keq)/(params.ZWF_KdG6P)/(1+G6P/params.ZWF_KdG6P+GL6P/params.ZWF_KdGL6P);
        end

        function r=PGL(params,GL6P, PGN, G6P)
            r=(params.PGL_Vmax)*(GL6P-PGN/params.PGL_Keq)/(params.PGL_KmGL6P)/(1+GL6P/params.PGL_KmGL6P+PGN/params.PGL_KmPGN+G6P/params.PGL_KiG6P);
        end

        function r = reaction_GND(params, PGN, RB5P)
            r=params.GND_Vmax*(PGN-RB5P/params.GND_Keq)/(params.GND_KmPGN)/(1+PGN/params.GND_KmPGN+RB5P/params.GND_KdRB5P);
        end

        function r = reaction_RPE(params,X5P,RB5P)
            r=funs12.Function_for_Revers1(X5P,RB5P,params.RPE_Vmax,params.RPE_KmX5P,params.RPE_KmRB5P, params.RPE_Keq);end
        
        function r = reaction_RPI(params,RB5P,R5P,E4P,GAP,PGA3,PGN)
            r=(params.RPI_Vmax/(1+PGN/params.RPI_KiR5P_InhPGN))*(RB5P-R5P/params.RPI_Keq)/(params.RPI_KmRB5P*(1+E4P/params.RPI_KiRB5P_inhE4P))/(1+RB5P/(params.RPI_KmRB5P*(1+E4P/params.RPI_KiRB5P_inhE4P))+R5P/(params.RPI_KmR5P*(1+E4P/params.RPI_KiR5P_inhE4P+GAP/params.RPI_KiR5P_inhGAP+PGA3/params.RPI_KiR5P_inhPGA3)));
        end

        function r = TAL(params, F6P, E4P, S7P, GAP)
            r = funs12.Function_for_Revers22(F6P, E4P, S7P, GAP,...
                params.TAL_Vmax, params.TAL_KmF6P, params.TAL_KmE4P,...
                params.TAL_KmS7P, params.TAL_KmGAP, params.TAL_Keq);end

        function r = TKT1(params,GAP, F6P, S7P, X5P, E4P, R5P)
            Den = funs12.Function_for_Den(GAP, F6P, S7P, X5P, E4P, R5P,...
                params.TKT1_KmGAP, params.TKT1_KmF6P, params.TKT2_KmS7P,... 
                params.TKT2_KmGAP, params.TKT2_KmP5P, params.TKT1_KmP5P);

            r = funs12.Function_for_TKT1(F6P,GAP,E4P, X5P, params.TKT1_Vmax,...
                params.TKT1_Keq, params.TKT2_KmP5P,...
                params.TKT1_KmP5P, Den);end

        function r = TKT2(params,GAP, F6P, S7P, X5P, E4P, R5P)
            Den = funs12.Function_for_Den(GAP, F6P, S7P, X5P, E4P, R5P,...
                params.TKT1_KmGAP, params.TKT1_KmF6P, params.TKT2_KmS7P,... 
                params.TKT2_KmGAP, params.TKT2_KmP5P, params.TKT1_KmP5P);

            r = funs12.Function_for_TKT2(S7P, GAP, R5P, X5P, params.TKT2_Vmax,...
                params.TKT2_Keq, params.TKT2_KmP5P,...
                params.TKT1_KmP5P, Den); end

        function r=PGM(params,R1P,R5P)
            r=funs12.Function_for_Revers1(R1P,R5P,params.PGM_Vmax,params.PGM_KmR1P,params.PGM_KmR5P,params.PGM_Keq);end

        %% Anaplerotic reactions
        function r = reaction_PCK(params, OAA)
            r=funs12.Function_for_1sub_ternary(OAA,params.PCK_KmOAA,params.PCK_Vmax);end
    
        function r=PEPC(params,PEP,MAL)
            PEPC_KmPEP_new = params.PEPC_KmPEP*(1+(MAL/params.PEPC_KiMAL));
            r=params.PEPC_Vmax*PEP/(PEP+PEPC_KmPEP_new);
        end
        
        %% TCA cycle
        
        function r = reaction_GLT(params, ACCOA, OAA,CIT,AKG)
            GLT_KmOAA_new = params.GLT_KmOAA*(1+OAA/params.GLT_KiOAA_inhOAA+AKG/params.GLT_KiOAA_inh2kg);
            r=funs12.Function_for_Revers21(ACCOA,OAA,CIT,params.GLT_Vmax,params.GLT_KmACCOA, ...
                GLT_KmOAA_new, params.GLT_KmCIT,params.GLT_Keq);
        end
               
        function r = reaction_ACN1(params, CIT,ACO)
            r=funs12.Function_for_Revers1(CIT,ACO,params.ACN1_Vmax,params.ACN1_KmCIT,params.ACN1_KmACO,params.ACN1_Keq);end
       
        function r = reaction_ACN2(params, ACO,ICIT)
            r=funs12.Function_for_Revers1(ACO,ICIT,params.ACN2_Vmax,params.ACN2_KmACO,params.ACN2_KmICIT,params.ACN2_Keq);end

        function r = reaction_ICD(params, ICIT, AKG)
            ICD_KmAKG_new=params.ICD_KmAKG*(1+ICIT/params.ICD_Ki2kg_inhICIT);
            r=funs12.Function_for_Revers1(ICIT, AKG, params.ICD_Vmax,params.ICD_KmICIT,ICD_KmAKG_new,params.ICD_Keq);
        end

        function r = reaction_LPD(params, AKG)
            r=funs12.Function_for_1sub_ternary(AKG,params.LPD_KmAKG,params.LPD_Vmax);end


        function r = reaction_SK(params, SUCCOA)
            r=funs12.Function_for_1sub_ternary(SUCCOA,params.SK_KmSUCCOA,params.SK_Vmax);
        end
        
        function r = reaction_SDH(params, SUC, FUM, OAA)
            SDH_KmSUC_new = params.SDH_KmSUC*(1+OAA/params.SDH_KiSUC_inhOAA);
            r=funs12.Function_for_Revers1(SUC,FUM,params.SDH_Vmax,SDH_KmSUC_new,params.SDH_KmFUM,params.SDH_Keq);
        end

        function r = reaction_FUMA(params, FUM,MAL,PEP,PYR,CIT,AKG,OAA)
            FUMA_KmFUM_new = params.FUMA_KmFUM*(1+PEP/params.KiFUM_inhPEP+PYR/params.KiFUM_inhPYR+CIT/params.KiFUM_inhCIT+AKG/params.KiFUM_inhAkg+OAA/params.KiFUM_inhOAA);
            r=funs12.Function_for_Revers1(FUM,MAL,params.FUMA_Vmax,FUMA_KmFUM_new,params.FUMA_KmMAL,params.FUMA_Keq);
        end
 
        function r = reaction_MDH(params,OAA,MAL)
             r=funs12.Function_for_Revers1(OAA,MAL,params.MDH_Vmax,params.MDH_KmOAA,params.MDH_KmMAL,params.MDH_Keq);end


        %% Glyoxylate shunt
        
        function r = reaction_ICL(params, ICIT)
            r=funs12.Function_for_1sub_ternary(ICIT,params.ICL_KmICIT,params.ICL_Vmax);end

        function r = reaction_MALS(params,ACCOA,GLX,MAL,PYR,OAA)
            MALS_KmGLX_new = params.MALS_KmGLX*(1+PYR/params.MALS_KiGLX_inhPYR+OAA/params.MALS_KiGLX_inhOAA);
            r=funs12.Function_for_Revers21(ACCOA,GLX,MAL,params.MALS_Vmax,params.MALS_KmACCOA,MALS_KmGLX_new,params.MALS_KmMAL,params.MALS_Keq);
        end
  
        %% Acetate metabolism
        
        function r = reaction_PDC(params, PYR)
            r=funs12.Function_for_1sub_ternary(PYR,params.PDC_KmPYR,params.PDC_Vmax);end
        
        function r = reaction_ALDH(params,ACLD,ACE)
            r=funs12.Function_for_Revers1(ACLD,ACE,params.ALDH_Vmax,params.ALDH_KmACLD,params.ALDH_KmACE,params.ALDH_Keq);end

        function r = reaction_ACS(params, ACE)
            r=funs12.Function_for_1sub_ternary(ACE,params.ACS_KmACE,params.ACS_Vmax);end
        
        %% GABA Shunt
        function r = reaction_SSADH(params, SSA)
            r=funs12.Function_for_1sub_ternary(SSA,params.SSADH_KmSSA,params.SSADH_Vmax);end

        function r = reaction_GABAPT1(params, GABA,GLX,SSA,GLY)
            r=funs12.Function_for_Revers22(GABA,GLX,SSA,GLY,params.GABAPT1_Vmax, ...
                params.GABAPT1_KmGABA,params.GABAPT1_KmGLX,params.GABAPT1_KmSSA, ...
                params.GABAPT1_KmGLY, params.GABAPT1_Keq);end

        function r = reaction_GABAPT2(params, GABA,PYR,SSA,ALN,GLY)
            GABAPT2_KmGABA_new = params.GABAPT2_KmGABA*(1+GLY/params.GABAPT2_KiGABA_inhGLY);
            r=funs12.Function_for_Revers22(GABA,PYR,SSA,ALN,params.GABAPT2_Vmax, ...
                GABAPT2_KmGABA_new,params.GABAPT2_KmPYR,params.GABAPT2_KmSSA, ...
                params.GABAPT2_KmALN, params.GABAPT2_Keq);end


        function r = reaction_GDC(params, GLU)
            r=funs12.Function_for_1sub_ternary(GLU,params.GDC_KmGLU,params.GDC_Vmax);end
               
        %% GOGAT
        
        function r = reaction_GD(params,GLU,AKG,NH4)
            GD_KmAKG_new = params.GD_KmAKG*(1+GLU/params.GD_Ki_2kg_inhGLU);
            r=funs12.Function_for_Revers12(GLU,AKG,NH4,params.GD_Vmax,params.GD_KmGLU,GD_KmAKG_new,params.GD_KmNH4,params.GD_Keq);
        end

        function r = reaction_GOGAT(params, GLN, AKG, GLU, OAA, ASP)
            r=params.GOGAT_Vmax*GLN*AKG/(params.GOGAT_KmAKG*(1+(GLU/params.GOGAT_KiAKG_inhGLU)+(OAA/params.GOGAT_KiAKG_inhOAA)+(ASP/params.GOGAT_KiAKG_inhASP))*GLN+params.GOGAT_KmGLN*AKG+GLN*AKG);
        end

        function r = reaction_GS(params, GLU, NH4,GLN,ALN,GLY)
            GS_Vmaxnew = params.GS_Vmax/(1+GLN/params.GS_Ki_GLU_inhGLN+GLN/params.GS_Ki_NH4_inhGLN);
            GS_KmGLNnew = params.GS_KmGLN*(1+ALN/params.GS_Ki_GLN_InhALN+GLY/params.GS_Ki_GLN_InhGLY);
            r=funs12.Function_for_Revers21(GLU,NH4,GLN,GS_Vmaxnew,params.GS_KmGLU,params.GS_KmNH4,GS_KmGLNnew,params.GS_Keq);
        end
       
        function r=GG(params,NH4,AKG)
            r=funs12.Function_for_2sub_pingpong(NH4,AKG,params.GG_Vmax,params.GG_KmNH4,params.GG_KmAKG);
        end

        %% De novo purine synthesis
        
        function r = reaction_PRS(params, R5P,PRPP)
            PRS_KmR5P_new = params.PRS_KmR5P*(1+PRPP/params.PRS_KisR5P_inhPRPP);
            PRS_R5P_new = R5P*(1+PRPP/params.PRS_KiiR5P_inhPRPP);
            r=funs12.Function_for_Revers1(PRS_R5P_new,PRPP,params.PRS_Vmax,PRS_KmR5P_new,params.PRS_KmPRPP, params.PRS_Keq);
        end

        function r = reaction_PRAT(params, PRPP, GLN, PRA, GLU, NH4,IMP)
            PRAT_KmPRPP_new = params.PRAT_KmPRPP*(1+IMP/params.PRAT_KiPRPP_inhIMP);%+XMP/params.PRAT_KiPRPP_inhXMP
            PRAT_KmGLN_new = params.PRAT_KmGLN*(1+NH4/params.PRAT_KiGLN_inhNH4+GLU/params.PRAT_KiGLN_inhGLU);
            r=funs12.Function_for_Revers22(PRPP,GLN,PRA,GLU,params.PRAT_Vmax,PRAT_KmPRPP_new,PRAT_KmGLN_new,params.PRAT_KmPRA, params.PRAT_KmGLU, params.PRAT_Keq);
        end

        function r = reaction_GARS(params, PRA, GLY, GAR)
            r=funs12.Function_for_Revers21(PRA, GLY, GAR,params.GARS_Vmax,params.GARS_KmPRA,params.GARS_KmGLY,params.GARS_KmGAR, params.GARS_Keq);
        end

        function r = reaction_GARTF(params, GAR,CHOTHF,THF)
            GARTF_KmCHOTHFnew = params.GARTF_KmCHOTHF*(1+THF/params.GARTF_CHOTHF_inhTHF);
            r=funs12.Function_for_2sub_pingpong(GAR,CHOTHF,params.GARTF_Vmax,params.GARTF_KmGAR,GARTF_KmCHOTHFnew);
        end
         
        function r = reaction_FGAMS(params, FGAR, GLN, GLU)
            FGAMS_KmFGAR_new = params.FGAMS_KmFGAR * (1+GLU/params.FGAMS_KisFGAR_inhGLU);
            FGAR_new = FGAR*(1+GLU/params.FGAMS_KiiFGAR_inhGLU);
            FGAMS_KmGLN_new = params.FGAMS_KmGLN*(1+GLU/params.FGAMS_KiGLN_inhGLU);
            r=funs12.Function_for_2sub_pingpong(FGAR_new,GLN,params.FGAMS_Vmax,FGAMS_KmFGAR_new,FGAMS_KmGLN_new);
        end

        function r = reaction_AIRS(params, FGAM)
            r=funs12.Function_for_1sub_ternary(FGAM,params.AIRS_KmFGAM,params.AIRS_Vmax);
        end       

        function r = reaction_CAIRS(params, AIR)
            r=funs12.Function_for_1sub_ternary(AIR,params.CAIRS_KmAIR,params.CAIRS_Vmax);
        end 

        function r = reaction_SS(params, CAIR,ASP)
            r=funs12.Function_for_2sub_pingpong(CAIR,ASP,params.SS_Vmax,params.SS_KmCAIR,params.SS_KmASP);end

        function r = reaction_ADSL(params, SAICAR,FUM,AICAR)
            r=funs12.Function_for_Revers12(SAICAR,FUM,AICAR,params.ADSL_Vmax,params.ADSL_KmSAICAR,params.ADSL_KmFUM,params.ADSL_KmAICAR,params.ADSL_Keq);end    
       

        function r = reaction_AICARTF(params,AICAR,CHOTHF)
            r=funs12.Function_for_2sub_pingpong(AICAR,CHOTHF,params.AICARTF_Vmax,params.AICARTF_KmAICAR,params.AICARTF_KmCHOTHF);
        end    

        function r = reaction_IMPCH(params, FAICAR,IMP)
            r=funs12.Function_for_Revers1(FAICAR,IMP,params.IMPCH_Vmax,params.IMPCH_KmFAICAR,params.IMPCH_KmIMP,params.IMPCH_Keq);
        end

        %% Purine to xanthine:

        % Reaction:id=HPTF
        % xanthine + PRPP <=> IMP+diphosphate
        function r=HPTF(params,XAN,PRPP,IMP)
            HPTF_KmIMPnew = params.HPTF_KmIMP*(1+PRPP/params.HPTF_KiIMP_inhPRPP);
            HPTF_KmXANnew = params.HPTF_KmXAN*(1+IMP/params.HPTF_KiXAN_inhIMP);
            r=funs12.Function_for_Revers21(XAN,PRPP,IMP,params.HPTF_Vmax, HPTF_KmXANnew, params.HPTF_KmPRPP, HPTF_KmIMPnew, params.HPTF_Keq);
        end

        function r = reaction_IMPDH(params, IMP)
            r=funs12.Function_for_1sub_ternary(IMP,params.IMPDH_KmIMP,params.IMPDH_Vmax);
        end  

        function r = reaction_NT(params, XMP)
            r=funs12.Function_for_1sub_ternary(XMP,params.NT_KmXMP,params.NT_Vmax);
        end  

        function r = reaction_PNP(params, XAO)
            r=funs12.Function_for_1sub_ternary(XAO,params.PNP_KmXao,params.PNP_Vmax);
        end

        function r = reaction_XOR(params, XAN, URATE)
            XOR_KmXAN_new = params.XOR_KmXan*(1+URATE/params.XOR_KiXAN_inhUrate);
            r=funs12.Function_for_1sub_ternary(XAN,XOR_KmXAN_new,params.XOR_Vmax);
        end  

        function r = reaction_UOD(params, URATE, XAN)
            UOD_KmURATE_new = params.UOD_KmURATE*(1+XAN/params.UOD_Ki_Urate_inhXAN);
            r=funs12.Function_for_1sub_ternary(URATE,UOD_KmURATE_new, params.UOD_Vmax);
        end
        
        function r = reaction_HIUHS(params, HIUH)
            r=funs12.Function_for_1sub_ternary(HIUH,params.HIUHS_KmHIUH,params.HIUHS_Vmax);
        end


        function r = reaction_OHCUD(params, OHCU)
            r=funs12.Function_for_1sub_ternary(OHCU,params.OHCUD_KmOHCU, params.OHCUD_Vmax);
        end


        %% Aspartate metabolism
        function r = reaction_ATS(params,OAA,GLU,ASP,AKG)
            ATS_KmASPnew = params.ATS_KmASP*(1+GLU/params.ATS_KiASP_inhGLU);
            ATS_KmAKGnew = params.ATS_KmAKG*(1+OAA/params.ATS_KiAKG_inhOAA);
            ATS_KmOAAnew = params.ATS_KmOAA*(1+OAA/params.ATS_KiOAA_inhOAA);
            r=funs12.Function_for_Revers22(OAA,GLU,ASP,AKG,params.ATS_Vmax,ATS_KmOAAnew,params.ATS_KmGLU,ATS_KmASPnew,ATS_KmAKGnew,params.ATS_Keq);
        end

        function r=reaction_ASNS(params,ASP,GLN)
            r=funs12.Function_for_2sub_pingpong(ASP,GLN,params.ASNS_Vmax, params.ASNS_KmASP, params.ASNS_KmGLN);
        end

        function r = reaction_ASNA(params,ASN)
            r=funs12.Function_for_1sub_ternary(ASN,params.ASNA_KmASN,params.ASNA_Vmax);
        end

        function r = ASPA(params,FUM,ASP,NH4)
            r=funs12.Function_for_Revers12(FUM,ASP,NH4,params.ASPA_Vmax,params.ASPA_KmFUM,params.ASPA_KmASP,params.ASPA_KmNH4,params.ASPA_keq);
        end

        %% Alanine metabolism
        function r = reaction_ALTSS(params,PYR,GLU,ALN,AKG)
            r=funs12.Function_for_Revers22(PYR,GLU,ALN,AKG,params.ALTSS_Vmax,params.ALTSS_KmPYR,params.ALTSS_KmGLU,params.ALTSS_KmALN,params.ALTSS_KmAKG,params.ALTSS_Keq);
        end

        %% Serine and glycine metabolism
        function r = reaction_PGDH(params,PGA3, PHP3)
            r=funs12.Function_for_Revers1(PGA3,PHP3,params.PGDH_Vmax,params.PGDH_Km3PG,params.PGDH_Km3PHP,params.PGDH_Keq);
        end

        function r = reaction_newPSTS(params,PHP3,GLU,SRN,AKG)
            r=funs12.Function_for_Revers22(PHP3,GLU,SRN,AKG,params.PSTS_Vmax,params.PSTS_Km3PHP,params.PSTS_KmGLU,params.PSTS_KmSerine,params.PSTS_KmAKG,params.PSTS_Keq);
        end

        function r = reaction_PSPH(params,PPSR)
            r=funs12.Function_for_1sub_ternary(PPSR,params.PSPH_KmPPSR,params.PSPH_Vmax);
        end

        function r = reaction_PATF(params, SRN,GLX,HDPYR,GLY)
            r=funs12.Function_for_Revers22(SRN,GLX,HDPYR,GLY,params.PATF_Vmax,params.PATF_KmSerine,params.PATF_KmGLX,params.PATF_KmHDPYR,params.PATF_KmGLY,params.PATF_Keq);
        end

        function r = reaction_SHMT(params,SRN,THF,GLY,CH2THF)
            SHMT_KmSerine_new = params.SHMT_KmSerine*(1+GLY/params.SHMT_KiSRN_inhGLY);
            SHMT_Vmax_new = params.SHMT_Vmax/(1+CH2THF/params.SHMT_KiTHF_inhCH2THF);
            r=funs12.Function_for_Revers22(SRN,THF,GLY,CH2THF,SHMT_Vmax_new,SHMT_KmSerine_new,params.SHMT_KmTHF,params.SHMT_KmGLY,params.SHMT_KmCH2THF,params.SHMT_Keq);
        end

        function r=reaction_MTHFR(params,CH2THF)
            MTHFR_KmCH2THF_new = params.MTHFR_KmCH2THF*(1+CH2THF/params.MTHFR_KiCH2THF);
            r=funs12.Function_for_1sub_ternary(CH2THF,MTHFR_KmCH2THF_new,params.MTHFR_Vmax);
        end

        function r=reaction_MTHFD(params,CHTHF)
            r=funs12.Function_for_1sub_ternary(CHTHF,params.MTHFD_KmCHTHF,params.MTHFD_Vmax);
        end        

        function r=reaction_GLYDC(params,GLY,THF)
            r=funs12.Function_for_2sub_pingpong(GLY,THF,params.GLYDC_Vmax, params.GLYDC_KmGLY, params.GLYDC_KmTHF);
        end

        function r = reaction_GLYDH(params, HDPYR, GLYCT)
            r=funs12.Function_for_Revers1(HDPYR, GLYCT,params.GLYDH_Vmax,params.GLYDH_KmHDPYR,params.GLYDH_KmGLYCT,params.GLYDH_Keq);
        end

        function r = reaction_PGAP(params, GLYCT)
            r=funs12.Function_for_1sub_ternary(GLYCT,params.PGAP_KmGLYCT,params.PGAP_Vmax);
        end

        function r = reaction_AMTF(params,GLY,THF,CHOTHF10,NH4)
            r=funs12.Function_for_Revers22(GLY,THF,CHOTHF10,NH4,params.AMTF_Vmax,params.AMTF_KmGLY,params.AMTF_KmTHF,params.AMTF_Km10CHOTHF,params.AMTF_KmNH4,params.AMTF_Keq);
        end

        function r = reaction_AGT(params,ALN,GLX,PYR,GLY)
            AGT_KmALN_new = params.AGT_KmALN*(1+PYR/params.AGT_KiALN_inhPYR_Kii);
            ALN_new = ALN*(1+PYR/params.AGT_KiALN_inhPYR_Kis);
            AGT_KmGLX_new = params.AGT_KmGLX*(1+PYR/params.AGT_KiGLX_inhPYR);
            r=funs12.Function_for_Revers22(ALN_new,GLX,PYR,GLY,params.AGT_Vmax,AGT_KmALN_new,AGT_KmGLX_new,params.AGT_KmPYR,params.AGT_KmGLY,params.AGT_Keq);
        end


        %% Metabolite transfer
        function r = ALTNOut(params,ALTN)
            r=funs12.Function_for_1sub_ternary(ALTN,params.ALTNOut_Km,params.ALTNOut_Vmax);
        end

        function r = MALOut(params,MAL)
            r=funs12.Function_for_1sub_ternary(MAL,params.SourceMAL_KmMAL,params.SourceMAL_Vmax);
        end

        function r = SUCOut(params,SUC)
            r=funs12.Function_for_1sub_ternary(SUC,params.SourceSUC_KmSUC,params.SourceSUC_Vmax);
        end
     
        function r = FUMOut(params,FUM)
            r=funs12.Function_for_1sub_ternary(FUM,params.SourceFUM_KmFUM,params.SourceFUM_Vmax);
        end



    end
end

