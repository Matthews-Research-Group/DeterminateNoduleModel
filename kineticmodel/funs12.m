%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script includes equations of different reaction types.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef funs12
    methods(Static)
        %% Inreversible reactions
        % One-substrate ternary reaction
        function z=Function_for_1sub_ternary(A,Km,Vmax)
            z=Vmax*A/(Km+A);end
        
        % Random Bi Bi reaction: A+B->C+D
        function z=Function_for_2sub(A,B,Vmax,KmA,KmB)
            z=Vmax*A*B/((A+KmA)*(B+KmB));end

        % 2-substrate ternary reaction 
        function z=Function_for_2sub_ternary(A,B,Vmax,Kd,KmA,KmB)
            z=Vmax*A*B/(Kd*KmB+KmB*A+KmA*B+A*B);end
        
        % ping-pong reaction 
        function z=Function_for_2sub_pingpong(A,B,Vmax,KmA,KmB)
            z=Vmax*A*B/(KmB*A+KmA*B+A*B);end

        %% Reversible reactions
        % 1-substrate-1-product: A<=>B
        function z=Function_for_Revers1(A,B,Vmax,KmA,KmB,Keq)
            z = Vmax*(A-B/Keq)/(KmA)/(1+A/KmA+B/KmB);end

        % 1-substrate-2-product: A<=>B+C
        function z=Function_for_Revers12(A,B,C,Vmax,KmA,KmB,KmC,Keq)
            z = Vmax*(A-B*C/Keq)/(KmA)/(1+A/KmA+B/KmB+C/KmC+B*C/(KmB*KmC));end
        
        % 2-substrate-1-product: A+B<=>C
        function z=Function_for_Revers21(A,B,C,Vmax,KmA,KmB,KmC,Keq)
            z = Vmax*(A*B-C/Keq)/(KmA*KmB)/((1+A/KmA)*(1+B/KmB)+(1+C/KmC)-1);end

        % 2-substrate-2-product: A+B<=>C+D
        function z=Function_for_Revers22(A,B,C,D,Vmax,KmA,KmB,KmC,KmD,Keq)
            z = Vmax*(A*B-C*D/Keq)/(KmA*KmB)/(1+A/KmA+B/KmB+C/KmC+D/KmD+A*B/(KmA*KmB)+C*D/(KmC*KmD));end

        % pingpong Bi Bi
        function z=PingpongBiBi(A,B,C,D,Keq,Vf,Vr,KmA,KmB,KmC,KmD,KiA,KiD)
            z=Vf*(A*B-C*D/Keq)/(A*B+KmB*A+KmA*B*(1+D/KiD)+Vf/(Vr*Keq)*(KmD*C*(1+A/KiA)+D*(KmC+C)));
        end
       
        %% PP pathway
        function z=Function_for_Den(GAP, F6P, S7P, X5P, E4P, R5P, TKT1_KmGAP, TKT1_KmF6P, TKT2_KmS7P, TKT2_KmGAP, TKT2_KmP5P, TKT1_KmP5P)
            z = 1+(1+GAP/TKT1_KmGAP)*(F6P/TKT1_KmF6P+S7P/TKT2_KmS7P)+...
                GAP/TKT2_KmGAP + 1/TKT2_KmP5P * (X5P*(1+E4P*R5P/TKT1_KmP5P)+E4P+R5P);end

        % TKT1:GAP+F6P <=> E4P+X5P
        function z=Function_for_TKT1(F6P,GAP,E4P, X5P, Vmax, Keq, TKT2_KmP5P, TKT1_KmP5P, Den)
            z=Vmax * (F6P*GAP*Keq-E4P*X5P)/(TKT2_KmP5P*TKT1_KmP5P*Den);end
        
        % TKT2:GAP+S7P <=> R5P+X5P
        function z=Function_for_TKT2(S7P, GAP, R5P, X5P, Vmax, Keq, TKT2_KmP5P, TKT1_KmP5P, Den) 
            z=Vmax * (S7P*GAP*Keq-X5P*R5P)/(TKT2_KmP5P*TKT1_KmP5P*Den);end
        
               
        function z= Function_Permeability(P, S1, S2, SA, D)
            z = P*(S1-S2)*SA/(4/3*3.14*((D*10^(-6))/2)^3);end


    end
end


