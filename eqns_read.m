function [dydt] = eqns_read(t,y,sp,pa)
%eqns_read Differential equations to model binding
%   All receptors (R1, R2) can bind all ligands (La,Lb,Lc) 
%   Once bound, they form complexes (e.g., R1_La, R1_Lb, etc.)
%   INPUTS:
%      t = time (passed in automatically by ode-solver) 
%      y = value of species 
%      pa = struct storing parameter values 
%      sp = struct storing species indices (e.g., R1 = 1, means R1 is the first index of y) 
% 
%% CODE SNIPPETS
bind_R1_La = pa.kon_R1_La*y(sp.R1)*y(sp.La);
unbind_R1_La = pa.koff_R1_La*y(sp.R1_La);

bind_R1_Lb = pa.kon_R1_Lb*y(sp.R1)*y(sp.Lb);
unbind_R1_Lb = pa.koff_R1_Lb*y(sp.R1_Lb);

bind_R1_Lc = pa.kon_R1_Lc*y(sp.R1)*y(sp.Lc);
unbind_R1_Lc = pa.koff_R1_Lc*y(sp.R1_Lc);

bind_R2_La = pa.kon_R2_La*y(sp.R2)*y(sp.La);
unbind_R2_La = pa.koff_R2_La*y(sp.R2_La);

bind_R2_Lb = pa.kon_R2_Lb*y(sp.R2)*y(sp.Lb);
unbind_R2_Lb = pa.koff_R2_Lb*y(sp.R2_Lb);

bind_R2_Lc = pa.kon_R2_Lc*y(sp.R2)*y(sp.Lc);
unbind_R2_Lc = pa.koff_R2_Lc*y(sp.R2_Lc);


%% Equations
dydt = zeros(length(y),1);
dydt(sp.R1) = - bind_R1_La + unbind_R1_La ...
              - bind_R1_Lb + unbind_R1_Lb ...
              - bind_R1_Lc + unbind_R1_Lc;
dydt(sp.R2) = - bind_R2_La + unbind_R2_La ...
              - bind_R2_Lb + unbind_R2_Lb ...
              - bind_R2_Lc + unbind_R2_Lc;
dydt(sp.La) = - bind_R1_La + unbind_R1_La ...
              - bind_R2_La + unbind_R2_La;
dydt(sp.Lb) = - bind_R1_Lb + unbind_R1_Lb ...
              - bind_R2_Lb + unbind_R2_Lb;
dydt(sp.Lc) = - bind_R1_Lc + unbind_R1_Lc ...
              - bind_R2_Lc + unbind_R2_Lc;
dydt(sp.R1_La) = + bind_R1_La - unbind_R1_La;
dydt(sp.R1_Lb) = + bind_R1_Lb - unbind_R1_Lb;
dydt(sp.R1_Lc) = + bind_R1_Lc - unbind_R1_Lc;
dydt(sp.R2_La) = + bind_R2_La - unbind_R2_La;
dydt(sp.R2_Lb) = + bind_R2_Lb - unbind_R2_Lb;
dydt(sp.R2_Lc) = + bind_R2_Lc - unbind_R2_Lc;
end

