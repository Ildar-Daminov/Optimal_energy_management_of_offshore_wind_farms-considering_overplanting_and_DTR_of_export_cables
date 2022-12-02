function I_adm=I_admissible(Wd,T1,T2,T3,T4,n,R,Lambda1,Lambda2,Tmax,Tamb)

% This function finds the admissible current in accordance with 
% IEC 60287-1-1, p11
I_adm=sqrt(((Tmax-Tamb)-Wd*(0.5*T1+n*(T2+T3+T4)))/(R*(T1+n*(1+Lambda1)*T2+n*(1+Lambda1+Lambda2)*(T3+T4))));

% Input:
% Wd - dielectric losses 
% T1 - thermal resistance in RC scheme 
% T2 - thermal resistance in RC scheme 
% T3 - thermal resistance in RC scheme 
% T4 - thermal resistance in RC scheme 
% n  - number of cores (3)
% R  - electrical resistance 
% Lambda1 - Ratio of the total losses in sheath respectively to the total
%           conductor losses (or losses in one sheath to the losses in one
%           conductor)(p24 IEC 60287-1-1)
% Lambda2,  Ratio of the total losses in armour respectively to the total
%           conductor  losses (or losses in one sheath or armour to the losses
%           in one conductor) p28 CEI 60287-1-1
% Tmax - maximal allowable temperature of cable, degC
% Tamb - ambient temperature, degC

% Output:
% I_adm - admissible current of cable, in A