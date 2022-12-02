function [Theta_ct,Theta_Et,Alpha]=losses_partial_transients_new(I,R_AC,Lambda1_new)
%% Function input/output
% This function is similar to the losses_partial_transients. m but it requires 
% 2 additional inputs (R_AC,Lambda1_new). These two parameters are changing
% as a function of temperature. Hence, this function takes into account
% their impact on Theta_ct,Theta_Et and Alpha

% This function is used only inside of temperature_correction_on_losses_new.m 

% INPUT
% I - current step, A
% R_AC - Elctrical resistance at given temperature
% Lambda1_new  Ratio of the total losses in sheath respectively to the total
% conductor losses (or losses in one sheath to the losses in one
% conductor)(p24 IEC 60287-1-1)

% OUTPUT
% Theta_ct - cable partial transient,K
% Theta_Et - cable environment partial transient, K
% Alpha    - conductor to cable surface attainment factor, pu
global  Lambda2  t De Diff_sol Depth
global Rho_sol Wd Ta Tb a b TA TB

%% Losses estimation
% p27 IEC60287-1-1
% 2.4.2.1 Single-core lead-sheathed cables – steel wire armour, bonded to
% sheath at both ends
% Wc is the power loss per unit length (W/m)in a conductor or an equivalent
% conductor based on the maximum conductor temperature attained.
% The power loss is assumed to be constant during the transient.
Wc=3*I^2*R_AC;                                                       % or Wc=I^2*R_AC90_eval (without multiplication by 3);
% Wc=I^2*R_AC90_eval;

% total I^2*R (ohmic) power loss of each cable
WI=Wc*(1+Lambda1_new+Lambda2);  % excluding Wd
% WI=Wc*(1+Lambda1+Lambda2)+3*Wd;  % including Wd

%  But (p12  Nielsen, Applied Sciences 2019)says that WI includes Wd...
% Total internal cable losses WI, including dielectric, sheath and armour
% losses 
%% 4.2.3 Calculation of cable partial transient (p 33 IEC 60853-2)
Theta_ct=Wc.*(Ta.*(1-exp(-a.*t))+Tb.*(1-exp(-b.*t)));
% Theta_ct=3*Wc.*(Ta.*(1-exp(-a.*t))+Tb.*(1-exp(-b.*t)));

% Where t is a time from start of application of heating, a general symbol
% for time, usually in seconds (t = 3 600i)

% The conductor to cable surface attainment factor Alpha is then obtained
% Alpha=Theta_ct/(Wc*(TA+TB));   %  iniial formula with Wc but see next

% However, if I==0 then Wc==0 therefore Alpha would turn to NaN since Wc
% is the denominator of alpha (i.e. dividing by 0). Hence,we suggest to 
% avoid using Wc in alpha as follows:
Alpha=(Ta.*(1-exp(-a.*t))+Tb.*(1-exp(-b.*t)))./(TA+TB);                                              

%% 4.2.4 Calculation of cable environment partial transient (p33 IEC853-2)

% argument of the first integral exponential
X=(De^2)./(16.*Diff_sol.*t);

% argument of the second integral exponential
Y=(Depth^2)./(Diff_sol.*t);

% NOTES TO FIGURES 6 AND 7 "SCALES FOR EXPONENTIAL INTEGRAL"(p101 IEC853-2)
% preallocation of X_1 and Y_1 to increase the speed 
% X_1=NaN(length(X),1);
% Y_1=NaN(length(Y),1);
% for i=1:length(X)
%     if (X(i)>=0 && X(i)<=1)
%         X_1(i)=-0.5772-log(X(i))+X(i)-0.2499*X(i)^2+0.0552*X(i)^3+...
%             0.0098*X(i)^4+0.0011*X(i)^5;                                    % contradiction between IEC 853-2 (+0.0098*X^4) p101 and Dorison (-0.0098*X^4) p142
%     elseif (X(i)>1 && X(i)<8)                                               % contradiction between IEC 853-2 (x>8, E=formula) p101 and Dorison (x>=8, E=formula) p142
%         X_1(i)=1/(X(i)*exp(X(i)))*((X(i)^2+2.3347*X(i)+0.2506)/(X(i)^2+...
%             3.3307*X(i)+1.6815));
%     elseif (X(i)<0.01)
%         X_1(i)=-log(X(i))-0.5772;
%     else
%         X_1(i)=0;
%     end
%     
%     if (Y(i)>=0 && Y(i)<=1)
%         Y_1(i)=-0.5772-log(Y(i))+Y(i)-0.2499*Y(i)^2+0.0552*Y(i)^3+...
%             0.0098*Y(i)^4+0.0011*Y(i)^5;                                        % contradiction between IEC 853-2 (+0.0098*X^4) p101 and Dorison (-0.0098*X^4) p142
%     elseif (Y(i)>1 && Y(i)<8)                                                   % contradiction between IEC 853-2 (x>8, E=formula) p101 and Dorison (x>=8, E=formula) p142
%         Y_1(i)=1/(Y(i)*exp(Y(i)))*((Y(i)^2+2.3347*Y(i)+0.2506)/(Y(i)^2+...
%             3.3307*Y(i)+1.6815));
%     elseif (Y(i)<0.01)
%         Y_1(i)=-log(Y(i))-0.5772;
%     else
%         Y_1(i)=0;
%     end
% end
% 
% % cable environment partial transient
% Theta_Et=((WI*Rho_sol)/(4*pi))*(X_1-Y_1);

X_1=expint(X);

Y_1=expint(Y);

% cable environment partial transient
Theta_Et=((WI*Rho_sol)/(4*pi))*(X_1-Y_1);

end