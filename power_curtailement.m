function [I_output]=power_curtailement(I_input)
%% Global variables
global optimization  days
%% What does this function do?
% This function simulates the power curtailements of offshore wind farm in
% case if cable's limit of temperature (90 degC) is exceeded

% Input
% I_input - current profile (in A)

% Output
% I_output _ current profile (in A) respecting the temperature constraint
%% Constants
Tmax_limit=90;  % Tmax limit is 90 deg C
delta=0.01;     % Delta for current reduction
iter=0;         % create a variable with zero value (needed for monitoring only)

%% Algorithm of power curtailement

% Change a variable
I=I_input;

% Save I as I_init
I_init=I;

% turn on the mode "optimization" in cable_thermal_model_IEC_60853_2.m
optimization=1;

% Calculate Tmax and T
[Tmax,T,~]=cable_thermal_model_IEC_60853_2(I);
disp('cable_thermal_model_IEC_60853_2 is used')

% delete T at 00:00;
T(1,:)=[];

% Save T as T_init
T_init=T;
Tmax=max(T);
Tmax_init=Tmax;
idx=1; % index of time instant where the temperature exceeds 90 degC
        % and when curtailement is applied

% Checking if Tmax is higher than T
if Tmax>Tmax_limit % if Tmax is greater than Tmax limit
    disp('Violation of temperature is found. Launching the algorithm...')
    
    while ~(Tmax<=Tmax_limit)% while temperature is higher than limit
        % Iteration monitoring
        iter=iter+1;
        
        % Save previous I profile
        I_old=I;
        
        % Find index of values greater than 90 degC
        idx_T=find(T>Tmax_limit);
        
        % Reduce I at T_descend(1,2) for I_delta
        I(idx_T(idx)+days*96-96)=I(idx_T(idx)+days*96-96)*(1-delta);
        
        if abs(I_old(idx_T(idx)+days*96-96)-I(idx_T(idx)+days*96-96))<0.1 % checking if the algorithm converges to the same solution
            while abs(I_old(idx_T(idx)+days*96-96)-I(idx_T(idx))+days*96-96)<0.1
                idx=idx+1; % shiting time of power shedding
                I(idx_T(idx)+days*96-96)=I(idx_T(idx)+days*96-96)*(1-delta); % power shedding
            end
        else % if algorithm converged to another solution
            idx=1; % return time to zero (power shedding shift)
        end
        % Calculate the thermal regime for corrected load profile
        [~,T,~]=cable_thermal_model_IEC_60853_2(I);
        disp('cable_thermal_model_IEC_60853_2 is used')
        
        % delete T==0
        T(1,:)=[];
        % Max T
        Tmax=max(T);

    end % end of while cycle
    
else % if Tmax<Tmax_limit
    disp('Temperature is OK ');
    
end % end if Tmax>Tmax_limit

% Results
I_output=I;

optimization=[]; % turn off the mode "optimization" in cable_thermal_model_IEC_60853_2.m

% figure('DefaultAxesFontSize',14)
% yyaxis left
% plot([I_init I_output],'linewidth',2)
% ylabel('Current, A')
% ylim([0 1500])
% yyaxis right
% plot([T_init T_output],'linewidth',2)
% ylabel('Temperature, degC')
% legend('I init','I optim','Temperature init', 'Temeprature optim' )
% xlabel('Time')
end % end of function
