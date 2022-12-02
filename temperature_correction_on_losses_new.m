function [theta_VECTOR_new,MAE]=temperature_correction_on_losses_new(theta_VECTOR,f,ks,kp,D_ame,s,I_ech,duration,...
    timestep,Decalage_down,Rso,Rao,Xcable,c,Dmoy_Arm,w,stainless_steel)
%% What does this dunction do? 
% The losses of cable change as a function of cable temperature. Thus, 
% this function adjusts the final profile of temperature rises as a 
% function of losses change due to a temperature variations. 
%% Global variables (needed from other functions)
global  Tamb R_AC90_eval Lambda1 Lambda2 Wd T1 T2 T3 T4 R_DC20 Alpha_E_Al
global  Theta_inf_dielectric memoring days optimization thermal_memory t preloading
global theta_VECTOR_global TEMPERATURE_D_1  Alpha_E_Pb Alpha_E_AG thermal_memory_preload

%% Starting the adjustments 
% Save previous theta_VECTOR
theta_VECTOR_previous=theta_VECTOR;

if optimization==1
    if isempty(thermal_memory) % no thermal memory
        theta_VECTOR=theta_VECTOR;
    else % there is a thermal memory
        theta_VECTOR=theta_VECTOR+thermal_memory(days*96-95:days*96,1); % J day
    end
end
if memoring==1 % if memoring is activated
    if isempty(thermal_memory) % if thermal_memory does not exist
        %         do nothing
    else % there is already a thermal memory
        theta_VECTOR=theta_VECTOR+thermal_memory(days*96-95:end,1); % save transients for period [studied day : end of horizon]
    end
end
if isempty(optimization)&& isempty(memoring)&& ~isempty(thermal_memory_preload)
    theta_VECTOR=theta_VECTOR+thermal_memory_preload; % thermal_memory_preload should be at studied horizon!
end

% Temperature profile 
Temp_profile=Tamb+Theta_inf_dielectric+theta_VECTOR-273.15;

% The d.c. resistance per unit length of the conductor (Ω/m)at given
% operating temperature θ is given by (p14 IEC 60287-1-1):
R_DC=R_DC20.*(1+Alpha_E_Al.*(Temp_profile-20));
R_DC90=R_DC20.*(1+Alpha_E_Al.*(90-20))

Difference=R_DC90-R_DC;
% plot([Difference])
% ylabel('R DC difference')

% 2.1.2 Skin effect factor Ys(p14 IEC 60287-1-1)
Xs=sqrt(8*pi*f./R_DC.*ks*1e-7);
Xs90=sqrt(((8*pi*f)./(R_DC90)).*ks*1e-7)

Difference=Xs90-Xs;
% plot([Difference])
% ylabel('Xs difference')

% The above formula is accurate providing xs does not exceed 2,8, and
% therefore applies to the majority of practical cases.
if Xs>=2.8
    error('Huston, we have a problem : Xs>=2.8')
end

% The skin effect factor ys is given by (p14 IEC 60287-1-1):
Ys=Xs.^4./(192+0.8.*Xs.^4);
Ys90=Xs90.^4./(192+0.8.*Xs90.^4)

Difference=Ys90-Ys;
% plot([Difference])
% ylabel('Ys difference')

% The Ys formula is accurate providing Xs does not exceed 2.8, and
% therefore applies to the majority of practical cases.

% 2.1.4 Proximity effect factor yp for three-core cables and for three
% single-core cables (p15 IEC 60287-1-1):
% 2.1.4.1 Circular conductor cables
Xp=sqrt(8*pi*f./R_DC*kp*1e-7);
Xp90=sqrt(((8*pi*f)./(R_DC90))*kp*1e-7)

Difference=Xp90-Xp;
% plot([Difference])
% ylabel('Xp difference')

% The above formula is accurate providing xp does not exceed 2,8, and
% therefore applies to the majority of practical cases. (p15 IEC 60287-1-1)
if Xp>=2.8
    error('Huston, we have a problem: Xp>=2.8')
end

% Only if Xp<2.8                                                            else Xp=>2.8?
% The proximity effect factor is given by:
Yp=(Xp.^4./(192+0.8.*Xp.^4)).*(D_ame/s)^2.*(0.312.*((D_ame/s)^2)...                D_ame ou dc?? (D_ame)
    +1.18./(Xp.^4./(192+0.8.*Xp.^4)+0.27));
Yp90=(Xp90.^4./(192+0.8.*Xp90.^4)).*(D_ame/s)^2.*(0.312.*((D_ame/s)^2)...  D_ame ou dc?? (D_ame)
    +1.18./(Xp90.^4./(192+0.8.*Xp90.^4)+0.27))

Difference=Yp90-Yp;
% plot([Difference])
% ylabel('Yp difference')

% Yp=(Xp^4/(192+0.8*Xp^4))*(dc/s)^2*(0.312*((dc/s)^2)...                   D_ame ou dc?? (D_ame)
%     +1.18/(Xp^4/(192+0.8*Xp^4)+0.27));
% disp('Yp formula is changed. dc is used instead of D_ame')

% 2.3.10 Cables with each core in a separate lead sheath (SL type)&armoured
% (p24 IEC 60287-1-1)
% For a three-core cable of which each core has a separate lead sheath λ′′
% is zero and the loss factor for the sheaths is given by:

% 2.1 AC resistance of conductor (p13 IEC 60287-1-1)
% The a.c. resistance per unit length of the conductor at its maximum
% operating temperature is given by the following formula
R_AC=R_DC.*(1+Ys+Yp);
R_AC90=R_DC90.*(1+Ys90+Yp90)

Difference=R_AC90-R_AC;
% plot([Difference])
% ylabel('R_AC difference')

%% Lambda1 and Lambda2
% maximum operating temperature of sheath or screen is given by
% p17 IEC 60287-1-1:
Temp_sc = Temp_profile-(I_ech.^2.*R_AC+0.5*Wd).*T1;
% The resistance of the sheath or screen at its maximum operating
% temperature is given by p17 IEC 60287-1-1:
Rs = Rso.*(1+Alpha_E_Pb.*(Temp_sc-20));
% Ratio of the total losses in sheath respectively to the total
% conductor losses (or losses in one sheath to the losses in one
% conductor)(p24 IEC 60287-1-1)
Lambda1_new=(Rs./R_AC).*(1.5./(1+(Rs./Xcable).^2));
% NB: Here Lambda1 is actually Lambda1'. In general Lambda1 is
% calculated as Lambda1=Lambda1'+Lambda1'' (see p16  IEC 60287-1-1)but
% in accordance with p24  IEC 60287-1-1 Lambda1''==0. Thus, we assume
% that Lambda1=Lambda1';

% The maximum operating temperature of the armour is given by p25 IEC
% 60287-1-1:
Temp_AG = Temp_profile-((I_ech.^2.*R_AC+0.5*Wd).*T1+(I_ech.^2.*...
    R_AC.*(1+Lambda1_new)+Wd).*3.*T2);

% The resistance of the armour at its maximum operating temperature is
% given by p25 IEC 60287-1-1:
Rarm=Rao.*(1+Alpha_E_AG.*(Temp_AG-20));

% 2.4.2.5 SL type cables (p29 CEI 60287-1-1)
% Where the armour is over a SL type cable, the screening effect of the
% sheath currents reduces the armour loss. The formula for λ2 given in
%  2.4.2.3.1 shall be multiplied by the factor

Lambda1_dash=(Rs./R_AC).*1./(1+(Rs./Xcable).^2);

Factor=(1-(R_AC./Rs).*Lambda1_dash);
% where Lambda1_dash is obtained from 2.3.1.                            % !!! 2.3.1 is for Two single-core cables, and three single-core cables

% Ratio of the total losses in armour respectively to the total
% conductor  losses (or losses in one sheath or armour to the losses
% in one conductor) p28 CEI 60287-1-1

Lambda2=1.23.*(Rarm./R_AC).*(2.*c./Dmoy_Arm).^2.*1./(1+(2.77.*Rarm...
    .*1e6./(w)).^2);
Lambda2=Factor.*Lambda2;
if stainless_steel==1
    Lambda2=0; % if stainless steel is used
end
if stainless_steel==1
    disp('Stainless steel is assumed for armour. Lambda2=0 ie no losses in armour')
end
%% Superposition principle is applied in this "for cycle":
if preloading==1
    I_ech=I_ech(1:35040);
end
for k=1:length(I_ech) % for each current value
    if memoring==1 % if memoring status of IEC60853_2 is switched on
        if k<=(days*96+1)-(days*96-95) % continue k until 96 values are passed
            % Define the time offset for starting point
            Decalage=zeros(duration/timestep*k-duration/timestep,1);
            
            % Define I as k-th current value from I_ech
            I=I_ech(k);
            
            
            % Find cable and environment trainsient as well as attainement factor
            [Theta_ct,Theta_Et,Alpha_t]=losses_partial_transients_new(I,R_AC(k),Lambda1_new(k));
            
            % 4.4 Calculation of the complete temperature transient (p43 IEC853-2)
            % After calculating separately the two partial transients and the
            % conductor to cable surface attainment factor the total transient rise
            % above an ambient temperature is obtained:
            % a) for buried configurations by simple addition of the cable and
            %    modified environment partial transients;
            
            % the complete temperature transient (p43 IEC853-2)
            Theta=Theta_ct+Alpha_t.*Theta_Et;
            
            if k==1 % first current value in load profile
                % Shift the complete temperature transient for 1 hour ahead and
                % make it negative
                Theta_down=-[Decalage_down;Theta];
                % save only studied interval
                Theta_down=Theta_down(1:length(t),1);
                % save to matrix for down step current
                theta_VECTOR_new=Theta+Theta_down;
                
            else % second and next current value in load profile
                % Shift the complete temperature transient (up step) in time to its
                % start point
                Theta=[Decalage;Theta];
                % save only studied interval
                Theta=Theta(1:length(t),1);
                % Shift the complete temperature transient for 1 hour ahead and
                % make it negative
                Theta_down=-[Decalage_down;Theta];
                % save only studied interval
                Theta_down=Theta_down(1:length(t),1);
                % save to matrix for down step current
                theta_VECTOR_new=theta_VECTOR_new+Theta+Theta_down;
                
            end % end of if k == 1
        else % if k>(days*24+1)-(days*24-23)
            break % stop for cycle and go to next section of IEC60853_2.m
        end %  end of  if k<=(days*96+1)-(days*96-95)
        
    else % if memoring = []
        % Define the time offset for starting point
        Decalage=zeros(duration/timestep*k-duration/timestep,1);
        
        % Define I as k-th current value from I_ech
        I=I_ech(k);
        
        % Find cable and environment trainsient as well as attainement factor
        [Theta_ct,Theta_Et,Alpha_t]=losses_partial_transients_new(I,R_AC(k),Lambda1_new(k));
        %     [Theta_ct,Theta_Et,Alpha_t]=losses_partial_transients_new(I,R_AC(k),Lambda1);
        %     [Theta_ct90,Theta_Et90,Alpha_t90]=losses_partial_transients_new(I,R_AC90,Lambda1);
        
        % 4.4 Calculation of the complete temperature transient (p43 IEC853-2)
        % After calculating separately the two partial transients and the
        % conductor to cable surface attainment factor the total transient rise
        % above an ambient temperature is obtained:
        % a) for buried configurations by simple addition of the cable and
        %    modified environment partial transients;
        
        % the complete temperature transient (p43 IEC853-2)
        Theta=Theta_ct+Alpha_t.*Theta_Et;
        %     Theta90=Theta_ct90+Alpha_t90.*Theta_Et90;
        
        if k==1 % first current value in load profile
            % Shift the complete temperature transient for 1 hour ahead and
            % make it negative
            Theta_down=-[Decalage_down;Theta];
            %         Theta_down90=-[Decalage_down;Theta90];
            % save only studied interval
            Theta_down=Theta_down(1:length(t),1);
            %         Theta_down90=Theta_down90(1:length(t),1);
            
            % save to matrix for down step current
            theta_VECTOR_new=Theta+Theta_down;
            %         theta_VECTOR90=Theta90+Theta_down90;
            
        else % second and next current value in load profile
            % Shift the complete temperature transient (up step) in time to its
            % start point
            Theta=[Decalage;Theta];
            %         Theta90=[Decalage;Theta90];
            % save only studied interval
            Theta=Theta(1:length(t),1);
            %         Theta90=Theta90(1:length(t),1);
            
            % Shift the complete temperature transient for 1 hour ahead and
            % make it negative
            Theta_down=-[Decalage_down;Theta];
            %         Theta_down90=-[Decalage_down;Theta90];
            
            % save only studied interval
            Theta_down=Theta_down(1:length(t),1);
            %         Theta_down90=Theta_down90(1:length(t),1);
            % save to matrix for down step current
            theta_VECTOR_new=theta_VECTOR_new+Theta+Theta_down;
            %         theta_VECTOR90=theta_VECTOR90+Theta90+Theta_down90;
            
        end % end of if k == 1
    end % end of if memoring==1
end % go to next k or end of for cycle

% Difference=theta_VECTOR90-theta_VECTOR_new;
% plot([Difference])
% ylabel('theta_VECTOR difference')
if optimization==1
    if isempty(thermal_memory) % no thermal memory
        theta_VECTOR2=theta_VECTOR_new;
    else % there is a thermal memory
        theta_VECTOR2=theta_VECTOR_new+thermal_memory(days*96-95:days*96,1); % J day
    end
    Temp_profile_new=Tamb+Theta_inf_dielectric+theta_VECTOR2-273.15;
    %     if min(Difference)<0 && max(Temp_profile_new)<=90
    %         disp('some value in theta_VECTOR is higher than correponding value theta_VECTOR90')
    %     end
end
% if memoring==1 % if memoring is activated
%     if isempty(thermal_memory) % if thermal_memory does not exist
% %         do nothing
%     else % there is already a thermal memory
%         theta_VECTOR_new=theta_VECTOR_new+thermal_memory(days*96-95:end,1); % save transients for period [studied day : end of horizon]
%     end
% end

%% Mean average error
MAE=abs(theta_VECTOR_new-theta_VECTOR_previous); % absolute value
MAE=mean(MAE); % mean value

if isnan(MAE)
    error('MAE is NaN')
end
end