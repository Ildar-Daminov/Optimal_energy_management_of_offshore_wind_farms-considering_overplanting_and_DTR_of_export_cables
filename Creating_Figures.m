clc
clear all
close all

%% Purpose
%  This script generates figures used in the paper [1]

% [1] I. Daminov et al, Optimal energy management of offshore wind farms 
% considering the combination of overplanting and dynamic rating in 
% CIGRE Paris proceedings, France  2022
%% Figure 1  Illustration of the electrical infrastructure of OWF. In our
% study, only the submarine section of export cable is considered

% This figure was ploted without using MATLAB
%% Figure 2 Equivalent thermal circuit with two cells

% This figure was ploted without using MATLAB
%% Figure 3 Power generation forecasts and actual measurements of the OWFs 
% connected to Belgian network

% Clear workspace
clear all

% Load the power profiles from Elia
load('Elia_Jan13_2016_Jan13_2017_powers.mat')

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Plot load factors of offshore wind farm(s)
plot(t_year,[Load_factor90_year,Load_factor50_year,Load_factor10_year,Load_factor_measur_year])

% Plot the legend 
legend('P90','P50','P10','Measured')

% Display the ylabel 
ylabel('Load Factor, pu')

% The Figure reprsents annual load factors o offshore wind farms (both 
% forecast and measurements. To see the Figure 3 as it is done in paper, just 
% zoom in the power profile at the needed interval  

%% Figure 4 Market prices in France and different situations where the DA 
%% price may be greater, inbetween or less than the imbalance prices

% Clear workspace
clear all

% Load day-ahead and imbalance prices
load('data_DA_IMB_years_2015 2020.mat')

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'WindowState','maximized');

% Plot day-ahead price and imbalance prices+-
plot(time_2018,[DA_2018,Cb_plus_2018,Cb_minus_2018])

% Plot the legend 
legend('Day-ahead','Imbalance+','Imbalance-')

% Display the ylabel 
ylabel('Price, EUR/MWh')

%% Figure 5 Annual revenue as a function of the overplanting rate for the 
% STR and DTR cases. In all cases, the P50 strategy is used as the 
% commitment strategy

% Clear workspace
clear all

% Extract and process the results from main_simulations_case.mat
for cases=1:2 % case 1 = STR case 2 = DTR

    % Load the data for given case 
    filename=sprintf('main_simulations_case%d.mat',cases);
    load(filename)

   
    if cases==1 % if STR

        for capacity_idx=1:length(Overplanting_capacity) % for each overplanting rate
            
            % Extract a daily revenue for P50 strategy  
            eval(['P50_Revenues_case1=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_P50_Elia_2018/4']); % /4 for MWh
           
            % Calculate the annual revenue for P50 strategy
            eval(['Annual_P50_Revenue_case1_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(P50_Revenues_case1)']);

            % Extract a daily revenue for Optimal power profile strategy  
            eval(['Optimal_power_profile_Revenues_case1=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Revenue_opt']);
            Optimal_power_profile_Revenues_case1=cell2mat(Optimal_power_profile_Revenues_case1)/4; % /4 for MWh
            
            % Calculate the annual revenue for Optimal power profile strategy 
            eval(['Annual_Optimal_power_profile_Revenues_case1_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Optimal_power_profile_Revenues_case1);']);
            
            % Note: 
            % - "P50 strategy": Using P50 over the entire year (Reference for industry, “business as usual�?)
            % - "Optimal power profile":Using optimal power profile st wind installed capacity. Optimal power profile is calculated by fmincon
            
        end % end of for cycle 

    elseif cases==2 % if DTR 

        for capacity_idx=1:length(Overplanting_capacity) % for each overplanting rate
            % Extract a daily revenue for P50 strategy  
            eval(['P50_Revenues_case2=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_P50_Elia_2018/4']); % /4 for MWh
            
            % Calculate the annual revenue for P50 strategy
            eval(['Annual_P50_Revenue_case2_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(P50_Revenues_case2)']);
            
            % Extract a daily revenue for Optimal power profile strategy  
            eval(['Optimal_power_profile_Revenues_case2=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Revenue_opt']); 
            Optimal_power_profile_Revenues_case2=cell2mat(Optimal_power_profile_Revenues_case2)/4;% /4 for MWh
            
            % Calculate the annual revenue for Optimal power profile strategy 
            eval(['Annual_Optimal_power_profile_Revenues_case2_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW=sum(Optimal_power_profile_Revenues_case2)']);

        end % end of for cycle
    end
end

% case1 (STR): preparing the annual revenue in % to reference case (without
% overplanting - Annual_P50_Revenue_case1_339_MW)
case1_STR=[Annual_P50_Revenue_case1_339_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_339_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_373_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_373_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_407_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_407_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_441_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_441_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_509_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_509_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_577_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_577_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_645_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_645_MW/Annual_P50_Revenue_case1_339_MW*100;...
    Annual_P50_Revenue_case1_679_MW/Annual_P50_Revenue_case1_339_MW*100,Annual_Optimal_power_profile_Revenues_case1_679_MW/Annual_P50_Revenue_case1_339_MW*100;];

% case2(DTR): preparing the annual revenue in % to reference case (without
% overplanting - Annual_P50_Revenue_case2_339_MW)
case2_DTR=[Annual_P50_Revenue_case2_339_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_339_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_373_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_373_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_407_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_407_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_441_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_441_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_509_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_509_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_577_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_577_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_645_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_645_MW/Annual_P50_Revenue_case2_339_MW*100;...
    Annual_P50_Revenue_case2_679_MW/Annual_P50_Revenue_case2_339_MW*100,Annual_Optimal_power_profile_Revenues_case2_679_MW/Annual_P50_Revenue_case2_339_MW*100;];

% create a figure 
fig5=figure('DefaultAxesFontSize',14)

% Plot bars (for P50 strategy only)
H=bar([case1_STR(:,1),case2_DTR(:,1)])


% Showing the bar values
hT=[];              % placeholder for text object handles
for ii=1:length(H)  % iterate over number of bar objects
    hT=[hT text(H(ii).XData+H(ii).XOffset,H(ii).YData,num2str(round(H(ii).YData.')), ...
        'VerticalAlignment','bottom','horizontalalign','center')];
end

% Change tick lables in x-axis
xticklabels({'1','1.1','1.2','1.3','1.5','1.7','1.9','2'})

% Show legend
legend('STR', 'DTR')

% Show x/y labels 
ylabel('Annual Revenue,reference %')
xlabel('Overplanting rate,pu')

% Adjusting the units and position of figures 
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

% Table 2 in paper 
Table_2=[case1_STR(:,1)';case1_STR(:,2)';case2_DTR(:,1)';case2_DTR(:,2)']
%% Figure 6 Annual generated and curtailed energy (normalized to the annual
% generated energy of a non-overplanted OWF) a function of the overplanting 
% rate (with both STR and DTR cable limits)

% Clear workspace
clear all

for cases=1:2 % case 1 = STR case 2 = DTR

    % Load the results of simulations
    filename=sprintf('main_simulations_case%d.mat',cases);
    load(filename)

    if cases==1 % if STR 

        % Extracting and processing the simulation results 
        for capacity_idx=1:length(Overplanting_capacity) % for each overplanting rate 
            
            % Prepare the variables  
            Pactual=NaN;
            Pmeasur=NaN;

            for days=1:365 % for each day 

                % Extact actual production of offshore wind farm
                % (considering cable constraint) for given day
                eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.Pactual(days);']);
                
                % Extact theoretically-possible production of offshore wind farm
                % (without considering cable constraints) for given day
                eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.CurrentConstraints.P_measur_day(days);']);
                
                % Extract power profile from cell and convert it to mat
                Pactual_interm=cell2mat(Pactual_interm);
                Pmeasur_interm=cell2mat(Pmeasur_interm);

                % Append the power profiles for the annual vectors
                Pactual=[Pactual;Pactual_interm];
                Pmeasur=[Pmeasur;Pmeasur_interm];
            end % end of for cycle 
    
            % Delete NAN values at the beginning 
            Pactual(1,:)=[];
            Pmeasur(1,:)=[];

            % Calculate the energy production 
            E_measur=trapz(Pmeasur)/4; % /4 for MWh
            E_actual=trapz(Pactual)/4; % /4 for MWh
            
            % Calculate the energy curtailement 
            E_curt=E_measur-E_actual;

            % Find the power profile for curtailements 
            Curtailements=Pmeasur-Pactual;
            
            % Save data
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pmeasur=Pmeasur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pactual=Pactual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Curtailements=Curtailements;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_measur=E_measur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual=E_actual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_curt=E_curt;'])
        
        end % end of for cycle 

    elseif cases==2 % if DTR


        for capacity_idx=1:length(Overplanting_capacity) % for each overplanting rate 
            
            % Prepare the variables  
            Pactual=NaN;
            Pmeasur=NaN;


            for days=1:365 % for each day

                % Extact actual production of offshore wind farm
                % (considering cable constraint) for given day
                eval(['Pactual_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.Pactual(days);']);
                
                % Extact theoretically-possible production of offshore wind farm
                % (without considering cable constraints) for given day
                eval(['Pmeasur_interm=Overplanting_results.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.TempConstraints.P_measur_day(days);']);
                
                % Extract power profile from cell and convert it to mat
                Pactual_interm=cell2mat(Pactual_interm);
                Pmeasur_interm=cell2mat(Pmeasur_interm);

                % Append the power profiles for the annual vectors
                Pactual=[Pactual;Pactual_interm];
                Pmeasur=[Pmeasur;Pmeasur_interm];
            end % end of for cycle 

            % Delete NAN values at the beginning 
            Pactual(1,:)=[];
            Pmeasur(1,:)=[];

            % Calculate the energy production 
            E_measur=trapz(Pmeasur)/4; % /4 for MWh
            E_actual=trapz(Pactual)/4; % /4 for MWh

            % Calculate the energy curtailement 
            E_curt=E_measur-E_actual;

            % Find the power profile for curtailements 
            Curtailements=Pmeasur-Pactual;

            % Save data
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pmeasur=Pmeasur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Pactual=Pactual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.Curtailements=Curtailements;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_measur=E_measur;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_actual=E_actual;'])
            eval(['PandE.case' num2str(cases) '.Capacity_' num2str(round(Overplanting_capacity(capacity_idx))) '_MW.E_curt=E_curt;'])
        end % end of for cycle
    end % end of if case
end % end of for cylce 


% Caclulating the energy actually generated (in % to reference) at each overplanting rate 
data_additional_Energy_generated=[PandE.case1.Capacity_339_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,... % case 1 (STR)
    PandE.case1.Capacity_373_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_407_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_441_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_509_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_577_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_645_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case1.Capacity_679_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100;... % note ; in this line
    PandE.case2.Capacity_339_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,... % case 2 (DTR)
    PandE.case2.Capacity_373_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_407_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_441_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_509_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_577_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_645_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100,...
    PandE.case2.Capacity_679_MW.E_actual/PandE.case1.Capacity_339_MW.E_actual*100];


% Caclulating the energy curtailed (in % to reference) at each overplanting rate 
data_all_Energy_curtail2=[PandE.case1.Capacity_339_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,... % case 1 (STR)
    PandE.case1.Capacity_373_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_407_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_441_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_509_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_577_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_645_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case1.Capacity_679_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100;...  % note ; in this line
    PandE.case2.Capacity_339_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...  % case 2 (DTR)
    PandE.case2.Capacity_373_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_407_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_441_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_509_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_577_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_645_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100,...
    PandE.case2.Capacity_679_MW.E_curt/PandE.case1.Capacity_339_MW.E_measur*100];


% Create a figure 
figure('DefaultAxesFontSize',14) % stacked 

% Data for plotting the bars 
nbars = 8;
x = 1;
xpoints=NaN;

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on

for i = 1:nbars % for each bars 

    % Prepate the data for plotting STR
    data_bar = [data_additional_Energy_generated(1,i);data_all_Energy_curtail2(1,i)];
    
    % Find x point 
    xpoint=x+3;
    x=xpoint;

    % Plot the bars
    H=bar(xpoint, data_bar,'stacked','b');

    % Showing the bar colors 
    H(2).FaceColor = 'flat';
    H(2).CData = [0.5 0.5 0.5];
    hT=[];              % placeholder for text object handles

    % Showing the bar values 
    for ii=1:length(H)  % iterate over number of bar objects
        if rem(ii,2)==0 % even number (2,4,6 etc) corresponding to curtailement bar 
            hT=[hT text(H(ii).XData+H(ii).XOffset,H(ii).YData+H(ii-1).YData-5,num2str(round(H(ii).YData.')), ...
                'VerticalAlignment','bottom','horizontalalign','center')];
        else % odd number 1,3,5 etc corresponding to energy generated bar 
            hT=[hT text(H(ii).XData+H(ii).XOffset,H(ii).YData-5,num2str(round(H(ii).YData.')), ...
                'VerticalAlignment','bottom','horizontalalign','center')];        
        end
    end

    % Prepate the data for plotting DTR
    data_bar = [data_additional_Energy_generated(2,i);data_all_Energy_curtail2(2,i)];
    
    % Find x point for the bar 
    xpoint=x+groupwidth;

    % Plot the bar
    H=bar(xpoint, data_bar,'stacked','FaceColor',[0.8500 0.3250 0.0980]);

    % Change the color 
    H(2).FaceColor = 'flat';
    H(2).CData = [0 0 0];
    hT=[];              % placeholder for text object handles

    % Showing the bar values 
    for ii=1:length(H)  % iterate over number of bar objects
        if rem(ii,2)==0 % even number (2,4,6 etc) corresponding to curtailement bar 
            hT=[hT text(H(ii).XData+H(ii).XOffset,H(ii).YData+H(ii-1).YData,num2str(round(H(ii).YData.')), ...
                'VerticalAlignment','bottom','horizontalalign','center')];
        else % odd number 1,3,5 etc corresponding to energy generated bar 
            hT=[hT text(H(ii).XData+H(ii).XOffset,H(ii).YData-5,num2str(round(H(ii).YData.')), ...
                'VerticalAlignment','bottom','horizontalalign','center')];        
        end
    end

    % Setting the x tick labels 
    xticks([xpoint]);
    xticklabels({num2str(Overplanting_rate(i))})

    % Save x points 
    xpoints(end+1)=xpoint;

end

% Delete first NaN
xpoints(:,1)=[];

% Adujst the xpoints 
xticks([xpoints-0.1]);

% Set the x tick labels 
xticklabels({'1','1.1','1.2','1.3','1.5','1.7','1.9','2'})

% Set axes labels
ylabel('Energy,% of reference')
xlabel('Overplanting rate')

% Show the legend
legend('Generated energy: Current', 'Curtailed Energy: Current','Generated energy: Temperature','Curtailed energy: Temperature','Location','northwest')

% Adjusting the units and position of figures 
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

%% Table 2 Annual revenue (in %) as a function of overplanting rate for the
% two commitment strategies with STR and DTR

%  For these results see Figure 5