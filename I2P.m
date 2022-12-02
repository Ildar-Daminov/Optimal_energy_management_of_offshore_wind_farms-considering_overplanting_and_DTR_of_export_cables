function [S]=I2P(I)
%% Purpose of the function
% This function converts I (in A) to power values (in MVA) for submarine
% cable 225 kV
Unom=225000; %line-to-line voltage of submarine cable
S=I*Unom*sqrt(3)/1000000;  % cos phi=1 power in MW
end