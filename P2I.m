function [I]=P2I(S)
%% Purpose of the function
% This function converts power S (in MVA) to current I (in A) for submarine
% cable 225 kV
Unom=225000; %line-to-line voltage of submarine cable, V
I=S*1000000/(Unom*sqrt(3));  % cos phi=1 current in A

end