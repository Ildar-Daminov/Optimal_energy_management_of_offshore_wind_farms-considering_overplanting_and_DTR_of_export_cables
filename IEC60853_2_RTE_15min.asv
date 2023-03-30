function [T_max,Temperature,I_adm]=cable_thermal_model_IEC_60853_2(Iech_load)
%% ************************************************************************
%--------------------------------------------------------------------------
% Input of cable_thermal_model_IEC_60853_2.m:
% Iech_load   - hourly current profile in A;

% Output of cable_thermal_model_IEC_60853_2.m:
% Temperature - hourly profile of cable temperature in degC from t=00:00 to
%               t_end=length(Iech_load)e.g.[25x1] 25 since 1 value for t=0;
% T_max       - maximal temperature of cable conductor, [1x1];
% I_adm       - steady-state current, calculated in accordance with IEC287
%               [1x1];
%% ***********************************************************************

%% ************************************************************************
%% Electrothermal modelling of submarine cable
%% ************************************************************************
% Global variables (used in different functions)
global S_ame_pertes Lcable Rho_0_E_Al Alpha_E_Al Tinit Lambda1 Lambda2 Alpha_E_Pb
global T1 T2 T3 T4 TA TB Ta Tb Qc Qi Qs Qf Qa Qj QA QB Wd t a b I_D_1
global pVW  Depth De Tamb Rho_sol Diff_sol R_DC20 R_AC90_eval Ys Yp Alpha_E_AG
global Theta_inf_dielectric optimization memoring theta_VECTOR_ISGT days
global thermal_memory theta_VECTOR_global thermal_memory_preload TEMPERATURE_D_1 preloading

%% ************************************************************************
%% Initial temperature of cable (dielectric losses are added later!)
%% ************************************************************************
Tamb=273.15+18;         % Ambient temperature in Kelvin (K)

% p10 IEC853-2: it has been assumed that the voltage had been applied for
% a sufficiently long time for the conductor temperature rise due to
% dielectric loss to have reached a steady state. The total temperature
% rise of the conductor above ambient is then the sum of the steady state
% temperature rise due to the dielectric loss (as given by IEC 60287) and
% the transient temperature via riations due to change of current.

% Tinit=Tamb;             % Initial temperature of cable in Kelvin (K)
% Tinit is updated in the end of the script!!!

% Here, Tinit = Tamb but it is also neccesary to add  the temperature rise
% due to dielectric losses. The dielectirc losses and temperature rise
% could be calculated by IEC 60287. But to do that , thermal parameters of
% cable should be firstly calculated. (see next parts of the script)

% Temperature rise from dielectric losses is be added in the end of script

%% ************************************************************************
%% System constants
%% ************************************************************************
U_0=225e3/sqrt(3) ;     % Phase to neutral voltage in V, 130 kV
f=50;                   % Frequency in Hz
w=2*pi*f;               % Network pulsation

%% ************************************************************************
%% Cable characteristics and geometry
%% ************************************************************************
Lcable=1 ;                   % cable length (m)
%--------------------------------------------------------------------------
% Core / conductor [for 1 core of three-core cable]
%--------------------------------------------------------------------------
D_ame =38.2e-3;              % core diameter                    (m)
rame=D_ame/2;                % core radius                      (m)
% S_ame_pertes=pi*rame^2;      % actual conductor cross-section   (m2)
S_ame_pertes=1000e-6;        % nominal conductor cross-section  (m2)

% Note: We use 1000 mm2 instead of 1146 mm2 because the rest is air
% (which does not transmit current). This is also done in the COMSOL model.
%--------------------------------------------------------------------------
% Semiconductor 1 (Semiconducting compound)
%--------------------------------------------------------------------------
Thick_SC1=1.5e-3;            % thickness of semiconductor (m)
D_SC1_int=D_ame;             % internal diameter          (m)
D_SC1_ext=D_ame+2*Thick_SC1; % external diameter          (m)

%--------------------------------------------------------------------------
% Insulation (XLPE)
%--------------------------------------------------------------------------
Thick_Isol=22e-3;                   % thickness of insulation (m)
D_Isol_int=D_SC1_ext;               % internal diameter       (m)
D_Isol_ext=D_SC1_ext+2*Thick_Isol;  % external diameter       (m)

%--------------------------------------------------------------------------
% Semiconductor 2 (Semiconducting compound)
%--------------------------------------------------------------------------
Thick_SC2=1.5e-3;                   % thickness of semiconductor(m)
D_SC2_int=D_Isol_ext;               % internal diameter         (m)
D_SC2_ext=D_Isol_ext+2*Thick_SC2;   % external diameter         (m)

%--------------------------------------------------------------------------
% Sealing mats (Semiconducting water-swellable tapes)
%--------------------------------------------------------------------------
Thick_Matelas=1e-3;                      % thickness of sealing mats(m)
D_Matelas_int=D_SC2_ext;                 % internal diameter        (m)
D_Matelas_ext=D_SC2_ext+2*Thick_Matelas; % external diameter        (m)

%--------------------------------------------------------------------------
% Metallic screen (Lead alloy)
%--------------------------------------------------------------------------
Thick_Pb_Screen=2.5e-3;                          % thickness of sreen(m)
D_Pb_Screen_int=D_Matelas_ext;                   % internal diameter (m)
D_Pb_Screen_ext=D_Matelas_ext+2*Thick_Pb_Screen; % external diameter (m)
D_moy_screen=D_Matelas_ext+Thick_Pb_Screen;      % average diameter  (m)

%--------------------------------------------------------------------------
% Internal protection sheath (Semiconducting PE)
%--------------------------------------------------------------------------
Thick_Gaine_PH=3e-3;                            % sheath thickness  (m)
D_Gaine_PH_int=D_Pb_Screen_ext;                 % internal diameter (m)
D_Gaine_PH_ext=D_Gaine_PH_int+2*Thick_Gaine_PH; % external diameter (m)

%--------------------------------------------------------------------------
% Filler (Plastic shaped elements)
%--------------------------------------------------------------------------
Dext_filler = 228e-3;      % D_Armure_Matelas_int (below)
%--------------------------------------------------------------------------
% Armour mats (PP yarn)
%--------------------------------------------------------------------------
Thick_Armure_Matelas=0.5e-3;                % thickness of armour mats (m)
D_Armure_Matelas_int=228e-3;                % internal diameter        (m)
D_Armure_Matelas_ext=D_Armure_Matelas_int...% external diameter        (m)
    +2*Thick_Armure_Matelas;

%--------------------------------------------------------------------------
% Armour (Steel)
%--------------------------------------------------------------------------
Ep_arm=6e-3;                   % armour thickness       (m)
Dfil_Arm=Ep_arm;               % diameter of round wire (m)
Dint_arm=D_Armure_Matelas_ext; % internal diameter      (m)
Dext_arm=Dint_arm+2*Ep_arm;    % external diameter      (m)
Dmoy_Arm=Dint_arm+1*Ep_arm;    % average diameter       (m)
Nfil_Arm=99;                   % wire number in armour  (-)
%--------------------------------------------------------------------------
% Outer serving (PP yarn)
%--------------------------------------------------------------------------
Ep_Gaine_Ext=5e-3;                      % thickness of outer serving (m)
Dgaine_ext_int=Dext_arm;                % internal diameter          (m)
Dgaine_ext_ext=Dext_arm+2*Ep_Gaine_Ext; % external diameter          (m)

%% ************************************************************************
%% Physical parameters of the materials (electrical and thermal)
%% ************************************************************************

%--------------------------------------------------------------------------
% Aluminium
%--------------------------------------------------------------------------
Rho_0_E_Al=2.8264e-8; % electrical resistivity of Alu at 20°C % (p30 Table1
                      % 60287-1-1)
Alpha_E_Al=4.03e-3;   % temperature coefficient(p30 Table1 IEC 60287-1-1)
Cth_Al=2.5e6;         % volumetric specific heat (J/m3.K)(p77 IEC 60853-2)
%--------------------------------------------------------------------------
% Semiconductor1 (Semiconducting compound)
%--------------------------------------------------------------------------
% Note: p19 IEC 60287-2-1:for thermal calculations semi-conducting layers
%(including metallized carbon paper tapes) are considered as part of
% the insulation.
%--------------------------------------------------------------------------
% Insulation - XLPE
%--------------------------------------------------------------------------
Rho_th_isol=3.5;      % thermal resistivity (K.m/W)      p75 IEC 60853-2
Cth_Isol=2.4e6;       % volumetric specific heat(J/m3.K) p75 IEC 60853-2
eps_r=2.5;            % relative permittivity (p64 CEI 60287-1-1)for
                      % cables greater than 18/30 (36) kV cables (unfilled)
tan_delta=0.001;      % loss factors (p64 CEI 60287-1-1)for
                      % cables greater than 18/30 (36) kV cables(unfilled)*

%*Safe values at maximum permissible temperature, applicable to the highest
% voltages normally specified for each type of cable.

%--------------------------------------------------------------------------
% SemiConductor2 (Semiconducting compound)
%--------------------------------------------------------------------------
% Note: p19 IEC 60287-2-1:for thermal calculations semi-conducting layers
%(including metallized carbon paper tapes) are considered as part of
% the insulation.
%--------------------------------------------------------------------------
% Water-swellable tape
%--------------------------------------------------------------------------
% Note: p19 IEC 60287-2-1:for thermal calculations semi-conducting layers
%(including metallized carbon paper tapes) are considered as part of
% the insulation.

% Note: Swelling tape is assumed to be a part of the insulation since its
% thickness is small, and it is assumed to have the same thermal
% resistivity. (p6 Nielsen, Applied Sciences 2019)

%--------------------------------------------------------------------------
% Screen (lead alloy)
%--------------------------------------------------------------------------
Cth_ecran_Pb=1.45e6; % 1.7e6; % volumetric specific heat(J/m3.K)

Rho_0_E_Pb=21.4e-8;         % electrical resistivity of lead alloy at20°C
                            % (p30 Table1 60287-1-1)

Alpha_E_Pb=4.0e-3;          % temperature coefficient of lead alloy at20°C
                            % (p30 Table1 IEC 60287-1-1)

% Note: Metallic layers are neglected, as their thermal resistance is
% negligible compared with that of polycomposite materials.
% (p6 Nielsen, Applied Sciences 2019)

% Another note: the thermal resistance of the metallic parts in the cable,
% even though not equal to zero, is so small that it is usually neglected
% in rating computations (p34 Anders Rating of electric power cables 1997)

% (IEC does not provide instructions about thermal resistances of
%  metallic parts?)

%--------------------------------------------------------------------------
% Internal sheath (Semiconducting PE)
%--------------------------------------------------------------------------
% Note: p19 IEC 60287-2-1:for thermal calculations semi-conducting layers
%(including metallized carbon paper tapes) are considered as part of
% the insulation.

Rho_th_gaine_PH=3.5;          % thermal resistivity (K.m/W)
% p75 IEC 60853-2
Cth_Gaine_PH=2.4e6;%1.7e6;    % volumetric specific heat(J/m3.K)
% p75 IEC 60853-2
%--------------------------------------------------------------------------
% Armour matelas (polypropylene yarn)
%--------------------------------------------------------------------------
Rho_th_armure_matelas=5.5;    % thermal resistivity (K.m/W)
% p55 IEC 60287-2-1
Cth_armure_matelas=1.8e6;     % volumetric specific heat(J/m3.K)
% not given in IEC. taken from p8 Nielsen Applied Sciences 2019

%--------------------------------------------------------------------------
% Filler (Plastic shaped elements)
%--------------------------------------------------------------------------
% Kappa_filler=Kappa_Gaine_tri_int*0.5;   % 0.5:article EWTEC 2017 (EMODI)
%CHB -> Plutot revue Elsevier renewable energy
% Cth_filler=Cth_armure_matelas*0.5;                                              


%--------------------------------------------------------------------------
% Armour (steel or stainless steel)
%--------------------------------------------------------------------------
% Cp_th_acier=1.006*500;
% Rho_m_acier=7900;
Cth_Armure=3.8e6;   % volumetric specific heat(J/m3.K)
% p77 TableE2 IEC 60853-2

% Note: Metallic layers are neglected, as their thermal resistance is
% negligible compared with that of poly-composite materials.
% (p6 Nielsen, Applied Sciences 2019)

% Another note: the thermal resistance of the metallic parts in the cable,
% even though not equal to zero, is so small that it is usually neglected
% in rating computations (p34 Anders Rating of electric power cables 1997)

stainless_steel=1; % stainless steel is 1;  steel is 0

% the resistance of the armour at 20°C (?/m). p30 IEC 60287-1-1
if stainless_steel==1
    Rho_E_AG=70e-8;      % Stainless Steel
else
    Rho_E_AG=13.8e-8;  % Steel
end
% Note that Temperature coefficient of resistance is negligible for
% stainless steel; thus,no correction is necessary in Sub-clause 4.4.2.     
% (see p77 IEC853-2)

% temperature coefficient of steel at20°C (p30 Table1 IEC 60287-1-1)
if stainless_steel==1
    Alpha_E_AG=0;        % Stainless Steel
else
    Alpha_E_AG=4.5e-3;     % Steel
end
%--------------------------------------------------------------------------
% Outer serving (polypropylene yarn)
%--------------------------------------------------------------------------
Rho_th_Gaine_ext=5.5;    % thermal resistivity (K.m/W)
% Rho_th_Gaine_ext=10;   % thermal resistivity (K.m/W) Nielsen 2019

% p55 IEC 60287-2-1
Cth_Gaine_ext=1.8e6;     % volumetric specific heat(J/m3.K)
                         % p77 TableE2 IEC 60853-2

%--------------------------------------------------------------------------
% Soil
%--------------------------------------------------------------------------
% Kappa_Sol=1/0.7;
% Cp_sol=2e6;
Rho_sol=0.8;        % thermal resistivity (K.m/W) p73 Table D1 IEC 60853-2
Diff_sol=0.6e-6;    % diffusivity of soil m2/s  p73 Table D1 IEC 60853-2

De=Dgaine_ext_ext;  % external diameter of cable
Depth=1.75;         % depth of cable burying m
Depth=1.87;         % Value used in simulations, m

%% ************************************************************************
%% Calculation of equivalent thermal restistance (T1,T2,T3) for a model of
%% three-core cable
%% ************************************************************************
% 2 Calculation of thermal resistances
% 2.1 Thermal resistance of the constituent parts of a cable, T1, T2 and T3

% Important note  (p19 IEC 60287-2-1):
% Where screening layers are present, for thermal calculations metallic
% tapes are considered to be part of the conductor or sheath while
% semi-conducting layers (including metallized carbon paper tapes) are
% considered as part of the insulation. The appropriate component
% dimensions must be modified accordingly.

% 4.2.1.2 Three-core cables (% p25 IEC 60853-2)
% The three-core cable is replaced by an equivalent single-core
% construction % dissipating the same total conductor losses.
% The equivalent single-core conductor has a diameter: dc but firstly
% p22 IEC 60287-2-1
%--------------------------------------------------------------------------
% T1,the thermal resistance between one conductor and the sheath:

%  2.1.1.5 SL and SA type cables
% The thermal resistance T1 is calculated in the same way as for
% single-core cables (see next 2.1.1.1).

% 2.1.1.1 Single-core cables (p19 IEC 60287-2-1)
% The thermal resistance between one conductor and the sheath T1 is given
% by:
T1 = (Rho_th_isol/(2*pi))*log(1+2*((D_Matelas_ext/2)-rame)/(2*rame));
%--------------------------------------------------------------------------
% T2, The thermal resistance of fillers and bedding under the armour:

% 2.1.2.2 SL and SA type cables (p.27 IEC 60287-2-1)

% Firstly, it is neccesary to calculate G the geometric factor given in
% figure 6 (p71 IEC 60287-2-1)

% Thickness between sheaths and armour
Ep_Gaine_Armure = Thick_Gaine_PH+Thick_Armure_Matelas;                         
% Ep_Gaine_Armure = Thick_Armure_Matelas;

% dans Dorison, au début du livre (p36, tome 2)                                facteur X p36 Dorison : contradiction avec la norme (p70 IEC 60287-2-1)
% Verification avec le dessin 6 (p70 IEC 60287-2-1)
% Thickness of material between sheaths and armour expressed as a fraction
% of the outer diamter of the sheath
% facteur_X = Ep_Gaine_Armure/D_Pb_Screen_ext;                                

facteur_X = Ep_Gaine_Armure/D_Armure_Matelas_int;                             
% facteur_X = Ep_Gaine_Armure/D_Gaine_PH_ext;                                 

if (0<facteur_X && facteur_X<=0.03)
    facteur_G = 2*pi*(0.00022619+2.11429*facteur_X-20.4762*facteur_X^2);      % upper curve p51 CEI 60287-2-1
    %     facteur_G = 2*pi*(0.000202380+2.03214*facteur_X-21.6667*facteur_X^2);   % Nielsen 2019 
    
elseif (0.03<facteur_X && facteur_X<=0.15)
    facteur_G = 2*pi*(0.0142108+1.17533*facteur_X-4.49737*facteur_X^2+...     % upper curvep 51 CEI 60287-2-1
        10.6352*facteur_X^3);
    %     facteur_G = 2*pi*(0.0126529+1.101*facteur_X-4.56104*facteur_X^2+...     % Nielsen 2019
    %         11.5093*facteur_X^3);
end

% Finally, the thermal resistance of fillers and bedding under the armour:
T2 = (Rho_th_gaine_PH/(6*pi))*facteur_G;                                      
% T2 = (Rho_th_armure_matelas/(6*pi))*facteur_G;                               


%--------------------------------------------------------------------------
% T3,the thermal resistance of outer covering (serving):

% 2.1.3 Thermal resistance of outer covering (serving) T3(p27 IEC60287-2-1)
% The external servings are generally in the form of concentric layers and
% the thermal resistance T3 is given by:
T3 = (Rho_th_Gaine_ext/(2*pi))*log(1+2*(Ep_Gaine_Ext)/(Dext_arm));
%--------------------------------------------------------------------------
% T4, the external thermal resistance:

% 2.2 External thermal resistance T4 (p29 IEC 60287-2-1)
% 2.2.2 Single isolated buried cable T4(p31 IEC 60287-2-1)
u = 2*Depth/De;
T4 = (Rho_sol/(2*pi))*log(u+sqrt(u^2-1));

%% ************************************************************************
%% Calculation of the equivalent thermal capacitances for a model of
%% three-core cable.
%% ************************************************************************
% 4.2.1.2 Three-core cables (p25 IEC 60853-2):
% The three-core cable is replaced by an equivalent single-core
% construction dissipating the same total conductor losses:

% external diameter of dielectricc. The same value of diameter over
% dielectric (under the sheath) as for the three-core cable
Di=D_Matelas_ext;
% Important note (p25 IEC 60853-2):
% Where screening layers are present: for thermal calculations metallic
% tapes are considered to be part of the conductor or sheath while
% semi-conducting layers (including metallized carbon paper tapes) are
% considered as part of the insulation. The appropriate component
% dimensions should be modified accordingly.

% T1 is one-third of the value for one of the cores of the three-core cable
% as given in IEC Publication 287,

% The equivalent single-core conductor has a diameter (dc):
dc=(Di*exp(-(2*pi*T1/3)/Rho_th_isol));

pVW = (1/(2*log(Di/dc)))-(1/((Di/dc)^2-1));

% Alternatives:
% % (Dorison, dès % la page 136 (pour câble tri-écranté et p156).
% % Coefficient de répartition p (long) et p* (court)
%
% % p156: "r1 et r2 sont les rayons de l'âme et sur isolant" ->
% % interprétation T1: rayon intérieur du diélectrique et T2 rayon extérieur
% % du diélectrique
% %
% r2_sur_r1=(D_CU_Screen_int/2)/rame;% internal diameter écran cuivre
% %(contenant âme+ 2 écrans semi-conducteurs + diélectrique)/rayon de l'âme
%
% %pVW=1/((r2_sur_r1)^2-1)+1/(2*log(r2_sur_r1));  % p
%
% % p (Dorison, p156) (corrigé pendant la révision)
% pVW=-1/((r2_sur_r1)^2-1)+1/(2*log(r2_sur_r1));                              
%
% % Alternative calculation of pVW (p10 Nielsen, Applied Sciences 2019)
% Dc_star=D_Isol_ext*e((-2*pi*T1)/(Rho_th_isol*3));
% pVW=1/(2*log(D_Isol_ext/Dc_star))-1/((D_Isol_ext/Dc_star)^2-1);
%--------------------------------------------------------------------------

% (p25 IEC 60853-2):
% Thermal capacitances are calculated on the following assumptions:
% a) The actual conductors are considered to be completely inside the
% diameter of the equivalent single-core conductor, the remainder of the
% equivalent conductor being occupied by insulation.
% (p27 IEC 60853-2):
% b) The space between the equivalent single-core conductor and the sheath
% is considered to be completely occupied by insulation

% (p27 IEC 60853-2):
% The factor p is then calculated using the dimensions of the equivalent
% single-core cable and is applied to the thermal capacitance of the
% insulation based on assumption b) above.

% thermal capacitance of conductor (or equivalent single-core conductor of
%  a three-core cable, see Sub-clause 4.2.1.2)(p25 IEC 60853-2):
Qc = pi*rame^2*Cth_Al;
Qc = S_ame_pertes*Cth_Al;
% Qc = pi*(dc/2)^2*Cth_Al;
% Qc is for equivalent single-core cable ie dc is used instead of rame

% total thermal capacitance of dielectric per conductor (or equivalent
% single-core conductor of a three-core cable, see Sub-clause 4.2.1.2)
Qi = Cth_Isol*(pi/4)*(D_Pb_Screen_int^2-(2*rame)^2);

% Total thermal capacitance of the screen
Qs=Cth_ecran_Pb*(pi/4)*(D_Pb_Screen_ext^2-D_Pb_Screen_int^2);

% Total thermal capacitance of the filler
Qf = Cth_filler*((pi/4*(D_Armure_Matelas_ext^2))-...                              
    3*(pi/4*(D_Gaine_PH_ext^2)));

% Total thermal capacitance of the armour
Qa  = Cth_Armure*Nfil_Arm*pi*Ep_arm^2/4;

% Total thermal capacitance of the outer serving
Qj = Cth_Gaine_ext*(pi*(Dgaine_ext_ext^2-Dgaine_ext_int^2)/(4));

% Note that IEC does not give precise instructions on how to calculate
% thermal capacitances of the equivalent single-core conductor.
% Nevertheless,G. Anders writes in "Transient analysis of 3-core SL-type
% submarine cables with jacket around each core", Jicable 2015:
% The thermal capacitances of the (equivalent single-core)  conductor, the
% insulation and the screen are equal to three times those of the SL cable
% Therefore we redefine:
Qc=3*Qc;
Qi=3*Qi;
Qs=3*Qs;

%% ************************************************************************
%% Loss factors Lambda_1 and Lambda_2
%% ************************************************************************
% 2 Calculation of losses
% 2.1 AC resistance of conductor (p13 IEC 60287-1-1)
% 2.1.1 DC resistance of conductor at 20 degC
% Note 287-1-1:The value of R_DC20 shall be derived directly from IEC 60228.
% Where the conductor size is outside the range covered by IEC 60228,
% the value of R_DC20 may be chosen by agreement between manufacturer and
% purchaser:

% R_DC20=Rho_0_E_Al*Lcable/S_ame_pertes;

% R_DC20 for 1000 mm2 (p9 IEC 60228)?/m!
R_DC20=0.0291e-3;
% Note that IEC 60228 provides R_DC20 in ?/km and IEC 60287-1-1 uses ?/m !

% The d.c. resistance per unit length of the conductor (?/m)at its maximum
% operating temperature ? is given by (p14 IEC 60287-1-1):
R_DC90=R_DC20*(1+Alpha_E_Al*(90-20));

% Type of conductor: Aluminium. Round, stranded.(p31 Table 2 IEC 60287-1-1)
ks=1;

% Although there are no accepted experimental results dealing specifically
% with the coefficient kp for aluminium conductors, it is recommended that,
% for stranded aluminium conductors, the values given for similar copper
% conductors are used. (p31 Table 2 IEC 60287-1-1)
kp=1;


% 2.1.2 Skin effect factor Ys(p14 IEC 60287-1-1)
Xs=sqrt(8*pi*f/R_DC90*ks*1e-7);

% The above formula is accurate providing xs does not exceed 2,8, and
% therefore applies to the majority of practical cases.
if Xs>=2.8
    error('Xs>=2.8')
end

% The skin effect factor ys is given by (p14 IEC 60287-1-1):
Ys=Xs^4/(192+0.8*Xs^4);

% The Ys formula is accurate providing Xs does not exceed 2.8, and
% therefore applies to the majority of practical cases.

% 2.1.4 Proximity effect factor yp for three-core cables and for three
% single-core cables (p15 IEC 60287-1-1):
% 2.1.4.1 Circular conductor cables
Xp=sqrt(8*pi*f/R_DC90*kp*1e-7);
% The above formula is accurate providing xp does not exceed 2,8, and
% therefore applies to the majority of practical cases. (p15 IEC 60287-1-1)
if Xp>=2.8
    error('Xp>=2.8')
end

% Distance between conductor axes (m) (equal to the external diameter of
% internal sheaths):
s=D_Gaine_PH_ext;

% the distance between the axis of a conductor and the cable centre (m)
c=s/sqrt(3);

% Only if Xp<2.8 else Xp=>2.8?
% The proximity effect factor is given by:
Yp=(Xp^4/(192+0.8*Xp^4))*(D_ame/s)^2*(0.312*((D_ame/s)^2)...                
    +1.18/(Xp^4/(192+0.8*Xp^4)+0.27));
% Yp=(Xp^4/(192+0.8*Xp^4))*(dc/s)^2*(0.312*((dc/s)^2)...                    
%     +1.18/(Xp^4/(192+0.8*Xp^4)+0.27));
% disp('Yp formula is changed. dc is used instead of D_ame')

% 2.3.10 Cables with each core in a separate lead sheath (SL type)&armoured
% (p24 IEC 60287-1-1)
% For a three-core cable of which each core has a separate lead sheath ???
% is zero and the loss factor for the sheaths is given by:

% d, mean diameter of sheath or screen
d=D_moy_screen;                                                             
Xcable=2*w*1e-7*log(2*s/d);

% 2.1 AC resistance of conductor (p13 IEC 60287-1-1)
% The a.c. resistance per unit length of the conductor at its maximum
% operating temperature is given by the following formula
R_AC90_eval=R_DC90*(1+Ys+Yp);

% 2.2 Dielectric losses (applicable to a.c. cables only) p16 IEC 60287-1-1
% the capacitance per unit length (F/m)
% C=eps_r*1e-9/(18*log(D_Matelas_ext/(2*rame)));                               
% It seems the formula with D_Isol_ext/(D_Isol_int is correct:
C=eps_r*1e-9/(18*log(D_Isol_ext/(D_Isol_int)));                               

% The dielectric loss per unit length in each phase is given by p16 IEC
% 60287-1-1:
Wd=C*w*(U_0)^2*tan_delta;

% The cross section of screen
S_ecran=(pi/4)*(D_Pb_Screen_ext^2-D_Pb_Screen_int^2);% Check no data in IEC
% Recran=Rho_0_E_Pb*Lcable/S_ecran;
% Rs=Recran;   % ou Rgaine? the same for Xcable above

%  the resistance of the cable sheath or screen at 20 °C (?/m).
%  p17 IEC 60287-1-1 but formula is taken from Cablesizer
%  https://www.cableizer.com/documentation/R_sh/ :
Rso=Rho_0_E_Pb/S_ecran;

%  the resistance of the armour at 20°C:
%  p25 IEC 60287-1-1 but formula is assumed based on the resistance of the
%  cable sheath or screen at 20 °C (?/m).
Rao=Rho_E_AG/(Nfil_Arm*pi*(Dfil_Arm/2)^2);                                  % Hypothesis based on SEMREV 
%--------------------------------------------------------------------------
% While cycle for calculation of lambda1 and lambda2
iteration=0;     % iteration number (used for debugging only)
Error=1;         % any arbitrary number >=1e-3
Tmax=90+273.15;  % maximal operating temperature of conductor

% First guess of steady-state current (needed for calculation of sheath and
% armour temperature in the while cycle):
I_adm_NEW=I_admissible(Wd,T1,T2,T3,T4,3,R_AC90_eval,0.001,0.001,Tmax,Tamb);
I_adm_OLD=I_adm_NEW;

while (Error>=1e-3)
    % maximum operating temperature of sheath or screen is given by
    % p17 IEC 60287-1-1:
    Temp_sc = 90-(I_adm_NEW^2*R_AC90_eval+0.5*Wd)*T1;
    % The resistance of the sheath or screen at its maximum operating
    % temperature is given by p17 IEC 60287-1-1:
    Rs = Rso*(1+Alpha_E_Pb*(Temp_sc-20));

    % Ratio of the total losses in sheath respectively to the total
    % conductor losses (or losses in one sheath to the losses in one
    % conductor)(p24 IEC 60287-1-1)
    Lambda1=(Rs/R_AC90_eval)*(1.5/(1+(Rs/Xcable)^2));
    % NB: Here Lambda1 is actually Lambda1'. In general Lambda1 is
    % calculated as Lambda1=Lambda1'+Lambda1'' (see p16  IEC 60287-1-1)but
    % in accordance with p24  IEC 60287-1-1 Lambda1''==0. Thus, we assume
    % that Lambda1=Lambda1';
    
    % The maximum operating temperature of the armour is given by p25 IEC
    % 60287-1-1:
    Temp_AG = 90-((I_adm_NEW^2*R_AC90_eval+0.5*Wd)*T1+(I_adm_NEW^2*...
        R_AC90_eval*(1+Lambda1)+Wd)*3*T2);
    
    % The resistance of the armour at its maximum operating temperature is
    % given by p25 IEC 60287-1-1:
    Rarm=Rao*(1+Alpha_E_AG*(Temp_AG-20));
    
    % 2.4.2.5 SL type cables (p29 CEI 60287-1-1)
    % Where the armour is over a SL type cable, the screening effect of the
    % sheath currents reduces the armour loss. The formula for ?2 given in
    %  2.4.2.3.1 shall be multiplied by the factor
    Lambda1_dash=(Rs/R_AC90_eval)*1/(1+(Rs/Xcable)^2);
    
    Factor=(1-(R_AC90_eval/Rs)*Lambda1_dash);
    % where Lambda1_dash is obtained from 2.3.1.                            % 2.3.1 is for Two single-core cables, and three single-core cables
    
    % Ratio of the total losses in armour respectively to the total
    % conductor  losses (or losses in one sheath or armour to the losses
    % in one conductor) p28 CEI 60287-1-1
    
    Lambda2=1.23*(Rarm/R_AC90_eval)*(2*c/Dmoy_Arm)^2*1/(1+(2.77*Rarm...
        *1e6/(w))^2);
    Lambda2=Factor*Lambda2;
    if stainless_steel==1
        Lambda2=0; % if stainless steel is used
    end
    % Steady-state rating of cable
    I_adm_NEW=I_admissible(Wd,T1,T2,T3,T4,3,R_AC90_eval,Lambda1,Lambda2,...
        Tmax,Tamb);
    Error=abs(I_adm_NEW-I_adm_OLD);
    iteration=iteration+1;
    I_adm_OLD=I_adm_NEW;
end
if stainless_steel==1
    disp('Stainless steel is assumed for armour. Lambda2=0 ie no losses in armour')
end
% Final (eventual) steady-state current
I_adm=I_admissible(Wd,T1,T2,T3,T4,3,R_AC90_eval,Lambda1,Lambda2,Tmax,Tamb);

%% ************************************************************************
%% Equivalent thermal resistances and capacitances
%% ************************************************************************
% Appendix A provides a detailed method for reducing a multiple thermal
% circuit to a two-cell circuit.

% 4.2.2.2 Representation of common types of cable (p27. IEC 60853-2)
% (f) Armoured cables with each core in a separate lead sheath (SL type)
% (p31 IEC 60853-2)

% Note T1  is the thermal resistance of the equivalent single-core cable.
% T1, is one-third of the value for one of the cores of the three-core
% cable as given in IEC Publication 287, i.e TA=T1/3

TA=T1/3;

% qs ratio: losses in conductors and screens/losses in conductors
qs=Lambda1+1;
% Not specified in IEC.
% Source: (https://www.cableizer.com/documentation/q_s/)

% Tf,thermal resistance of filling i.e.thermal resistance of filling
% between cores of SL type cable
Tf=T2; % no reference on how Tf is calculated thus assumed equal to T2
Tf=0.020344899802; 

% qa, ratio: losses in (conductors + sheaths + armour)/losses in conductors
qa=Lambda1+Lambda2+1;
% Not specified in IEC standard.
% Source qa=Lambda1+Lambda2+1(https://www.cableizer.com/documentation/q_a/)
% Another source confirming it is ANDERS, George J, George, GEORGALLIS.
% "Transient analysis of 3-core SL-type submarine cables with jacket around
% each core" 2015

TB=qs*Tf+qa*T3;    % intiial IEC formula
% TB=3*(qs*Tf+qa*T3);  % Anders modification in "Transient analysis of
% % 3-core SL-type submarine cables with jacket around
% % each core" 2015
% disp('Anders modification of TB is used')


QA=Qc+pVW*Qi;        % Attention Qc and Qi of equivalent conductor or not
% QA=3*(Qc+pVW*Qi);    % Attention Qc and Qi of equivalent conductor or not
% disp('QA is multiplied by 3. Hypothesis')

Ts=0;

QB=(1-pVW)*Qi+((Qs+0.5*Qf)/qs)+((qa*T3/(qs*T2+qa*T3))^2)*((0.5*Qf/qs)+...
    ((Qa+Qj)/qa));  % intiial IEC formula

% QB=(1-pVW)*Qi+((Qs+Qf/6)/qs)+((3*qa*T3/(qs*T2+3*qa*T3))^2)*(((Qf/6)/qs)+...
%     ((Qa+Qj)/(3*qa))); % Anders modification in "Transient analysis of
% %                        % 3-core SL-type submarine cables with jacket around
% %                         % each core" 2015
% disp('Anders modification of QB is used')
%% Partials transients calculation
% 4.2.3 Calculation of cable partial transient (p33 IEC 853-2)
% The transient response of a cable circuit to a step function of load
% current, considered in isolation, that is with the right-hand pair of
% terminals in Figure 2 short-circuited, is found as follows:

M0=0.5*(QA*(TA+TB)+QB*TB);
N0=QA*TA*QB*TB;

a=(M0+sqrt(M0^2-N0))/N0;

b=(M0-sqrt(M0^2-N0))/N0;

Ta=(1/(a-b))*((1/QA)-b*(TA+TB));

Tb=TA+TB-Ta;
%--------------------------------------------------------------------------

% Ensuring the column vector
if (size(Iech_load,2) > size(Iech_load,1)) % If number of columns > number of rows
    Iech_load = Iech_load';
end

if isempty(optimization)
    % extract  all values of load
    I_ech=Iech_load;
elseif optimization==1  %if optimization is in progress
    % extract last 96 values of load for J day
    I_ech=Iech_load(days*96-95:end,1);
else % something wrong
    error('Check the optimization status');
end

if memoring==1
    % extract J day and the rest (zero arrays)
    I_ech=Iech_load(days*96-95:end,1);
end

% Duration of one current value in seconds
duration=900;

%-------------------------Initial code-------------------------------
tstart=900;    % choosing tstart=3600 instead of 0.001 allows avoiding
% the problem with temperature_correction_on_losses.m later
% timestep=5;   % time step, 5 seconds
timestep=900;  % time step, 1 hour = 3600 s
tend=duration;
t=tstart:timestep:tend*length(I_ech);     % time IEC in seconds
t=t';
Decalage_down=zeros(duration/timestep,1); % time offset for down current
%--------------------------------------------------------------------------
if preloading==1
    I_ech=I_ech(1:35040);
end

%--------------------------------------------------------------------------
% Superposition principle is applied in this "for cycle":
for k=1:length(I_ech) % for each current value
    if memoring==1 % if memoring status of IEC60853_2 is switched on
        if k<=(days*96+1)-(days*96-95) % continue k until 96 values are passed
            % Define the time offset for starting point
            Decalage=zeros(duration/timestep*k-duration/timestep,1);
            
            % Define I as k-th current value from I_ech
            I=I_ech(k);
            
            
            % Find cable and environment trainsient as well as attainement factor
            [Theta_ct,Theta_Et,Alpha_t]=losses_partial_transients(I);
            
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
                theta_VECTOR=Theta+Theta_down;
                
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
                theta_VECTOR=theta_VECTOR+Theta+Theta_down;
                
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
        [Theta_ct,Theta_Et,Alpha_t]=losses_partial_transients(I);
        
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
            theta_VECTOR=Theta+Theta_down;
            
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
            theta_VECTOR=theta_VECTOR+Theta+Theta_down;
            
        end % end of if k == 1
    end % end of if memoring==1
end % go to next k or end of for cycle

%% Thermal memory consideration in the optimization problem
% if optimization==1
%     if isempty(thermal_memory) % no thermal memory
%         %do nothing
%     else % there is a thermal memory
%         theta_VECTOR=theta_VECTOR+thermal_memory(days*96-95:days*96,1); % J day
%     end
% end

%% Correction of R as a function of T
Error=1e-3;
MAE=100;

if preloading==1
    I_ech=Iech_load;
end
% The conductor steady-state temperature rise above ambient:
% Source: 1.4.1 Buried cables where drying out of the soil does not occur
% % or cables in air (1.4.1.1 AC cables) (p10 IEC 60287-1-1)
Theta_inf_dielectric=(0.^2.*R_AC90_eval+0.5*Wd)*T1+(0.^2.*R_AC90_eval*(1+Lambda1)+...
    Wd)*3*T2+(0.^2.*R_AC90_eval*(1+Lambda1+Lambda2)+Wd)*3*(T3+T4);

while (MAE>=Error) % Mean absolute error in degrees of C
    if preloading==1
        I_ech=Iech_load;
    end
    [theta_VECTOR,MAE]=temperature_correction_on_losses_new(theta_VECTOR,f,ks,kp,D_ame,s,I_ech,duration,...
        timestep,Decalage_down,Rso,Rao,Xcable,c,Dmoy_Arm,w,stainless_steel);
end

theta_VECTOR_ISGT=theta_VECTOR;

if isempty(optimization)&& isempty(memoring)&& ~isempty(thermal_memory_preload)
    theta_VECTOR=theta_VECTOR+thermal_memory_preload; % thermal_memory_preload should be at studied horizon!
end

if optimization==1
    if isempty(thermal_memory) % no thermal memory
        %do nothing
    else % there is a thermal memory
        theta_VECTOR=theta_VECTOR+thermal_memory(days*96-95:days*96,1); % J day
    end
end

Theta_at=theta_VECTOR;

if memoring==1 % if memoring is activated
    theta_VECTOR_global=theta_VECTOR; % prepare the transient for thermal memory
    if isempty(thermal_memory) % if thermal_memory does not exist
        thermal_memory=theta_VECTOR_global; % save transient as thermal_memory
    else % there is already a thermal memory
        thermal_memory(days*96-95:end,1)=theta_VECTOR_global+thermal_memory(days*96-95:end,1); % save transients for period [studied day : end of horizon]
    end
    Theta_at=thermal_memory(days*96-95:end,1);
    %     if isempty(TEMPERATURE_D_1)
    %         TEMPERATURE_D_1{1,1}=thermal_memory;
    %     else
    %         TEMPERATURE_D_1{end+1,1}=thermal_memory;
    %     end
end

%% Memoring
% if memoring==1 % if memoring is activated
%     theta_VECTOR_global=theta_VECTOR; % prepare the transient for thermal memory
%     if isempty(thermal_memory) % if thermal_memory does not exist
%         thermal_memory=theta_VECTOR_global; % save transient as thermal_memory
%     else % there is already a thermal memory
%         thermal_memory(days*96-95:end,1)=theta_VECTOR_global+thermal_memory(days*96-95:end,1); % save transients for period [studied day : end of horizon]
%     end
% end % end of if memoring== 1


% if isempty(optimization)&& isempty(memoring)
%     MAE=100;
%     while (MAE>=Error) % Mean absolute error in degrees of C
%         [theta_VECTOR,MAE]=temperature_correction_on_losses_new(theta_VECTOR,f,ks,kp,D_ame,s,I_ech,duration,timestep,Decalage_down);
%     end
%     if ~isempty(thermal_memory_preload)
%         theta_VECTOR=theta_VECTOR+thermal_memory_preload(377*96+1:377*96+days*96,1);
%     end
%     Theta_at=theta_VECTOR;
% elseif optimization==1 % if optimization is activated
%     theta_VECTOR1=thermal_memory;
%     theta_VECTOR1(days*96-95:days*96,1)=theta_VECTOR; % theta_VECTOR considers thermal_memory
%     MAE=100;
%     while (MAE>=Error) % Mean absolute error in degrees of C
%         [theta_VECTOR,MAE]=temperature_correction_on_losses_new(theta_VECTOR1(days*96-95:days*96,1),f,ks,kp,D_ame,s,I_ech,duration,timestep,Decalage_down);
%         disp('theta_VECTOR1(days*96-95:days*96,1)???')
%         theta_VECTOR1(days*96-95:days*96,1)=theta_VECTOR; % theta_VECTOR after correction R(T)
%         Theta_at=theta_VECTOR1(1:days*96,1);
%     end
% elseif memoring==1 %if memoring is activated
%     theta_VECTOR1=thermal_memory(days*96-95:end,1);
%     MAE=100;
%     while (MAE>=Error) % Mean absolute error in degrees of C
%         [theta_VECTOR1,MAE]=temperature_correction_on_losses_new(theta_VECTOR1,f,ks,kp,D_ame,s,I_ech,duration,timestep,Decalage_down);
%         Theta_at=theta_VECTOR1;
%     end
% end
%% Preparing the temperature profile
% Initial temperature in Kelvin (K) considering dielectric losses
Tinit=Tamb+Theta_inf_dielectric;
% Where Tinit was defined earlier as Tinit=Tamb (see beginning of  script)

if isempty(thermal_memory)&& isempty(optimization)&& isempty(memoring)&& isempty(thermal_memory_preload) %if thermal_memory
    % Temperature in degC at the beginning of transient t=0
    Temperature=Tinit-273.15;
elseif isempty(optimization)&& isempty(memoring)&& ~isempty(thermal_memory_preload)
    Temperature=Tinit+thermal_memory_preload(1)-273.15;
elseif optimization==1 %& ~isempty(thermal_memory)
    if days>1
        Temperature=Tinit+thermal_memory(days*96-96)-273.15;
    else % days==1 or other
        Temperature=Tinit-273.15;
    end
elseif memoring==1 %& ~isempty(thermal_memory)
    if days>1
        Temperature=Tinit+thermal_memory(days*96-96)-273.15;
    else % days==1 or other
        Temperature=Tinit-273.15;
    end
elseif isempty(optimization) && isempty(memoring)
    Temperature=Tinit-273.15;
else% otherwise
    if days>1
        Temperature=Tinit+thermal_memory(days*96-96)-273.15;
    else % days==1 or other
        Temperature=Tinit-273.15;
    end
end

% Temperature profile in degC (final result)t>0
% Temperature(end+1:length(Theta_at)+1,1)=Theta_at+Tamb+Theta_inf_dielectric-273.15;
Temperature(end+1:length(Theta_at)+1,1)=Theta_at+Tamb+Theta_inf_dielectric-273.15;
disp('Temperature is changed !!!')
disp(['MAE ',num2str(MAE)])
% Note that if Iech_load had a size [24x1] then Temperature will have
% the size [25x1]. The fisrt value of Temperature  corresponds to
% temperature value at t=0 and the second one to t=3600s whereas the first
% value of Iech_load corresponds to t=3600s

% The maximal temperature of submarine cable
T_max=max(Temperature);
% if optimization==1
%     if days==4
%         if isempty(TEMPERATURE_D_1)
%             TEMPERATURE_D_1{1,1}=Temperature;
%             I_D_1{1,1}=I_ech;
%         else
%             TEMPERATURE_D_1{end+1,1}=Temperature;
%             I_D_1{end+1,1}=I_ech;
%         end
%     end
% end
end % end of function