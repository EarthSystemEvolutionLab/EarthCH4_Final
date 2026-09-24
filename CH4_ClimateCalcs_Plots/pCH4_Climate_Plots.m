%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climate and Atmospheric Chemistry Plots for Revised Phanerozoic pCH4
% 
% Receives and plots pCH4 and pO3 data from Phanerozoic *photochem* model
% solutions and computes paleoclimatic radiative forcings and GMST for
% comparison against the independent Judd et al. (2024) reconstruction.
% Wetland CH4 emissions are computed based on Beerling et al. (2009) in a
% separate script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial Constants

% scenario codes
Constant.invalT = 1; % NOMINAL CASE
Constant.invalCO2 = 1;
Constant.invalO2 = 1;

% plot colors
%colorCO2 = [120,94,240]./255;
%colorCH4 = [254,97,0]./255;
%colorsolar = [255,176,0]./255;

colorN2O = [220,38,127]./255;
colorN2Orange = 1.1.*colorN2O;
colorCO2 = [120,94,240]./255;
colorCH4 = [254,97,0]./255;
colorsolar = [255,176,0]./255;
colorRFtot = [0.8, 0, 0];%[0.5, 0.1, 0.2];
colorRFtotrange = colorRFtot; % 1.3*
colorT = [0 0 0];
darkgreen = [0,130,0]./255;
controlColor = [45,0,0]./255;
Emcolor = [255,50,70]./255;
KPgcolor = [255, 190, 0]./255;

% Thermodynamic constants
Constant.R = 8.3144; % J./mol*K; ideal gas constant = n*kB
Constant.N_avogadro = 6.022e23; % Avogadro's number (molecules./mol)

% From EONS
v.S_Pref         = 1361;             % W/m2; present day solar constant (Fs)
v.const.bol      = 5.67e-8;          % W/m^2K^4; Stefan-Boltzmann constant
v.ea.alb         = 0.29; % 0.3;      % earth albedo (~0.3 to 0.22) - Judd + 2024 value of 0.29 used here

T_ref = 287.15; % 14 C preindustrial modern reference temperature (all GHGs at PIM reference values, modern S(t))
rf.Albedo = v.ea.alb; %PGC.Albedo; % v.ea.alb; % can be variable?
% -0.0000618557.*(4500-time./1e6) + 0.29; % 0.3;              % earth albedo (~0.3 to 0.22)

% spatiotemporal scaling factors!
s_yr = (60.*60.*24.*365.25); % seconds per year
SA_Earth = 5.101e14; % m^2; earth surface area, from EONS (Horne and Goldblatt 2024), essentially invariant over time (cf. Scotese et al. 2021a)
photochem_fluxscaling = s_yr.*1e4.*SA_Earth./Constant.N_avogadro; % s*cm^2*mol CH4/(yr*global area*molecules CH4); converts photochem flux units (molec/cm^2/s) to mol/yr

%% Load in and extract *photochem* output data for pCH4 and other photochemical state variables over deep time

% Primary input structure
In_data = load('biome_output.mat');

% OLD Beerling et al. (2009) pCH4 reconstruction for comparison:
TimePhanCH4 = In_data.PhanBiomes.CH4PhanTime; % years BP
PhanCH4vals = In_data.PhanBiomes.CH4Phanppb; % in ppb: Beerling et al. 2009 pCH4 estimate 
% cf. Fig 6b, corrected from ppb to mixing ratio

% pCO2

TimePhanCO2 = 1e6.*(In_data.PhanBiomes.CO2PhanTime); % converted from Ma to yrs BP
TimePhanCO2 = cat(1,0,rmmissing(TimePhanCO2));
%if Constant.invalCO2 == 1 % nominal (median, 50th percentile)
PhanCO2vals = 1e-6.*In_data.PhanBiomes.CO2Phanppm; % transformation from ppm to mixing ratio abundance accounted for
%elseif Constant.invalCO2 == 2 % low (16th percentile, -sigma)
PhanCO2valsMin = 1e-6.*In_data.PhanBiomes.CO2Phanppm_16perc; 
%elseif Constant.invalCO2 == 3 % high (84th percentile, +sigma)
PhanCO2valsMax = 1e-6.*In_data.PhanBiomes.CO2Phanppm_84perc; 
%end
PhanCO2vals1 = rmmissing(PhanCO2vals);
PhanCO2valsmin = rmmissing(PhanCO2valsMin);
PhanCO2valsmax = rmmissing(PhanCO2valsMax);

PhanCO2valsF = cat(1,278e-6,PhanCO2vals1); % preindustrial value of 278 ppm assumed per Byrne and Goldblatt 2014 GRL paper
PhanCO2valsMinF = cat(1,278e-6,PhanCO2valsmin);
PhanCO2valsMaxF = cat(1,278e-6,PhanCO2valsmax);

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% Load and transform original CH4 emission model data and forcings

%AnchorFluxes = load('HighRes_InputFluxesCH4/time_GMST_pO2_FCH4_photochemInputAnchors.mat');
%AnchorTime = AnchorFluxes.timeUVb.*1e9; % years BP
%AnchorGMST = AnchorFluxes.Anchor_GMST_C; % degC
%AnchorpO2 = AnchorFluxes.AnchorpO2; % PAL
%AnchorFluxes.FCH4_tot_lo/hi.*(0.4.*26.7e12) % rescaled to CH4 mol/yr   time_GMST_pO2_FCH4_photochemInputHiRes
HiResFluxes = load('time_GMST_pO2_FCH4_photochemInputHiRes_updated.mat'); % HighRes_InputFluxesCH4/
PleistoFluxes = load('time_GMST_pO2_FCH4_photochemInputPleisto.mat'); % HighRes_InputFluxesCH4/

%F_CH4_min_photochem = HiResFluxes.FCH4_tot_hires_BC89_min.*(0.4.*1e11); % converts to molec/cm^2/s
%F_CH4_std_photochem = HiResFluxes.FCH4_tot_hires_BC89_std.*(0.4.*1e11); % converts to molec/cm^2/s
%F_CH4_max_photochem = HiResFluxes.FCH4_tot_hires_BC89_max.*(0.4.*1e11); % converts to molec/cm^2/s
FluxCH4_mid = HiResFluxes.FCH4_tot_hires_BC89_std.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
timeFluxB09 = HiResFluxes.time_CH4FluxB09.*1e9; % rescaled to yrs from Ga BP
FluxCH4_B09 = HiResFluxes.FluxCH4_B09.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
FluxCH4_max = HiResFluxes.FCH4_tot_hires_BC89_max.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
FluxCH4_min = HiResFluxes.FCH4_tot_hires_BC89_min.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% Load and transform *photochem* outputs for pCH4 and other state variables

% HiRespCH4_MinCalc = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_min_gammaT.mat');
% PleistopCH4_MinCalc = load('CH4_O3_outputs_Pleisto5ka_minGammaT_CH4flux_PHOTOCHEMcalc_PItuned_highestPIflux715_soilCons.mat'); % CH4_O3_outputs_Pleisto_hiT_photochemCalc.mat   CH4_O3_outputs_Pleisto5ka_strongT_CH4flux_PHOTOCHEMcalc_PItuned_higherPIflux

% GMST, pO2, Time all are identical between high and low T-sensitivity scenarios
% Total_pCH4_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);
% Total_Tau_CH4_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2));
% Total_pO3_DU_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,3); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3));
% Total_pO3_tropo_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,4);
% Total_O2flux_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7);
% Total_SoilDep_min = HiRespCH4_MinCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8);

% pCH4_Pleisto_min = PleistopCH4_MinCalc.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,1);


% Low EMISSIONS
%AnchorpCH4 = load('CH4_O3_outputs_PhaneroAnchors.mat');
%HiRespCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_weak_gammaT.mat'); % _5pt5MaPI
%PleistopCH4 = load('CH4_O3_outputs_Pleisto5ka_weakGammaT_CH4flux_PHOTOCHEMcalc_PItuned.mat'); %  CH4_O3_outputs_Pleisto.mat   CH4_O3_outputs_Pleisto5ka_weakT_CH4flux_PHOTOCHEMcalc_PItuned_higherPIflux

TotalTime = HiResFluxes.timeUVb_hires.*1e9; %cat(1,AnchorFluxes.timeUVb.*1e9,HiResFluxes.timeUVb_hires.*1e9);
TotalGMST = HiResFluxes.GMST_C_hires; %cat(1,AnchorFluxes.AnchorGMST_C,HiResFluxes.GMST_C_hires);
TotalGMSTMax = HiResFluxes.GMST_C84;
TotalGMSTMin = HiResFluxes.GMST_C16;

%Total_pCH4 = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,1),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,1));
%Total_Tau_CH4 = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2));
%Total_pO3_DU = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,3); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3));
%Total_pO3_tropo = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,4); % cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,4),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,4));
%Total_O2flux = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7);
%Total_SoilDep = HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8);

%time_Pleisto = PleistoFluxes.timeUVb_Pleisto5ka.*1e9;
%pCH4_Pleisto = PleistopCH4.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,1); % for Pleistocene only

% NOMINAL EMISSIONS
HiRespCH4_HiCalc = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_STDREF.mat'); % _5pt5MaPI
PleistopCH4_HiCalc = load('CH4_O3_outputs_Pleisto5ka_nominal_CH4flux_PHOTOCHEMcalc_PItuned.mat'); % CH4_O3_outputs_Pleisto_hiT_photochemCalc.mat   CH4_O3_outputs_Pleisto5ka_strongT_CH4flux_PHOTOCHEMcalc_PItuned_higherPIflux

% GMST, pO2, Time all are identical between high and low T-sensitivity scenarios
Total_pCH4_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);
Total_Tau_CH4_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2));
Total_pO3_DU_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,3); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3));
Total_pO3_tropo_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,4);
Total_O2flux_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7);
Total_SoilDep_high = HiRespCH4_HiCalc.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8);

pCH4_Pleisto_high = PleistopCH4_HiCalc.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,1);
Tau_Pleisto_high = PleistopCH4_HiCalc.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,2);
pOH_Pleisto_high = PleistopCH4_HiCalc.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,3);
time_Pleisto = PleistoFluxes.timeUVb_Pleisto5ka.*1e9;
%pCH4_Pleisto = PleistopCH4_HiCalc.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound(:,1); % for Pleistocene only
% Pleisto_high has same times exactly

% REFERENCE pCH4 (constant PI Emissions over Phanerozoic, tests influence of GMST/pO2/TauSun assumptions)
% HiRespCH4_Ref = load('CH4_O3_outputs_PhaneroHiRes_constFCH4.mat'); % 
% Total_pCH4_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);
% Total_Tau_CH4_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,2));
% Total_pO3_DU_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,3); %cat(1,AnchorpCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3),HiRespCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound(:,3));
% Total_pO3_tropo_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,4);
% Total_O2flux_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7);
% Total_SoilDep_ref = HiRespCH4_Ref.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8);


MaxEnvelope_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxEnvelope.mat'); % 
MaxEnv_pCH4 = MaxEnvelope_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MinEnvelope_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minEnvelope.mat'); % 
MinEnv_pCH4 = MinEnvelope_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MaxCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxCoal.mat'); % 
MaxCoal_pCH4 = MaxCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MinCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minCoal.mat'); % 
MinCoal_pCH4 = MinCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MaxTs_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxT.mat'); % 
MaxT_pCH4 = MaxTs_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MinTs_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minT.mat'); % 
MinT_pCH4 = MinTs_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MaxEas_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxEa.mat'); % 
MaxEa_pCH4 = MaxEas_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MinEas_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minEa.mat'); % 
MinEa_pCH4 = MinEas_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MaxO2s_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxO2.mat'); % 
MaxO2_pCH4 = MaxO2s_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

MinO2s_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minO2.mat'); % 
MinO2_pCH4 = MinO2s_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

% New Sensitivity Tests for Fermentation-Limitation and no Gamma_coal
Ferms_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_Ferm.mat'); % 
Ferm_pCH4 = Ferms_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

maxFerms_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxFerm.mat'); % 
maxFerm_pCH4 = maxFerms_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

minFerms_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minFerm.mat'); % 
minFerm_pCH4 = minFerms_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

NoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_StdNoCoal.mat'); % 
NoCoal_pCH4 = NoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

maxNoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxEnvelopeNoCoal.mat'); % 
maxNoCoal_pCH4 = maxNoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

minNoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minEnvelopeNoCoal.mat'); % 
minNoCoal_pCH4 = minNoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

FermNoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_FermNoCoal.mat'); % 
FermNoCoal_pCH4 = FermNoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

maxFermNoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxFermNoCoal.mat'); % 
maxFermNoCoal_pCH4 = maxFermNoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

minFermNoCoals_pCH4 = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_minFermNoCoal.mat'); % 
minFermNoCoal_pCH4 = minFermNoCoals_pCH4.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1);

% timesteps from 375 to 0 Ma at 1 Ma resolution
MegaTime = flip(0:1e6:375e6).';

%pCH4_Beerling = interp1(TimePhanCH4,PhanCH4vals,MegaTime,'pchip'); % ALWAYS SAME


% CO2 range for testing
CO2_scenario.Min = interp1(TimePhanCO2,PhanCO2valsMinF,MegaTime,'pchip');
CO2_scenario.Mid = interp1(TimePhanCO2,PhanCO2valsF,MegaTime,'pchip');
CO2_scenario.Max = interp1(TimePhanCO2,PhanCO2valsMaxF,MegaTime,'pchip');

% CH4 range
CH4_scenario.Min = MinEnv_pCH4;
CH4_scenario.Mid = Total_pCH4_high;
CH4_scenario.Max = MaxEnv_pCH4;

% T Range
GMST.Min = TotalGMSTMin;
GMST.Mid = TotalGMST;
GMST.Max = TotalGMSTMax;


%% Compute Radiative Forcing (up to 100 ppm pCH4) and compare CO2-CH4-solar irradiation climate model against Judd et al. (2024) GMST reconstruction
% modified in places from EONS (Horne and Goldblatt 2024)
    
% Solar Luminosity
    
%timeSolar = 4.5e9 - (flip(0:1e5:4e8)).';
% cf. agreement with UV flux from Sun since beginning of Earth, with time-power law evolution of flux intensity/surface area based on Zahnle and Walker 1982 Fig 9 caption equation, 
rf.F_solar = 1./(1+(0.4.*(4.57e9 - ((4.57e9-MegaTime)))./4.57e9)); % cf. Gough 1981, Feulner 2012, with modern irradiance (W/m^2) of 1361 (irrelevant due to self-normalization);  
    
rf.sc  = rf.F_solar.*v.S_Pref; % W/m2; solar luminosity over time ; (1 - (0.38.*(nt - 4e9) ./4.55e9)).^(-1) 
    
rf.Fs = (rf.sc./4).*(1-rf.Albedo);  % W/m2; solar energy influx to Earth's surface, drives climate warming, accounts for potentially dynamic global albedo
    
rf.Fs_o = (v.S_Pref./4).*(1-v.ea.alb); % preindustrial modern solar energy influx with PIM albedo
    
rf.Delta_Fs = rf.Fs - rf.Fs_o; % differential solar radiative forcing (negative, reaches 0 at reference PIM)

% % % 

% Derived from data file (fig5b.yaml) shared by Eric Wolf on May 7, 2026
%
% # Time-integrated sensitivity analysis
% # Project: Wolf et al. 2018, Figure 5b data
% # Generated: Retreived May 7th, 2026
% # Model: a proto version of ExoCAM, using the n28archean ExoRT scheme
% # Contact: eric.wolf@colorado.edu
% 
% metadata:
%   description: "Time-integrated climate sensitivity, Fig. 5b, Wolf et al. 2018"
%   units:
%     x: "Surface Temperature (K)"
%     y: "Climate sensitivity K/(W/m^2)"
%   notes: |
%     S### notation corresponds to:
%       S100: 1 × solar constant
%       S75:  0.75 × solar constant
%       S85:  0.875 x solar constant
%       S110: 1.1 x solar constant
% 
%   S100:
%     x: [283.402, 285.571, 288.150, 291.460, 295.238, 299.164, 303.994, 312.421, 325.541, 334.836, 338.165, 341.396, 345.624, 351.221, 358.742, 369.322]
%     y: [0.561218, 0.580272, 0.717566, 0.867801, 0.854644, 0.831606, 1.08929, 1.95090, 2.33870, 0.460613, 0.319745, 0.247341, 0.268273, 0.249032, 0.382138, 0.630249]
%     n_points: 16

% % % 

% Wolf et al. 2018 T-dependence climates sensitivity at equilibrium (from
% their Fig. 5b) for modern solar luminosity (So = 1) and for So = 0.875
% (Proterozoic)
dLuminosity = [0.875; 1]; % relative solar luminosity range
dTstep = (284:1:350).'; % spans Phanerozoic GMST range from <11 C up to >75 C
%Wolf_875_ClimateS = readtable('All_CSV_files/Wolf+2018_pt875So_climateSensitivity_vs_surfaceT_WPD.csv');
%Wolf875_T = Wolf_875_ClimateS{:,1}; % K
%Wolf875_G = Wolf_875_ClimateS{:,2};
%Wolf875_Grefined = interp1(Wolf875_T,Wolf875_G,dTstep,"pchip"); % linear

%Wolf_1000_ClimateS = readtable('All_CSV_files/Wolf+2018_1So_climateSensitivity_vs_surfaceT_WPD.csv');
Wolf1_T = [283.402; 285.571; 288.150; 291.460; 295.238; 299.164; 303.994; 312.421; 325.541; 334.836; 338.165; 341.396; 345.624; 351.221; 358.742; 369.322];
%Wolf_1000_ClimateS{:,1}; % ;
Wolf1_G = [0.561218; 0.580272; 0.717566; 0.867801; 0.854644; 0.831606; 1.08929; 1.95090; 2.33870; 0.460613; 0.319745; 0.247341; 0.268273; 0.249032; 0.382138; 0.630249];
%Wolf_1000_ClimateS{:,2}; % K/(W/m2) = C/(W/m2)
Wolf1_Grefined = interp1(Wolf1_T,Wolf1_G,dTstep,"pchip"); % linear

%WolfGc_matrix = cat(2,Wolf875_Grefined,Wolf1_Grefined); % forms an array


% Interpolate from Byrne and Goldblatt (2014) Fig. 5 red (1-bar) CH4
% radiative forcing curve at pCH4 > 100 ppm (up to 1% pCH4), assuming no
% N2O overlap (pN2O marginal ~ 270 ppb assumed)
BG14_highpCH4_RF = readtable('All_CSV_files/BG14_Fig5_1bar_RF_CH4_Wm-2.csv');
BG14_pCH4 = BG14_highpCH4_RF{:,1};
BG14_RF_CH4 = BG14_highpCH4_RF{:,2};
% NOT USED in nominal results


hilo = {'Min','Mid','Max'}; % ,'Min' ,'Low'
for ihl = 1:length(hilo)
    HL = hilo{ihl};

    % Interpolate climate sensitivity over time from So and GMST
    %ClimateSensitivityK_withLuminosityEffects.(HL) = interpn(dTstep,dLuminosity,WolfGc_matrix,(GMST.(HL) + 273.15),rf.F_solar,"linear"); 
    % unsure what other interpolation would be better

    ClimateSensitivityKWolf.(HL) = interp1(dTstep,Wolf1_Grefined,(GMST.(HL) + 273.15),'pchip');  % ,"linear",'extrap'
    % linear interp of modern solar luminosity case,
    % assuming that the <4% difference in solar luminosity since the late
    % Devonian will have only a marginal effect on the climate sensitivity
    % (especially in the late Cretaceous, where the effect will be closer 1%
    % less than modern luminosity) - we thus assume that the climate
    % sensitivity is roughly applicable over the mid-late Phanerozoic, to avoid
    % additional uncertainty from interpolating (linearly?) between the
    % Proterozoic and modern luminosity cases from Wolf+2018.
    ClimateSensitivityKJudd = 1.9453; % corresponds to 8 degC per CO2 doubling per Judd+2024
    % 1.7; % 0.7742; % roughly should equal 3 deg/CO2 doubling per BG14 function (Table 2) %1.5; 
    %1.9    2; % 1.5; % 3; % 3; % climate sensitivity (deg K/(W/m^2))
    % constant for now to reproduce 1 deg C change in T per W/m^2 - Goldblatt, personal comm.
    % can make this a T-dependent power-law (or pCO2/RF dependent - cf. He et al. 2023)
    % could go as high as 8 per Judd et al. 2024

    
   % CH4_scenario.Low = Total_pCH4;
    % CH4_scenario.Min = Total_pCH4_min;

    mr.CO2 = CO2_scenario.(HL); %interp1(TimePhanCO2,PhanCO2valsF,MegaTime,'pchip'); % ALWAYS SAME
    mr.CH4 = CH4_scenario.(HL); % Total_pCH4;
    %PGC.CH4; % 7.15e-7; % PIM constant test-case 
    % Flux.pCH4_max; % Flux.pCH4_wetlands; % Flux.pCH4_marine; %
    mr.N2O = 270e-9; % ppb pN2O at PIM, assumedly invariant for now:  all values account for moles of molecules (not of N, C, or O atoms in molecules)
    
    % Calculate Radiative Forcings relative to Preindustrial Modern
    % uses Byrne and Goldblatt 2014 (GRL) radiative forcing functions accounting for
    % spectral overlap to estimate radiative forcings of CO2, CH4, and N2O
    % relative to their preindustrial magnitudes - complements above
    % predictions using widely used expression for relative (not absolute) RF
    NM1 = 2.5e-6; % critical transition pX for X = N2O, CH4, change in RF function mode
    
    % CO2
    rf.CO2_atm = mr.CO2; % concentration (prescribed) 
    % use this as equivalent to ppmv
    % mixing ratios 
    %rf.CO2_fpp = mr.CO2./(Flux.atmoGas_tot); % concentration as partial pressure normalized to total atmospheric pressure
    rf.CO2_pm = rf.CO2_atm; % rf.CO2_fpp;
    CO2_o = 278e-6; % preindustrial modern CO2, Byrne and Goldblatt 2014
    rf.relConc_CO2 = rf.CO2_pm./CO2_o;
    
    % CH4
    rf.CH4_atm = mr.CH4;
    %rf.CH4_fpp = mr.CH4./(Flux.atmoGas_tot);
    rf.CH4_pm = rf.CH4_atm; %rf.CH4_fpp;
    CH4_o = 715e-9;
    rf.relConc_CH4 = rf.CH4_pm./CH4_o;
    
    % N2O
    rf.N2O_atm = mr.N2O;
    %rf.N2O_fpp = mr.N2O./(Flux.atmoGas_tot);
    rf.N2O_pm = rf.N2O_atm; %rf.N2O_fpp;
    N2O_o = 270e-9;
    rf.relConc_N2O = rf.N2O_pm./N2O_o;
    
    % NH3 not present in BG14, so not directly comparable
    
    % Calculate overlap for N2O-CO2 RF (reduction (-) added to RF of each gas)
    % EQUALS 0 OVER ALL TIME! CAN IGNORE EFFECTIVELY though included below
    % for BG14 functions
    rf.dRF_CO2_N2O = -16.16.*exp(-0.036.*(log(rf.CO2_pm - CO2_o) - 0.0024).^2 - 0.05.*(log(rf.N2O_pm - N2O_o) + 6.5).^2);
    % Calculate overlap for N2O-CH4 RF (reduction (-) added to RF of each gas)
    rf.dRF_CH4_N2O = -24.*exp(-0.02.*(log(rf.CH4_pm - CH4_o) - 0.01).^2 - 0.044.*(log(rf.N2O_pm - N2O_o) + 7.73).^2);
    
    % calculate RF for each gas depending on level of abundance in atmosphere
    %if rf.CO2_pm < 200e-6
    rf.gCO2 = log(1 + 1.2.*(rf.CO2_pm.*1e6) + 0.005.*((rf.CO2_pm.*1e6).^2) + 1.4e-6.*((rf.CO2_pm.*1e6).^3)); % note that pCO2 is in ppmv here, not ppv as in BG14
    rf.gCO2_o = log(1 + 1.2.*(CO2_o.*1e6) + 0.005.*((CO2_o.*1e6).^2) + 1.4e-6.*((CO2_o.*1e6).^3)); % note that pCO2 is in ppmv here, not ppv as in BG14
    rf.RF_CO2_1 = 3.35.*(rf.gCO2 - rf.gCO2_o); % + rf.dRF_CO2_N2O; % using WMO 1999 radiative forcing function, Table 1 Byrne and Goldblatt 2014 - no overlap at such low concentrations of CO2 typically
        % this is intended only for use when BG14 boundaries are violated (<200 ppm), otherwise use BG14
        % unsure how to apply CO2-N2O overlap, so none used at this low CO2 level - BG14 summative overlap
        % function ignored bc = 0 always anyway
    %elseif rf.CO2_pm >= 200e-6
    rf.RF_CO2_2 = 5.32.*log(rf.relConc_CO2) + 0.39.*(log(rf.relConc_CO2)).^2 + rf.dRF_CO2_N2O;
    %end
    
    %if rf.CH4_pm < 0.1e-6
    rf.RF_CH4_1 = 0.036.*(sqrt(rf.CH4_pm.*1e9) - sqrt(CH4_o.*1e9)) - (0.47.*log(1 + 2.01e-5.*(rf.CH4_pm.*1e9.*N2O_o.*1e9).^(0.75) + 5.31e-15.*(rf.CH4_pm.*1e9).*(rf.CH4_pm.*1e9.*N2O_o.*1e9).^(1.52)) ...
            - 0.47.*log(1 + 2.01e-5.*(CH4_o.*1e9.*N2O_o.*1e9).^(0.75) + 5.31e-15.*(CH4_o.*1e9).*(CH4_o.*1e9.*N2O_o.*1e9).^(1.52))); % uses ppbv, IPCC 1990 fits with overlap, Table 1 BG14
    %elseif (rf.CH4_pm < 2.5e-6) % (rf.CH4_pm >= 0.1e-6) && 
    rf.RF_CH4_2 = 1173.*(sqrt(rf.CH4_pm) - sqrt(CH4_o)) - 71636.*(sqrt(rf.CH4_pm) - sqrt(CH4_o)).^2 + rf.dRF_CH4_N2O; % BG14 Table 2
    %elseif (rf.CH4_pm < 100e-6) %(rf.CH4_pm >= 2.5e-6) &&  up to 100 ppmv
    rf.RF_CH4_3 = 0.824 + 0.8.*log(rf.CH4_pm./NM1) + 0.2.*(log(rf.CH4_pm./NM1)).^2 + rf.dRF_CH4_N2O;
    %end
    %rf.RF_CH4_4 = 0.824 + 0.8.*log(100e-6./NM1) + 0.2.*(log(100e-6./NM1)).^2 + rf.dRF_CH4_N2O; % capped at highest resolvable value for BG14 polynomial equations
    
    %if rf.N2O_pm < 0.1e-6
    % rf.RF_N2O_1 = 0.12.*(sqrt(rf.N2O_pm.*1e9) - sqrt(N2O_o.*1e9)) - (0.47.*log(1 + 2.01e-5.*(CH4_o.*1e9.*rf.N2O_pm.*1e9).^(0.75) + 5.31e-15.*(CH4_o.*1e9).*(CH4_o.*1e9.*rf.N2O_pm.*1e9).^(1.52)) ...
    %         - 0.47.*log(1 + 2.01e-5.*(CH4_o.*1e9.*N2O_o.*1e9).^(0.75) + 5.31e-15.*(CH4_o.*1e9).*(CH4_o.*1e9.*N2O_o.*1e9).^(1.52))); % IPCC 1990 fits with overlap, Table 1 BG14
    % %elseif (rf.N2O_pm < 2.5e-6) %  (rf.N2O_pm >= 0.1e-6) &&
    % rf.RF_N2O_2 = 3899.*(sqrt(rf.N2O_pm) - sqrt(N2O_o)) + 38256.*(sqrt(rf.N2O_pm) - sqrt(N2O_o)).^2 + rf.dRF_CO2_N2O + rf.dRF_CH4_N2O; % BG14 Table 2
    % %elseif (rf.N2O_pm < 100e-6) % (rf.N2O_pm >= 2.5e-6) && 
    % rf.RF_N2O_3 = 4.182 + 3.*log(rf.N2O_pm./NM1) + 0.5469.*(log(rf.N2O_pm./NM1)).^2 + rf.dRF_CO2_N2O + rf.dRF_CH4_N2O;
    % %elseif (rf.N2O_pm >= 100e-6)
    % rf.RF_N2O_4 = 4.182 + 3.*log(100e-6./NM1) + 0.5469.*(log(100e-6./NM1)).^2 + rf.dRF_CO2_N2O + rf.dRF_CH4_N2O; 
    %     % capped at maximum resolved value if pN2O exceeds max domain value of 100 ppm 
    %     % Note that N2O has overlap interference from both CO2 and CH4
    % %end
    
    % SETS CO2 RF for each for-loop case
    rf.RF_CO2 = rf.RF_CO2_2;
    rf.RF_CO2(rf.CO2_pm < 200e-6) = rf.RF_CO2_1(rf.CO2_pm < 200e-6); % accounts for N2O overlap
    RF_CO2.(HL) = rf.RF_CO2;
    
    rf.RF_CH4 = rf.RF_CH4_2;
    rf.RF_CH4(rf.CH4_pm < 0.1e-6) = rf.RF_CH4_1(rf.CH4_pm < 0.1e-6);
    rf.RF_CH4(rf.CH4_pm >= 2.5e-6) = rf.RF_CH4_3(rf.CH4_pm >= 2.5e-6); % assumed that CH4 never exceeds 100 ppm! Reasonable, never above 40 ppm in Phanerozoic
    rf.RF_CH4(rf.CH4_pm >= 100e-6) = interp1(log10(BG14_pCH4),BG14_RF_CH4,log10(rf.CH4_pm(rf.CH4_pm >= 100e-6)),'linear'); %; % interpolated from Byrne and Goldblatt at >100 ppm pCH4
    
    % SETS CH4 RF for each for-loop case
    RF_CH4.(HL) = rf.RF_CH4;
    %if ihl == 1
    %    RF_CH4.High = rf.RF_CH4;
    %elseif ihl == 2
    %    RF_CH4.Low = rf.RF_CH4;
    %else
    %   
    %end

    % rf.RF_N2O = rf.RF_N2O_2;
    % rf.RF_N2O(rf.N2O_pm < 0.1e-6) = rf.RF_N2O_1(rf.N2O_pm < 0.1e-6);
    % rf.RF_N2O(rf.N2O_pm >= 2.5e-6) = rf.RF_N2O_3(rf.N2O_pm >= 2.5e-6);
    % rf.RF_N2O(rf.N2O_pm >= 100e-6) = rf.RF_N2O_4(rf.N2O_pm >= 100e-6); % accounts for CO2-CH4 overlap
    
    %rf.deltaT_CO2 = ClimateSensitivityK.*rf.RF_CO2;
    %rf.deltaT_CH4 = ClimateSensitivityK.*rf.RF_CH4;
    %rf.deltaT_N2O = ClimateSensitivityK.*rf.RF_N2O;
    % No NH3 temperature sensitivity resolved by BG14, so must be assessed by
    % comparison of RF calculated in EONS script above to RF of N2O/CH4
    
    
    % % CASE: where atmospheric pN2O is negligible and hence no overlap!
    % %if rf.CO2_pm < 200e-6
    % %rf.gCO2 = log(1 + 1.2.*(rf.CO2_pm.*1e6) + 0.005.*((rf.CO2_pm.*1e6).^2) + 1.4e-6.*((rf.CO2_pm.*1e6).^3)); % note that pCO2 is in ppmv here, not ppv as in BG14
    % %rf.gCO2_o = log(1 + 1.2.*(CO2_o.*1e6) + 0.005.*((CO2_o.*1e6).^2) + 1.4e-6.*((CO2_o.*1e6).^3)); % note that pCO2 is in ppmv here, not ppv as in BG14
    % rf.RF_CO2_1_noN2O = 3.35.*(rf.gCO2 - rf.gCO2_o); % NO N2O overlap + rf.dRF_CO2_N2O; % using WMO 1999 radiative forcing function, Table 1 Byrne and Goldblatt 2014 - no overlap at such low concentrations of CO2 typically
    %     % this is intended only for use when BG14 boundaries are violated (<200 ppm),
    %     % otherwise use BG14
    %     % unsure how to apply CO2-N2O overlap - BG14 summative overlap function assumed
    %     % (for now)
    % %elseif rf.CO2_pm >= 200e-6
    % rf.RF_CO2_2_noN2O = 5.32.*log(rf.relConc_CO2) + 0.39.*(log(rf.relConc_CO2)).^2; % NO N2O + rf.dRF_CO2_N2O;
    % %end
    % 
    % %if rf.CH4_pm < 0.1e-6
    % rf.RF_CH4_1_noN2O = 0.036.*(sqrt(rf.CH4_pm.*1e9) - sqrt(CH4_o.*1e9)); % NO N2O - (0.47.*log(1 + 2.01e-5.*(rf.CH4_pm.*1e9.*N2O_o.*1e9).^(0.75) + 5.31e-15.*(rf.CH4_pm.*1e9).*(rf.CH4_pm.*1e9.*N2O_o.*1e9).^(1.52)) ...
    %       %  - 0.47.*log(1 + 2.01e-5.*(CH4_o.*1e9.*N2O_o.*1e9).^(0.75) + 5.31e-15.*(CH4_o.*1e9).*(CH4_o.*1e9.*N2O_o.*1e9).^(1.52))); % uses ppmv, IPCC 1990 fits with overlap, Table 1 BG14
    % %elseif (rf.CH4_pm < 2.5e-6) % (rf.CH4_pm >= 0.1e-6) && 
    % rf.RF_CH4_2_noN2O = 1173.*(sqrt(rf.CH4_pm) - sqrt(CH4_o)) - 71636.*(sqrt(rf.CH4_pm) - sqrt(CH4_o)).^2; % No N2O + rf.dRF_CH4_N2O; % BG14 Table 2
    % %elseif (rf.CH4_pm < 100e-6) %(rf.CH4_pm >= 2.5e-6) &&  assume that CH4 never exceeds 100 ppm
    % rf.RF_CH4_3_noN2O = 0.824 + 0.8.*log(rf.CH4_pm./NM1) + 0.2.*(log(rf.CH4_pm./NM1)).^2; % No N2O + rf.dRF_CH4_N2O;
    % %end
    % rf.RF_CH4_4_noN2O = 0.824 + 0.8.*log(100e-6./NM1) + 0.2.*(log(100e-6./NM1)).^2; % No N2O + rf.dRF_CH4_N2O; % capped at highest resolvable value for BG14 polynomial equations
    % 
    % rf.RF_CO2_noN2O = rf.RF_CO2_2_noN2O;
    % rf.RF_CO2_noN2O(rf.CO2_pm < 200e-6) = rf.RF_CO2_1_noN2O(rf.CO2_pm < 200e-6);
    % 
    % rf.RF_CH4_noN2O = rf.RF_CH4_2_noN2O;
    % rf.RF_CH4_noN2O(rf.CH4_pm < 0.1e-6) = rf.RF_CH4_1_noN2O(rf.CH4_pm < 0.1e-6);
    % rf.RF_CH4_noN2O(rf.CH4_pm >= 2.5e-6) = rf.RF_CH4_3_noN2O(rf.CH4_pm >= 2.5e-6); % assumed that CH4 never exceeds 100 ppm!
    % % IF CH4 exceeds 100 ppm -
    % rf.RF_CH4_noN2O(rf.CH4_pm >= 100e-6) = interp1(log10(BG14_pCH4),BG14_RF_CH4,log10(rf.CH4_pm(rf.CH4_pm >= 100e-6)),'linear'); % linear interpolation in semilogx space should be reasonably accurate (rf.CH4_pm >= 100e-6);
    % % This resolves the case where pN2O is negligible, hence no overlap at all
    % % (endmember case, allows for somewhat enhanced CH4-CO2 greenhouse with less
    % % overlap)
    
    
    %rf.RF_tot_CO2_CH4 = rf.RF_CH4_noN2O + rf.RF_CO2_noN2O; % no N2O, no N2O-CO2/CH4 overlap
    %rf.RF_tot_all = rf.RF_CH4 + rf.RF_CO2 + rf.RF_N2O; % accounts for overlap by subtracting overlap diminution from both RFs in each pair
    
    
    % Calculate the equilibrium temperature from all combined GHGs excluding N2O




    % Calculate equilibrium T with RF (solar, CO2, CH4)

    % Judd+2024 Constant Climate Sensitivity
    Teq_CH4_CO2.(HL) = T_ref + ClimateSensitivityKJudd.*(RF_CH4.(HL) + RF_CO2.(HL) + rf.Delta_Fs); % no N2O overlap, N2O at 270 ppb reference level so no overlap

    Teq_CO2.(HL) = T_ref + ClimateSensitivityKJudd.*(RF_CO2.(HL) + rf.Delta_Fs); % _noN2O Teq estimate with only CO2 as GHG (no CH4, no N2O, no spectral overlap)

    % Wolf+2018 T-dependent Climate Sensitivity
    Teq_CH4_CO2Wolf.(HL) = T_ref + ClimateSensitivityKWolf.(HL).*(RF_CH4.(HL) + RF_CO2.(HL) + rf.Delta_Fs); % no N2O overlap, N2O at 270 ppb reference level so no overlap

    Teq_CO2Wolf.(HL) = T_ref + ClimateSensitivityKWolf.(HL).*(RF_CO2.(HL) + rf.Delta_Fs); % _noN2O Teq estimate with only CO2 as GHG (no CH4, no N2O, no spectral overlap)
    

end





%% Generate Relevant Plots
polyx = cat(1,MegaTime./1e6,flip(MegaTime./1e6)); % for polyshape ranges

polyK = [33.9, 33.9, 149.24, 149.24]; % constant, so can plot as rectangle
polyconst = [33.9, 33.9, 149.24, 149.24]; % constant, so can plot as rectangle
polyKPg = cat(1,0, 48, 48, 0).'; % constant, so can plot as rectangle
polyKPgrf = cat(1,-20, 1000, 1000, -20).'; % constant, so can plot as rectangle
polyKPglog = cat(1,0.00001, 1000, 1000, 0.00001).'; % constant, so can plot as rectangle

figure(201);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(201),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
nexttile

yyaxis left

size = 45;

scatter(-5,0.715,'o','MarkerFaceColor','k')
hold on
%scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
%hold on
%plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
%hold on
%plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
%hold on
%plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on

plot(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k')
hold on
plot(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','^','LineStyle','-','MarkerEdgeColor',colorCH4,'MarkerFaceColor',colorCH4,'LineWidth',2,'MarkerSize',3)
hold on

%plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
%plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on

xline([0],'-k')
hold on
scatter(-5,0.565,'o','MarkerFaceColor','k')
hold on
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[-4,65]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('CH_4 (ppm)')
%ylim([-20,320]) % 350
pbaspect([2 1 1])
fontsize(24,"points") % 14
L = legend('Preindustrial CH_4 (0.565-0.715 ppm)','Phanerozoic CH_4 per Beerling+(2009)','Revised Phanerozoic CH_4','FontSize',24); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
%,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
L.AutoUpdate = 'off';
title('CH_4 Mixing Ratio (ppm)');

box on

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(24,"points") % 14
%set(gca,'YTickLabel',[]);
yticks([])


%yyaxis right

%geotimescale_Mills_JFHmod;
%hold on
%PhanTransitions;
%set(gca,'YTickLabel',[]);
%yticks([]);

%set(gca,'XDir','reverse');
%set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14


%figure(301);clf

%plot(AnchorTime,AnchorGMST_C)
%hold on
%plot(timeslices,GMST_C_hires)


%figure(401);clf

%plot(AnchorTime,CoalDepoRate)
%hold on
%plot(timeslices,CoalDepoRate_hires)


figure(501);clf

tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
nexttile

% Ice Core data over past 800 ka from Loulergue+2008 (as ref'd and shown in Schilt+2010) no filtering applied,
% taken from Supplemental Material Table - all time in years BP from 1950
% (no correction applied, trivial), all mean pCH4 values in ppb
IceCoreCH4 = readtable("All_CSV_files/Loulergue+2008_pCH4_800ka_icecore_cfSchilt+2010.xlsx");
IceCoreTime_yrsBP1950 = IceCoreCH4{:,1};
IceCorepCH4_ppb = IceCoreCH4{:,2};

%size = 45;

plot(IceCoreTime_yrsBP1950./1e3,IceCorepCH4_ppb,'o','MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none','MarkerFaceColor','none')
hold on

%plot(time_Pleisto./1e3,pCH4_Pleisto_hi(1:5:length(timeslices2)).*1e9,'^','MarkerFaceColor','m') % 'none' for no marker
%hold on
%plot(time_Pleisto./1e3,pCH4_Pleisto_min.*1e9,'d','MarkerFaceColor','none','MarkerEdgeColor','k')
%hold on

%plot(time_Pleisto./1e3,pCH4_Pleisto.*1e9,'v','MarkerFaceColor','b','MarkerEdgeColor','b')
%hold on
plot(time_Pleisto./1e3,pCH4_Pleisto_high.*1e9,'^','MarkerFaceColor',colorCH4,'MarkerEdgeColor',colorCH4,'MarkerSize',10)
hold on
%plot(timeslices./1e3,Phan_pCH4_hires_noT.*1e9,'g','Marker','none')
%hold on
%yline([100],'--k')

%plot(timeslices2./1e3,Phan_pCH4_Pleisto_hi.*1e9,'-m','Marker','none') % 'none' for no marker
%hold on
%plot(timeslices2./1e3,Phan_pCH4_Pleisto_lo.*1e9,'-b','Marker','none')
%hold on
scatter(-5,565,'square','MarkerFaceColor','k')
hold on
scatter(-5,715,'square','MarkerFaceColor','k')
hold on
errorbar(-5,565,NaN,150,'Color','k','CapSize',1,'LineWidth',4)
hold on

xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,800]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (ka)'); ylabel('CH_4 (ppb)')
%title('Pleistocene yCH_4 Data vs. Model Output')
%ylim([-20,350])
pbaspect([3 1 1])
fontsize(24,"points") % 14
L = legend('Ice Core CH_4',...
    'Modeled CH_4','FontSize',28); % 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)', 14,'Modern pN_2O (337 ppb)'
%     %'Revised Phanerozoic yCH_4 (Weak \gamma_T, per Zhu+2014)',

L.AutoUpdate = 'off';

box on



figure(5011);clf

tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
nexttile

scatter(-5,7,'square','MarkerFaceColor','k')
hold on
scatter(-5,10,'square','MarkerFaceColor','k')
hold on
plot(time_Pleisto./1e3,Tau_Pleisto_high,'^','MarkerFaceColor','m','MarkerEdgeColor','m')
hold on


xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,800]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]

xlabel('Age Before Present (ka)'); ylabel('CH_4 Lifetime (yr)')

pbaspect([3 1 1])
fontsize(24,"points") % 14


L.AutoUpdate = 'off';

box on


%figure(55555); clf

%scatter(Tau_Pleisto_high,pOH_Pleisto_high);




%figure(502);clf

%tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
%nexttile

%scatter(IceCoreTime_yrsBP1950./1e3,IceCorepCH4_ppb,'+','MarkerEdgeColor','k')
%hold on

%xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
%set(gca,'XDir','reverse');
%set(gca,'xlim',[0,1000],'ylim',[300,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (ka)'); ylabel('yCH_4 (ppb)')
%ylim([-20,350])
%pbaspect([3 1 1])
%fontsize(24,"points") % 14
%L = legend('Preindustrial Modern yCH_4 (715 ppb, ~565 ppb pre-agrarian)','Ice Core yCH_4 (last 800 ka, Loulergue+2008)',...
%    'Revised Phanerozoic yCH_4 (Weak \gamma_T, per Zhu+2014)','Revised Phanerozoic yCH_4 (Strong \gamma_T, per Conrad 2023)','FontSize',18); % 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)', 14,'Modern pN_2O (337 ppb)'
%     %'Revised Phanerozoic yCH_4 (Strong \gamma_T)',

%L.AutoUpdate = 'off';

%box on











figure(3303);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(3303),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(3,1,"TileSpacing","compact","Padding","compact");


nexttile

yyaxis left

plot(polyshape(polyconst,polyKPgrf),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%plot(TotalTime./1e6,rf.Delta_Fs,'color',colorsolar)
%hold on
%scatter(TotalTime(1:16)./1e6,RF_Flux.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
%plot(TotalTime(17:end)./1e6,RF_Flux.Teq_o.t(17:end) - 273.15,'color',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CH4(1:16),size,'v','MarkerFaceColor',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CO2(1:16),size,'v','MarkerFaceColor',colorCO2)
%hold on
%
plot(MegaTime./1e6, rf.Delta_Fs,'color',colorsolar,'LineStyle','-','LineWidth',2,'Marker','none'); % 
hold on 
plot(MegaTime./1e6, RF_CO2.Mid,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
hold on
%plot(TotalTime./1e6, rf.RF_CH4,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
%hold on
%plot(MegaTime./1e6,RF_CH4.Min,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on
%plot(MegaTime./1e6,RF_CH4.Low,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b')
%hold on
plot(MegaTime./1e6,RF_CH4.Mid,'LineStyle','-','color',colorCH4,'LineWidth',2,'Marker','none'); %'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m')
hold on
%plot(MegaTime./1e6,(RF_CH4.Low + rf.RF_CO2 + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none')
%hold on
plot(MegaTime./1e6,(RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',3,'Marker','none')
hold on

%annotation('textbox',[.4 .55-0.022 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(-5,2.7e-7,'ok','MarkerFaceColor','k','MarkerSize',5)
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.RF_N2O,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.RF_N2O,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
annotation('textbox',[.245 .84 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .55-.03 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline(0,'k')
hold on
plot(-5,0,'o','MarkerFaceColor',colorRFtot,'MarkerSize',5)
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[-10,15]) %-9, 12 % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(16,"points") % 14
pbaspect([3 1 1])
%xlabel('Age Before Present (Ma)');
%xticklabels([]);
ylabel('Radiative Forcing (W/m^2 Relative to PIM)')
% 'CH_4 (Weak \gamma_T)', 
L = legend('','Solar','CO_2','CH_4',...
        'Total','FontSize',12);
% IN CAPTION, DISCUSS HOW CO2 and CH4 do not account for N2O RF but N2O RF (and total RF) do
%legend('Radiative Forcing from N_2O (low crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (high crustal E_a, high emission scenario)','Radiative Forcing from N_2O (high crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
title('Solar and GHG Radiative Forcings')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])

box on



nexttile

yyaxis left

%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 

L1 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(MegaTime./1e6,Teq_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(MegaTime./1e6,Teq_CH4_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %,'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,45]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
annotation('textbox',[.245 .55-.03 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([3 1 1])
%xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
%L = legend([L1, L2, L4, L5],'GMST Forcing (Judd+2024)','Solar + CO_2 only','Solar + CO_2 + CH_4','Preindustrial GMST (14^oC)','FontSize',20); % 'CO_2 + Low CH_4',
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
%L.AutoUpdate = 'off';
title('Constant Phanerozoic Climate Sensitivity (Judd+2024)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(20,"points") % 14

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])


box on

nexttile

yyaxis left

%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 

L1 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(MegaTime./1e6,Teq_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(MegaTime./1e6,Teq_CH4_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %^,'MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,45]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([3 1 1])
xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
L = legend([L1, L2, L4],'GMST Prior','Solar + CO_2 only','Solar + CO_2 + CH_4','FontSize',12); % 'CO_2 + Low CH_4', ,'PI GMST (14^oC)'
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
title('GMST-Dependent Climate Sensitivity (Wolf+2018)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(12,"points") % 14

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])
box on




figure(20101);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(20101),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
nexttile

yyaxis left

size = 45;

L1 = semilogy(-5,0.565,'o','MarkerFaceColor','k');
hold on
%L2 = plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
%hold on 
L3 = semilogy(-5,0.715,'o','MarkerFaceColor','k');
hold on
%scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
%hold on
%plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
%hold on
%plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
%hold on
%plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on

%L4 = semilogy(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k'); % Beerling for comparison
%hold on

L5 = plot(polyshape(polyx,cat(1,MaxEnv_pCH4.*1e6,flip(MinEnv_pCH4.*1e6))),'FaceColor',colorCH4,'EdgeColor','none');
hold on
L6 = plot(polyshape(polyx,cat(1,MaxEa_pCH4.*1e6,flip(MinEa_pCH4.*1e6))),'FaceColor','c','EdgeColor','none');
hold on
L7 = plot(polyshape(polyx,cat(1,MaxT_pCH4.*1e6,flip(MinT_pCH4.*1e6))),'FaceColor','r','EdgeColor','none');
hold on
L8 = plot(polyshape(polyx,cat(1,MaxCoal_pCH4.*1e6,flip(MinCoal_pCH4.*1e6))),'FaceColor',darkgreen,'EdgeColor','none');
hold on

L9 = semilogy(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','none','LineStyle','-','LineWidth',3); % nominal!
hold on
%semilogy(MegaTime./1e6,MaxCoal_pCH4.*1e6,'Color',darkgreen,'Marker','none','LineStyle','-.','LineWidth',2)
%hold on
%semilogy(MegaTime./1e6,MaxT_pCH4.*1e6,'r','Marker','none','LineStyle','--','LineWidth',2)
%hold on
%semilogy(MegaTime./1e6,MaxEa_pCH4.*1e6,'c','Marker','none','LineStyle',':','LineWidth',2)
%hold on
L10 = semilogy(MegaTime./1e6,MaxO2_pCH4.*1e6,'b','Marker','none','LineStyle',':','LineWidth',2);
hold on
%semilogy(MegaTime./1e6,MaxEnv_pCH4.*1e6,'Color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2)
%hold on

%semilogy(MegaTime./1e6,MinCoal_pCH4.*1e6,'Color',darkgreen,'Marker','none','LineStyle','-.','LineWidth',2)
%hold on
%semilogy(MegaTime./1e6,MinT_pCH4.*1e6,'r','Marker','none','LineStyle','--','LineWidth',2)
%hold on
%semilogy(MegaTime./1e6,MinEa_pCH4.*1e6,'c','Marker','none','LineStyle',':','LineWidth',2)
%hold on
L11 = semilogy(MegaTime./1e6,MinO2_pCH4.*1e6,'b','Marker','none','LineStyle',':','LineWidth',2);
hold on
%semilogy(MegaTime./1e6,MinEnv_pCH4.*1e6,'Color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2)
%hold on

% semilogy(MegaTime./1e6,Ferm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','-','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,maxFerm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minFerm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,NoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','-','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,maxNoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minNoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,FermNoCoal_pCH4.*1e6,'Color',Emcolor,'Marker','none','LineStyle','-','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,maxFermNoCoal_pCH4.*1e6,'Color',Emcolor,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minFermNoCoal_pCH4.*1e6,'Color',Emcolor,'Marker','none','LineStyle','--','LineWidth',2)
% hold on

%plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
%plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on

xline([0],'-k')
hold on
errorbar(-5,0.565,NaN,0.15,'Color','k','CapSize',0)
hold on
annotation('textbox',[.15 .825 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[0.01,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('CH_4 (ppm)')
yticks([0.1, 1, 10, 100]);
yticklabels({'0.1','1','10','100'});
%yticklabels({'0.1','','','','0.5','','','','','1','','','','','','','','','10','','','','50','','','','','100','','','','500','','','','','1000'});
%ylim([-20,320]) % 350
pbaspect([3 1 1])
fontsize(24,"points") % 14
L = legend([L9, L10, L8, L7, L6, L5],'Nominal','pO_2','\Gamma_{coal}','GMST','E_a',...
    'Total','FontSize',24,'NumColumns', 2); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
%,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
L.AutoUpdate = 'off';
%title('Sensitivity Test for yCH_4');

box on

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(24,"points") % 14
%set(gca,'YTickLabel',[]);
yticks([])


nexttile

yyaxis left

size = 45;

semilogy(-5,0.565,'o','MarkerFaceColor','k')
hold on
%plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
%hold on 
semilogy(-5,0.715,'o','MarkerFaceColor','k')
hold on
%scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
%hold on
%plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
%hold on
%plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
%hold on
%plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on

%semilogy(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k') % Beerling for comparison
%hold on
L1 = plot(polyshape(polyx,cat(1,MaxEnv_pCH4.*1e6,flip(MinEnv_pCH4.*1e6))),'FaceColor',colorCH4,'EdgeColor','none');
hold on
semilogy(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','none','LineStyle','-','LineWidth',3) % nominal!
hold on
% semilogy(MegaTime./1e6,MaxEnv_pCH4.*1e6,'Color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,MinEnv_pCH4.*1e6,'Color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2)
% hold on

semilogy(MegaTime./1e6,Ferm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','-','LineWidth',2)
hold on
L2 = plot(polyshape(polyx,cat(1,maxFerm_pCH4.*1e6,flip(minFerm_pCH4.*1e6))),'FaceColor','g','EdgeColor','none');
hold on
% semilogy(MegaTime./1e6,maxFerm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minFerm_pCH4.*1e6,'Color','g','Marker','none','LineStyle','--','LineWidth',2)
% hold on
semilogy(MegaTime./1e6,NoCoal_pCH4.*1e6,'Color','b','Marker','none','LineStyle','-','LineWidth',2)
hold on
L3 = plot(polyshape(polyx,cat(1,maxNoCoal_pCH4.*1e6,flip(minNoCoal_pCH4.*1e6))),'FaceColor','b','EdgeColor','none');
hold on
% semilogy(MegaTime./1e6,maxNoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minNoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
semilogy(MegaTime./1e6,FermNoCoal_pCH4.*1e6,'Color',colorRFtotrange,'Marker','none','LineStyle','-','LineWidth',2)
hold on
L4 = plot(polyshape(polyx,cat(1,maxFermNoCoal_pCH4.*1e6,flip(minFermNoCoal_pCH4.*1e6))),'FaceColor',colorRFtotrange,'EdgeColor','none');
hold on
% semilogy(MegaTime./1e6,maxFermNoCoal_pCH4.*1e6,'Color',Emcolor,'Marker','none','LineStyle','--','LineWidth',2)
% hold on
% semilogy(MegaTime./1e6,minFermNoCoal_pCH4.*1e6,'Color',Emcolor,'Marker','none','LineStyle','--','LineWidth',2)
% hold on

%plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
%plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on

xline([0],'-k')
hold on
errorbar(-5,0.565,NaN,0.15,'Color','k','CapSize',0)
hold on
annotation('textbox',[.15 .35 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[0.01,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('CH_4 (ppm)')
yticks([0.1, 1, 10, 100]);
yticklabels({'0.1','1','10','100'});
%yticklabels({'0.1','','','','0.5','','','','','1','','','','','','','','','10','','','','50','','','','','100','','','','500','','','','','1000'});
%ylim([-20,320]) % 350
pbaspect([3 1 1])
fontsize(24,"points") % 14
L = legend([L1, L2, L3, L4],'Nominal Range',...
    'Fermentation','No \Gamma_{coal}','Fermentation, No \Gamma_{coal}','FontSize',24,'NumColumns', 2); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
%,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
L.AutoUpdate = 'off';
%title('Sensitivity Test for yCH_4');

box on

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
%PhanTransitions;

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(20,"points") % 14
%set(gca,'YTickLabel',[]);
yticks([])





polypCH4_minmax = cat(1,MaxEnv_pCH4.*1e6,flip(MinEnv_pCH4.*1e6));

figure(20111);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(20111),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
nexttile

yyaxis left

semilogy(-5,0.565,'o','MarkerFaceColor','k')
hold on
plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
size = 45;

%scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
%hold on
%plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
%hold on
%plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
%hold on
%plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on

semilogy(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k')
hold on
plot(polyshape(polyx,polypCH4_minmax),'FaceColor',colorCH4,'EdgeColor','none');
hold on
semilogy(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2);%,'MarkerSize',3)
hold on

%plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
%plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
annotation('textbox',[.255 .85-0.031 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
hold on
semilogy(-5,0.715,'o','MarkerFaceColor','k')
hold on
yticks([0.1, 1, 10, 100, 1000]);
yticklabels({'0.1','1','10','100','1000'});
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[0.02,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('CH_4 (ppm)')
%ylim([-20,320]) % 350
pbaspect([2 1 1])
fontsize(16,"points") % 14 Preindustrial yCH_4 (0.565-0.715 ppm)
L = legend('','','CH_4 per Beerling+(2009)','','Revised CH_4','FontSize',16); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
%,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
L.AutoUpdate = 'off';
%title('CH_4 Mixing Ratio (ppm)');

box on

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(16,"points") % 14
%set(gca,'YTickLabel',[]);
yticks([])

nexttile

yyaxis left

plot(polyshape(polyconst,polyKPgrf),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%plot(TotalTime./1e6,rf.Delta_Fs,'color',colorsolar)
%hold on
%scatter(TotalTime(1:16)./1e6,RF_Flux.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
%plot(TotalTime(17:end)./1e6,RF_Flux.Teq_o.t(17:end) - 273.15,'color',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CH4(1:16),size,'v','MarkerFaceColor',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CO2(1:16),size,'v','MarkerFaceColor',colorCO2)
%hold on
%
plot(MegaTime./1e6, rf.Delta_Fs,'color',colorsolar,'LineStyle','-','LineWidth',2,'Marker','none'); % 
hold on 
plot(MegaTime./1e6, RF_CO2.Mid,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
hold on
%plot(TotalTime./1e6, rf.RF_CH4,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
%hold on
%plot(MegaTime./1e6,RF_CH4.Min,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on
%plot(MegaTime./1e6,RF_CH4.Low,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b')
%hold on
plot(MegaTime./1e6,RF_CH4.Mid,'LineStyle','-','color',colorCH4,'LineWidth',2,'Marker','none'); %'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m')
hold on
%plot(MegaTime./1e6,(RF_CH4.Low + rf.RF_CO2 + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none')
%hold on
plot(MegaTime./1e6,(RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',3,'Marker','none')
hold on

%annotation('textbox',[.4 .55-0.022 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(-5,2.7e-7,'ok','MarkerFaceColor','k','MarkerSize',5)
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.RF_N2O,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.RF_N2O,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
annotation('textbox',[.255 .55-0.212 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline(0,'k')
hold on
plot(-5,0,'o','MarkerFaceColor',colorRFtot,'MarkerSize',5)
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[-10,15]) %-9, 12 % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(16,"points") % 14
pbaspect([2 1 1])
xlabel('Age Before Present (Ma)');
%xticklabels([]);
ylabel('Radiative Forcing (W/m^2 Relative to PIM)')
% 'CH_4 (Weak \gamma_T)', 
L = legend('','Solar','CO_2','CH_4',...
        'Total','FontSize',16);
% IN CAPTION, DISCUSS HOW CO2 and CH4 do not account for N2O RF but N2O RF (and total RF) do
%legend('Radiative Forcing from N_2O (low crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (high crustal E_a, high emission scenario)','Radiative Forcing from N_2O (high crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
%title('Modeled Solar and GHG Radiative Forcings')

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(16,"points") % 14
yticks([])

box on








% plot(MegaTime./1e6,(TotalGMST - 14)./((RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs)./4.1124)) % calculates kappa_c needed for GMST prior-posterior match in Eocene








figure(330332);clf
polypGMSTminmax = cat(1,GMST.Max,flip(GMST.Min));
polypTCO2minmax = cat(1,Teq_CO2.Max - 273.15,flip(Teq_CO2.Min - 273.15));
polypTCH4minmax = cat(1,Teq_CH4_CO2.Max - 273.15,flip(Teq_CH4_CO2.Min - 273.15));
polypTCO2minmaxWolf = cat(1,Teq_CO2Wolf.Max - 273.15,flip(Teq_CO2Wolf.Min - 273.15));
polypTCH4minmaxWolf = cat(1,Teq_CH4_CO2Wolf.Max - 273.15,flip(Teq_CH4_CO2Wolf.Min - 273.15));
%polypGMSTminmax = cat(1,GMST.Max,flip(GMST.Min));
polyRFCO2 = cat(1,RF_CO2.Max,flip(RF_CO2.Min));
polyRFCH4 = cat(1,RF_CH4.Max,flip(RF_CH4.Min));
polyRFtot = cat(1,RF_CO2.Max + RF_CH4.Max + rf.Delta_Fs,flip(RF_CO2.Min + RF_CH4.Min + rf.Delta_Fs));

%figure(33055);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(330332),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(3,1,"TileSpacing","compact","Padding","compact");
nexttile

yyaxis left

plot(polyshape(polyconst,polyKPgrf),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%plot(TotalTime./1e6,rf.Delta_Fs,'color',colorsolar)
%hold on
%scatter(TotalTime(1:16)./1e6,RF_Flux.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
%plot(TotalTime(17:end)./1e6,RF_Flux.Teq_o.t(17:end) - 273.15,'color',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CH4(1:16),size,'v','MarkerFaceColor',colorCH4)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.RF_CO2(1:16),size,'v','MarkerFaceColor',colorCO2)
%hold on
%
plot(MegaTime./1e6, rf.Delta_Fs,'color',colorsolar,'LineStyle','-','LineWidth',2,'Marker','none'); % no uncertainty resolved, so just a line
hold on 
plot(polyshape(polyx,polyRFCO2),'FaceColor',colorCO2,'EdgeColor',colorCO2)
%plot(MegaTime./1e6, RF_CO2.Mid,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
hold on
plot(MegaTime./1e6, RF_CO2.Mid,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
hold on
%plot(TotalTime./1e6, rf.RF_CH4,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % includes N2O spectral overlap
%hold on
%plot(MegaTime./1e6,RF_CH4.Min,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on
%plot(MegaTime./1e6,RF_CH4.Low,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b')
%hold on
plot(polyshape(polyx,polyRFCH4),'FaceColor',colorCH4,'EdgeColor',colorCH4)
%plot(MegaTime./1e6,RF_CH4.Mid,'LineStyle','-','color',colorCH4,'LineWidth',2,'Marker','none'); %'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m')
hold on
plot(MegaTime./1e6,RF_CH4.Mid,'LineStyle','-','color',colorCH4,'LineWidth',2,'Marker','none'); %'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m')
hold on
%plot(MegaTime./1e6,(RF_CH4.Low + rf.RF_CO2 + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none')
%hold on
plot(polyshape(polyx,polyRFtot),'FaceColor',colorRFtot,'EdgeColor',colorRFtot)
%plot(MegaTime./1e6,(RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',3,'Marker','none')
hold on
plot(MegaTime./1e6,(RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none')
hold on

annotation('textbox',[.245 .84 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .55-.03 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
%annotation('textbox',[.4 .55-0.022 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(-5,2.7e-7,'ok','MarkerFaceColor','k','MarkerSize',5)
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.RF_N2O,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.RF_N2O,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
xline(0,'k')
hold on
plot(-5,0,'o','MarkerFaceColor',colorRFtot,'MarkerSize',5)
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[-10,20]) %-9, 12 % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14
pbaspect([3 1 1])
%xlabel('Age Before Present (Ma)');
%xticklabels([]);
ylabel('Radiative Forcing (W/m^2 Relative to PIM)')
% 'CH_4 (Weak \gamma_T)', 
L = legend('','Solar','CO_2','','CH_4','',...
        'Total','','FontSize',12);
% IN CAPTION, DISCUSS HOW CO2 and CH4 do not account for N2O RF but N2O RF (and total RF) do
%legend('Radiative Forcing from N_2O (low crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (high crustal E_a, high emission scenario)','Radiative Forcing from N_2O (high crustal E_a, low emission scenario)',...
%    'Radiative Forcing from N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
title('Modeled Solar and GHG Radiative Forcings')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])

box on


nexttile

yyaxis left
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)

L1 = plot(polyshape(polyx,polypGMSTminmax),'FaceColor',colorT,'EdgeColor',colorT); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L11 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(polyshape(polyx,polypTCO2minmax),'FaceColor',colorCO2,'EdgeColor',colorCO2);%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
L21 = plot(MegaTime./1e6,Teq_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(polyshape(polyx,polypTCH4minmax),'FaceColor',colorRFtot,'EdgeColor',colorRFtot); %,'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
L41 = plot(MegaTime./1e6,Teq_CH4_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %,'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
%hold on 
%plot(MegaTime./1e6,Teq_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
%hold on
%plot(MegaTime./1e6,Teq_CH4_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %^,'MarkerEdgeColor','m','MarkerFaceColor','m');
%hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 

%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on

xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,50]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%annotation('textbox',[.245 .84 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
annotation('textbox',[.245 .55-.03 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([3 1 1])
%xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
%L = legend([L1, L2, L4, L5],'GMST Forcing (Judd+2024)','Solar + CO_2 only','Solar + CO_2 + CH_4','Preindustrial GMST (14^oC)','FontSize',20); % 'CO_2 + Low CH_4',
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
%L.AutoUpdate = 'off';
title('Constant Phanerozoic Climate Sensitivity (Judd+2024)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(12,"points") % 14

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])

box on


nexttile

yyaxis left
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)

L1 = plot(polyshape(polyx,polypGMSTminmax),'FaceColor',colorT,'EdgeColor',colorT); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L12 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(polyshape(polyx,polypTCO2minmaxWolf),'FaceColor',colorCO2,'EdgeColor',colorCO2);%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
L22 = plot(MegaTime./1e6,Teq_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(polyshape(polyx,polypTCH4minmaxWolf),'FaceColor',colorRFtot,'EdgeColor',colorRFtot); %^,'MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
L42 = plot(MegaTime./1e6,Teq_CH4_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %^,'MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 


xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,50]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%annotation('textbox',[.245 .84 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%annotation('textbox',[.245 .55-.03 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
annotation('textbox',[.245 .2 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([3 1 1])
xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
L = legend([L1, L2, L4, L5],'GMST Prior','Solar + CO_2 only','Solar + CO_2 + CH_4','FontSize',12); % 'CO_2 + Low CH_4', ,'PI GMST (14^oC)'
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
title('GMST-Dependent Climate Sensitivity (Wolf+2018)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(12,"points") % 14

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(12,"points") % 14
yticks([])


box on





figure(330345);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(330345),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(2,1,"TileSpacing","compact","Padding","compact");


nexttile

yyaxis left

%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 

L1 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(MegaTime./1e6,Teq_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(MegaTime./1e6,Teq_CH4_CO2.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %,'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,45]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
annotation('textbox',[.22 .85-0.031 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([5 2 2])
%xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
%L = legend([L1, L2, L4, L5],'GMST Forcing (Judd+2024)','Solar + CO_2 only','Solar + CO_2 + CH_4','Preindustrial GMST (14^oC)','FontSize',20); % 'CO_2 + Low CH_4',
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
%L.AutoUpdate = 'off';
title('Constant Phanerozoic Climate Sensitivity (Judd+2024)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(20,"points") % 14

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(20,"points") % 14
yticks([])


box on

nexttile

yyaxis left

%scatter(TotalTime(1:16)./1e6,rf.Teq_CO2.t(1:16) - 273.15,size,'^','MarkerFaceColor',colorCO2)
%hold on
%scatter(TotalTime(1:16)./1e6,rf.Teq_o.t(1:16) - 273.15,size,'v','MarkerFaceColor',colorCH4)
%hold on
% rf.Teq_CH4_CO2.(HL)
L5 = plot(-2.5,14,'o','MarkerFaceColor',colorT,'MarkerSize',5);
hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 

L1 = plot(TotalTime./1e6, TotalGMST,'color',colorT,'LineStyle','-','LineWidth',2,'Marker','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
L2 = plot(MegaTime./1e6,Teq_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorCO2,'Marker','none');%d,'MarkerEdgeColor',colorCO2,'MarkerFaceColor',colorCO2);
hold on
%L44 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Min - 273.15,'LineStyle','-','color','k','Marker','d','MarkerEdgeColor','k','MarkerFaceColor','none');
%hold on
%L3 = plot(MegaTime./1e6,rf.Teq_CH4_CO2.Low - 273.15,'LineStyle','-','color','b','Marker','v','MarkerEdgeColor','b','MarkerFaceColor','b');
%hold on
L4 = plot(MegaTime./1e6,Teq_CH4_CO2Wolf.Mid - 273.15,'LineStyle','-','LineWidth',2,'color',colorRFtot,'Marker','none'); %^,'MarkerEdgeColor','m','MarkerFaceColor','m');
hold on
%L3 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_CO2.t - 273.15,'color',colorCO2,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%L4 = plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCH4,'LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(Gamma0.time{1}./1e6),(Gamma0.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(4500-(Gamma100.time{1}./1e6),(Gamma100.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',1,'Marker','none');
%hold on
%plot(polyshape(polyx,polyClim),'FaceColor',colorRFtotrange)
%hold on
%L1 = plot(4500-(GammaRef.time{1}./1e6),(GammaRef.RF_Flux.Teq_N.t - 273.15),'color',colorRFtot,'LineStyle','-','LineWidth',2,'Marker','none');
%hold on
%annotation('textbox',[.4 .25-0.033-0.0045 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on
%plot(4500-(Gamma0.time{1}./1e6), Gamma0.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(Gamma100.time{1}./1e6), Gamma100.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','-','LineWidth',2,'Marker','none'); % 
%hold on
%plot(4500-(GammaMaxMin.time{1}./1e6), GammaMaxMin.RF_Flux.Teq_N.t - 273.15,'c','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
%plot(4500-(GammaMinMax.time{1}./1e6), GammaMinMax.RF_Flux.Teq_N.t - 273.15,'r','LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on 
xline(0,'k')
%hold on
%plot(4500-(GammaRef.time{1}./1e6), GammaRef.RF_Flux.Teq_o.t - 273.15,'color',colorCO2,'LineStyle','--','LineWidth',2,'Marker','none'); % 
%hold on
set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375],'ylim',[0,45]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
annotation('textbox',[.22 .55-0.212 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on

pbaspect([5 2 2])
xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)') % ,'Solar + CO_2 + Weak CH_4 \gamma_T',
L = legend([L1, L2, L4],'GMST Prior','Solar + CO_2 only','Solar + CO_2 + CH_4','FontSize',20); % 'CO_2 + Low CH_4', ,'PI GMST (14^oC)'
%L = legend([L2, L3, L4, L1],'GMST Forcing','CO_2 only','CO_2 and CH_4 only',...
%    'CO_2, CH_4, and N_2O','FontSize',12,'Position',[0.5+0.01 0.1+0.01 0.1 0.08]); % 
%L = legend('','','','CO_2, CH_4, and N_2O','GMST Forcing Prior (J+24)','CO_2 only','CO_2 and CH_4 only','FontSize',10); % ,...
 %   'Climate History with CO_2, CH_4, and N_2O (high crustal E_a, high emission scenario)','Climate History with CO_2, CH_4, and N_2O (high crustal E_a, low emission scenario)',...
 %   'Climate History with CO_2, CH_4, and N_2O (low crustal E_a, high emission scenario)','FontSize',12)
L.AutoUpdate = 'off';
title('GMST-Dependent Climate Sensitivity (Wolf+2018)')
% (\Gamma_c = 1.7 K/(W/m^2)), assuming Beerling et al. [2009] yCH_4)
fontsize(20,"points") % 14

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-5,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(20,"points") % 14
yticks([])
box on


















%polypCH4_minmax = cat(1,MaxEnv_pCH4.*1e6,flip(MinEnv_pCH4.*1e6));

figure(201115);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(201115),'defaultAxesColorOrder',[lc; rc]);

tiledlayout(2,1,"TileSpacing","compact","Padding","compact");

polyFCH4em_minmax = cat(1,16.04.*FluxCH4_max./1e3,flip(16.04.*FluxCH4_min./1e3));

nexttile

yyaxis left

size = 45;

semilogy(-5,16.04.*1.31359e13./1e15,'o','MarkerFaceColor','k')
hold on
plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'k','Marker','d','LineWidth',1.25)
%Flux_CH4_emissions_BC89_hires_std./1e12
%hold on timeslices./1e6,FluxCH4_mid
LA = semilogy(timeFluxB09./1e6,16.04.*FluxCH4_B09./1e3,'k','Marker','o','LineStyle','-');
hold on
%[Nx,Ny] = boundary(polyshape(polyx2,polyFCH4em_minmax));
%patch(Nx,Ny,Emcolor) % ,'FaceColor',Emcolor,'EdgeColor',Emcolor
plot(polyshape(polyx,polyFCH4em_minmax),'FaceColor',colorCH4,'EdgeColor','none');
hold on
LB = semilogy(MegaTime./1e6,16.04.*FluxCH4_mid./1e3,'color',colorCH4,'Marker','none','LineWidth',2,'Linestyle','-'); %,'MarkerSize',3 'none' for no marker
hold on
%plot(polyshape([48, 48, 56, 56],cat(1,0.656, 0.909, 0.909, 0.656).'),'FaceColor','k','EdgeColor','none'); % cf. Wilton+2019 abstract, Ypresian Eocene (56-48 Ma) wetland fraction is 2-2.5x higher than reference modern value 
LC = semilogy(52,0.656,'MarkerFaceColor','k','Marker','^','LineStyle','none'); % low bounds from Wilton+2019
hold on 
semilogy(52,0.909,'MarkerFaceColor','k','Marker','v','LineStyle','none') % high bounds from Wilton+2019
hold on 
LD = semilogy(55,1.1426,'MarkerFaceColor','k','Marker','square','LineStyle','none'); % From Table S1, Beerling+2011 4x CO2 Eocene wetland CH4 flux total
hold on 
semilogy(90,0.857,'MarkerFaceColor','k','Marker','square','LineStyle','none') % From Table S1, Beerling+2011 4x CO2 Cretaceous wetland CH4 flux total
hold on 
%
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':')
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_max./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':') % 'none' for no marker
%hold on
%
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on

xline([0],'-k')
yticks([0.1, 1, 10, 100, 1000]);
yticklabels({'0.1','1','10','100','1000'});
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('CH_4 Emissions (Pg CH_4/yr)')
annotation('textbox',[.25 .85-0.031 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
%ylim([-80,1400]) % [-35,650]
%ylim([1,7000]) % [-35,650]
ylim([0.02,100])
%title('Global CH_4 Emissions')
pbaspect([2 1 1])
fontsize(16,"points") % 14   'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T per Rubisco, GMST per J+24)',  'Weak \gamma_T per Zhu+(2014) & Judd+(2024)',
L = legend([LB, LA, LD, LC],'EarthCH_4','Beerling+(2009)',...
    'Beerling+(2011)','Wilton+(2019)',...
    'FontSize',16); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
L.AutoUpdate = 'off';


yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(16,"points") % 14
box on


nexttile

yyaxis left

semilogy(-5,0.565,'o','MarkerFaceColor','k')
hold on
plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
size = 45;

%scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
%hold on
%plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
%hold on
%plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
%hold on
%plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
%hold on

semilogy(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k')
hold on
plot(polyshape(polyx,polypCH4_minmax),'FaceColor',colorCH4,'EdgeColor','none');
hold on
semilogy(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2);%,'MarkerSize',3)
hold on

semilogy(55,3.614,'MarkerFaceColor','k','Marker','square','LineStyle','none') % From Table 1, Beerling+2011 4x CO2 Eocene pCH4 (not PI isoprene)
hold on 
semilogy(90,3.304,'MarkerFaceColor','k','Marker','square','LineStyle','none') % From Table 1, Beerling+2011 4x CO2 Cretaceous pCH4 (not PI isoprene)
hold on 

%plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
%plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
%hold on
annotation('textbox',[.25 .55-0.212-.007 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
hold on
semilogy(-5,0.715,'o','MarkerFaceColor','k')
hold on
yticks([0.1, 1, 10, 100, 1000]);
yticklabels({'0.1','1','10','100','1000'});
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[0.02,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); 
ylabel('CH_4 (ppm)')
%ylim([-20,320]) % 350
pbaspect([2 1 1])
fontsize(16,"points") % 14 Preindustrial yCH_4 (0.565-0.715 ppm)
%L = legend('','','yCH_4 per Beerling+(2009)','','Revised yCH_4','FontSize',16); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
%,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
%L.AutoUpdate = 'off';
%title('CH_4 Mixing Ratio (ppm)');

box on

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
fontsize(16,"points") % 14
%set(gca,'YTickLabel',[]);
yticks([])














% % Duplicate
% figure(2011145);clf
% 
% lc = [0 0 0];
% rc = [0 0 0];
% set(figure(2011145),'defaultAxesColorOrder',[lc; rc]);
% 
% tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
% 
% polyFCH4em_minmax = cat(1,16.04.*FluxCH4_max./1e3,flip(16.04.*FluxCH4_min./1e3));
% 
% nexttile
% 
% yyaxis left
% 
% size = 45;
% 
% semilogy(-5,16.04.*1.31359e13./1e15,'o','MarkerFaceColor','k')
% hold on
% plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
% hold on 
% %scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
% %hold on
% %scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
% %hold on
% %scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
% %hold on
% 
% %plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
% %hold on
% %plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
% %hold on
% %plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'k','Marker','d','LineWidth',1.25)
% %Flux_CH4_emissions_BC89_hires_std./1e12
% %hold on timeslices./1e6,FluxCH4_mid
% semilogy(timeFluxB09./1e6,16.04.*FluxCH4_B09./1e3,'k','Marker','o','LineStyle','-')
% hold on
% %[Nx,Ny] = boundary(polyshape(polyx2,polyFCH4em_minmax));
% %patch(Nx,Ny,Emcolor) % ,'FaceColor',Emcolor,'EdgeColor',Emcolor
% plot(polyshape(polyx,polyFCH4em_minmax),'FaceColor',colorCH4,'EdgeColor','none');
% hold on
% semilogy(MegaTime./1e6,16.04.*FluxCH4_mid./1e3,'color',colorCH4,'Marker','none','LineWidth',2,'Linestyle','-') %,'MarkerSize',3 'none' for no marker
% hold on
% %
% %plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':')
% %hold on
% %plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_max./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':') % 'none' for no marker
% %hold on
% %
% %errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
% %hold on
% 
% xline([0],'-k')
% yticks([0.1, 1, 10, 100, 1000]);
% yticklabels({'0.1','1','10','100','1000'});
% %plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
% set(gca,'XDir','reverse');
% set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
% %set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
% %set(gca,'yaxislocation','left')
% %xlabel('Age Before Present (Ma)'); 
% ylabel('CH_4 Emissions (Pg CH_4/yr)')
% annotation('textbox',[.25 .85-0.031 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
% hold on
% %ylim([-80,1400]) % [-35,650]
% %ylim([1,7000]) % [-35,650]
% ylim([0.02,100])
% %title('Global CH_4 Emissions')
% pbaspect([2 1 1])
% fontsize(16,"points") % 14   'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T per Rubisco, GMST per J+24)',  'Weak \gamma_T per Zhu+(2014) & Judd+(2024)',
% L = legend('','','CH_4 Emissions no Gcoal','',...
%     'Revised CH_4 Emissions',...
%     'FontSize',16); % 14,'Modern pN_2O (337 ppb)'     
% % 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
% L.AutoUpdate = 'off';
% 
% 
% yyaxis right
% 
% %geotimescale_Mills_JFHmod_375Ma;
% %hold on
% PhanTransitions;
% set(gca,'YTickLabel',[]);
% yticks([]);
% 
% set(gca,'XDir','reverse');
% set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
% fontsize(16,"points") % 14
% box on
% 
% 
% nexttile
% 
% yyaxis left
% 
% semilogy(-5,0.565,'o','MarkerFaceColor','k')
% hold on
% plot(polyshape(polyconst,polyKPglog),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
% hold on 
% size = 45;
% 
% %scatter(TotalTime(1:16)./1e6,Total_pCH4(1:16).*1e6,size,'o','MarkerFaceColor','b')
% %hold on
% %scatter(AnchorTime./1e6,Phan_pCH4_hi.*1e6,size,'^','MarkerFaceColor','m')
% %hold on
% %scatter(AnchorTime./1e6,Phan_pCH4_noT.*1e6,size,'d','MarkerFaceColor','g')
% %hold on
% %plot(timeslices./1e6,Phan_pCH4_hires_hi.*1e6,'m','Marker','none') % 'none' for no marker
% %hold on
% %plot(TotalTime(17:end)./1e6,Total_pCH4(17:end).*1e6,'b','Marker','none','LineStyle','-')
% %hold on
% %plot(MegaTime./1e6,Total_pCH4_min.*1e6,'k','Marker','d','LineStyle','-','MarkerEdgeColor','k','MarkerFaceColor','none')
% %hold on
% 
% semilogy(TimePhanCH4./1e6,PhanCH4vals./1e3,'k','Marker','o','LineStyle','-','MarkerEdgeColor','k')
% hold on
% plot(polyshape(polyx,polypCH4_minmax),'FaceColor',colorCH4,'EdgeColor','none');
% hold on
% semilogy(MegaTime./1e6,Total_pCH4_high.*1e6,'color',colorCH4,'Marker','none','LineStyle','-','LineWidth',2);%,'MarkerSize',3)
% hold on
% 
% %plot(MegaTime./1e6,MaxEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
% %hold on
% %plot(MegaTime./1e6,MinEnv_pCH4.*1e6,'m','Marker','none','LineStyle',':','LineWidth',2)
% %hold on
% annotation('textbox',[.25 .55-0.212-.007 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
% hold on
% xline([0],'-k')
% hold on
% semilogy(-5,0.715,'o','MarkerFaceColor','k')
% hold on
% yticks([0.1, 1, 10, 100, 1000]);
% yticklabels({'0.1','1','10','100','1000'});
% %plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
% set(gca,'XDir','reverse');
% set(gca,'xlim',[-10,375],'ylim',[0.02,1000]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
% %set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
% %set(gca,'yaxislocation','left')
% xlabel('Age Before Present (Ma)'); 
% ylabel('yCH_4 (ppm)')
% %ylim([-20,320]) % 350
% pbaspect([2 1 1])
% fontsize(16,"points") % 14 Preindustrial yCH_4 (0.565-0.715 ppm)
% L = legend('','','yCH_4 no Gcoal','','Revised yCH_4','FontSize',16); % 14,'Modern pN_2O (337 ppb)' 'yCH_4 (Weak \gamma_T per Zhu+2014)',
% %,'Revised Phanerozoic yCH_4 (Strong \gamma_T)',...
% %    'Revised Phanerozoic yCH_4 (No \gamma_T)' 'Revised Phanerozoic yCH_4 (Weak \gamma_T for C-cycle)',
% L.AutoUpdate = 'off';
% %title('CH_4 Mixing Ratio (ppm)');
% 
% box on
% 
% yyaxis right
% 
% geotimescale_Mills_JFHmod_375Ma;
% hold on
% PhanTransitions;
% 
% set(gca,'XDir','reverse');
% set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
% %fontsize(18,"points") % 14
% fontsize(16,"points") % 14
% %set(gca,'YTickLabel',[]);
% yticks([])
% 
% 






polyf_RF_CH4_minmax = cat(1,100.*RF_CH4.Max./(RF_CH4.Max + RF_CO2.Max + rf.Delta_Fs),flip(100.*RF_CH4.Min./(RF_CH4.Min + RF_CO2.Min + rf.Delta_Fs)));

figure(32303235);clf

plot(polyshape(polyx,polyf_RF_CH4_minmax),'FaceColor','m','EdgeColor','none');
hold on
plot(MegaTime./1e6,100.*RF_CH4.Mid./(RF_CH4.Mid + RF_CO2.Mid + rf.Delta_Fs),'color','m','LineStyle','-','LineWidth',3,'Marker','none')
hold on
% plot(MegaTime./1e6,100.*RF_CH4.Max./(RF_CH4.Max + RF_CO2.Max + rf.Delta_Fs),'color','m','LineStyle','--','LineWidth',2,'Marker','none')
% hold on
% plot(MegaTime./1e6,100.*RF_CH4.Min./(RF_CH4.Min + RF_CO2.Min + rf.Delta_Fs),'color','m','LineStyle','--','LineWidth',2,'Marker','none')
% hold on

set(gca,'XDir','reverse');
set(gca,'xlim',[30,150]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
xlabel('Age Before Present (Ma)'); 
ylabel('CH_4 Fraction of Total Radiative Forcing (%)')
ylim([0,100]) % 350
pbaspect([2 1 1])
fontsize(20,"points") % 14
box on
