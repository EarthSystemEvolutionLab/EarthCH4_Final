%% %% Revised Wetland CH4 Emission Reconstruction For the Phanerozoic %% %%
% Revising Beerling et al. (2009) method for estimating wetland CH4
% emissions over the Phanerozoic, we apply a recent GMST forcing (Judd et
% al. 2024) and a revised T-sensitivity for methane emissions based on
% methanogenesis metabolic sensitivities (cf. Conrad 2023; Ho et al. 2025;
% Liu et al. 2025) to compute wetland CH4 emission fluxes over the Phanerozoic.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Introductory Material

% scenario codes
Constant.invalT = 1; % NOMINAL CASE
Constant.invalCO2 = 1;
Constant.invalO2 = 1;

% plot colors
colorCO2 = [120,94,240]./255;
colorCH4 = [254,97,0]./255;
colorsolar = [255,176,0]./255;
darkgreen = [0,130,0]./255;
KPgcolor = [255, 190, 0]./255;

% basic constants

% Thermodynamic constants
Constant.R = 8.3144; % J./mol*K; ideal gas constant = n*kB
Constant.N_avogadro = 6.022e23; % Avogadro's number (molecules./mol)

T_ref = 287.15; % 14 C preindustrial modern reference temperature (all GHGs at PIM reference values, modern S(t))

% spatiotemporal scaling factors!
s_yr = (60.*60.*24.*365.25); % seconds per year
SA_Earth = 5.101e14; % m^2; earth surface area, from EONS (Horne and Goldblatt 2024), essentially invariant over time (cf. Scotese et al. 2021a)
photochem_fluxscaling = s_yr.*1e4.*SA_Earth./Constant.N_avogadro; % s*cm^2*mol CH4/(yr*global area*molecules CH4); converts photochem flux units (molec/cm^2/s) to mol/yr

%% %% Load Forcing Inputs
% Judd et al. (2024): pCO2, GMST (both back to 485 Ma)

% Beerling et al. (2009): relative coal depositional rates, Pliocene and
% Preindustrial Modern (PIM) wetland CH4 emissions and wetland areas

% Tierney et al. (2025): mid-Pliocene reference GMST roughly 4 deg C warmer
% than PIM (assumedly 14 C, so 18 C in mid Pliocene)

% Loads primary input structure
In_data = load('biome_output.mat');

% Additional: two scenarios to proxy wetland area globally over Phanerozoic
%Macrostrat_fWetData = readmatrix("All_CSV_files\Husson_Macrostrat_coalpeat_fraction_seds_NAm_ROUGHwpd.csv");
%Macrostrat_fWetTime = (Macrostrat_fWetData(:,1)).*1e6; % Ma to years
%Macrostrat_fWet_Area = Macrostrat_fWetData(:,2); % fraction of seds in N America which are coal/peat, extrapolated global FOR NOW - Husson 2026, personal comm.
%Macrostrat_fWet_Area(Macrostrat_fWet_Area < 0) = 0; % avoids negative vals

BC89_fWetData = readmatrix("All_CSV_files/BernerCanfield1989_midAge_vs_coalFractionContClastics.xlsx");
BC89_fWetTime = (BC89_fWetData(:,1)).*1e6; % Ma to years, (midpoints, with 0 Ma tacked on to duplicate Pliocene values assuming marginal wetland area change - cf. Hopcroft + 2020
BC89_fWet_Area = BC89_fWetData(:,2); % f_wetland of land area, assuming that coal fraction of total continental seds (clastics + coal) is representative of land distribution
BC89_fWet_Area(BC89_fWet_Area < 0) = 0; % avoids negative vals

%% GMST
PhanAirTTime = In_data.PhanBiomes.JuddT_time(4:85).*1e6; % in yrs BP;
% (4:85) removes NaN and < 1 Ma BP data for replacement with Snyder data (815 ka to present)
PhanAirTTime = rmmissing(PhanAirTTime);
% values are in deg C (temperature)
if Constant.invalT == 1
    PhanAirT_C = In_data.PhanBiomes.JuddT_C(4:85);
elseif Constant.invalT == 2
    PhanAirT_C = In_data.PhanBiomes.JuddT_C_16perc(4:85);
elseif Constant.invalT == 3
    PhanAirT_C = In_data.PhanBiomes.JuddT_C_84perc(4:85);
end
PhanAirT_C = rmmissing(PhanAirT_C);

% Last 815 ka from Snyder 2016 differential GAST reconstruction
Snydertime = 1e3.*In_data.PhanBiomes.SnyderTime; % converted to years BP
SnyderTemp = 14 + In_data.PhanBiomes.SnyderdT; % in deg C;  + 273.15; % 14 C; 
%  resolves reasonable ~287 degK at 1000 years ago, similar to PIM assumedly - could do 13.8 as reference but not very different

GMST_t = cat(1,PhanAirTTime,Snydertime,0); % years BP
GMST_C = cat(1,PhanAirT_C,SnyderTemp,14); % deg C - added PIM with 14 C reference T
GMST_traw = cat(1,PhanAirTTime,0); % years BP

% 16th and 84th percentiles of GMST reconstruction (high and low
% uncertainty bounds)
GMST_C16raw = cat(1,rmmissing(In_data.PhanBiomes.JuddT_C_16perc(4:85)),14);
GMST_C84raw = cat(1,rmmissing(In_data.PhanBiomes.JuddT_C_84perc(4:85)),14);

%dtime = flip(0:1e6:500e6).';

%GMST_t = dtime;
%GMST_C = interp1(GMST_t_init,GMST_C_init,dtime,'pchip'); % interpolates GMST at 1 Ma intervals for use hereafter

% equatorial temperatures from Judd+2024 for interpolation and input into
% Phanerozoic pCH4 scripts - will need to use best-fit regression to
% estimate tropical SAT during Pleistocene, using function fit to Judd+2024
% data over Phanerozoic (375 Ma)
PhanTeq_t = cat(1,0,In_data.PhanBiomes.JuddTeq_time.*1e6);
PhanTeq_C = cat(1,In_data.PhanBiomes.JuddTeq_C(1),In_data.PhanBiomes.JuddTeq_C); 
% repeats final value of early Holocene to avoid unconstrained extrapolation
% NOT USED - GMST preferred herein

% pCO2

%TimePhanCO2 = 1e6.*(In_data.PhanBiomes.CO2PhanTime); % converted from Ma to yrs BP
%TimePhanCO2 = rmmissing(TimePhanCO2);
%if Constant.invalCO2 == 1 % nominal (median, 50th percentile)
%    PhanCO2vals = 1e-6.*In_data.PhanBiomes.CO2Phanppm; % transformation from ppm to atm abundance accounted for
%elseif Constant.invalCO2 == 2 % low (16th percentile, -sigma)
%    PhanCO2vals = 1e-6.*In_data.PhanBiomes.CO2Phanppm_16perc; 
%elseif Constant.invalCO2 == 3 % high (84th percentile, +sigma)
%    PhanCO2vals = 1e-6.*In_data.PhanBiomes.CO2Phanppm_84perc; 
%end
%PhanCO2vals = rmmissing(PhanCO2vals);
% not currently used

%% pO2 
% Mills et al. 2023 Phanerozoic pO2 consensus curve, assuming 1 atm = total atmo pressure 
% (thus can treat % atm in results as equal to pO2 in atm units)
pO2Phantime = In_data.PhanBiomes.MillspO2_time.*1e6; % converts to yrs BP

if Constant.invalO2 == 1 
    pO2_atm = 0.01.*In_data.PhanBiomes.MillspO2_vals; % rescales from % atm to fractional atm (partial pressure or pO2)
elseif Constant.invalO2 == 2
    pO2_atm = 0.01.*In_data.PhanBiomes.MillspO2_Min_vals; % - NeoprotO2scaling.*Constant.Neoprot_pO2_frac;
elseif Constant.invalO2 == 3
    pO2_atm = 0.01.*In_data.PhanBiomes.MillspO2_Max_vals; % - NeoprotO2scaling.*Constant.Neoprot_pO2_frac;
end
% converts from % atm to atm fraction 

% pO2 max and min
pO2_minraw = 0.01.*In_data.PhanBiomes.MillspO2_Min_vals;
pO2_maxraw = 0.01.*In_data.PhanBiomes.MillspO2_Max_vals;

%% B+09 CH4 Emissions

% extracted using figure digitization with WebPlotDigitizer (Rohatgi),
% hereafter called WPD throughout these scripts

In_data.PhanBiomes.CH4flux(1) = 0; % erases unrealistic negative flux from imperfect WPD
%In_data.PhanBiomes.CH4flux(end) = 210.7; % Tg CH4/yr simulated for PIM from Beerling et al. 2009 Table 2, replaces value from terminal flux reconstruction
CH4fluxtime = (flip(0:1e7:390e6)).'; % years BP: every 10 Ma, approximates x-values from WPD
%cat(1,0,4.5e9-In_data.PhanBiomes.CH4fluxTime); % reorients time to correct direction from Hadean
CH4flux = In_data.PhanBiomes.CH4flux.*(1e12./16.04); % rescaled to moles CH4/yr from Tg/yr
% avoids interpolation problems
%CH4_emissions = interp1(CH4fluxtime,CH4flux,Land_forcings.xq,'pchip').*(1e12./16.04); 
% interpolates CH4 wetland emissions through time and rescales from Tg
% CH4/yr to moles CH4/yr using MM = 16.04 g/mol
%Land_forcings.CH4_emissions(Land_forcings.CH4_emissions < 0) = 0; % avoids negative emissions

time_CH4FluxB09 = CH4fluxtime(2:end)./1e9; % converts to Ga, removes initial 0 which will cause model bugs later
FluxCH4_B09 = CH4flux(2:end)./(0.4.*1e11.*photochem_fluxscaling); % converts into multiples of base-flux (4e10 molec/cm^2/s in photochem = roughly 10.68 Tmol/yr herein)
FluxCH4_B09(end) = 1.23; % 1.22995414 = 1.23 * (0.4.*26.7e12) mol CH4/yr from wetlands global
% adds in reasonable approximation of PI wetland flux of Beerling+2009
% Table 2 - prior value unreasonably tuned down

%% %% Compute Revised CH4 Wetland Emission Flux over Phanerozoic

% We ignore marine and gut emissions herein, only track dominant wetland CH4
% sources over Phanerozoic since evolution of land plants capable of
% forming coal swamps

% T-DEPENDENCE FUNCTION PARAMETERS

Q10_CH2O_photo = 2; % typical metabolic Q10 formulation for photosynthesis/CH2O fixation by Rubisco/respiration and decay, 
% for substrate-limited CH4genesis with plant T-response controlling
% substrate availability across climate changes - cf. Galmes+2016, Maranon+2018
% NOT USED 

Q10_CH4 = 2.48; % 4.1; % 4.1; % Zhu et al. (2014) mean global wetland Q10
% for CH4 emissions; see Ho et al. 2025 ref. to Zhu et al. 2014 ESM model mean value
% ultimately scaled up by mid Pliocene-referenced (18 C, not 14.3 C as in
% original paper) T-dependence (Q10 = 5.4!) 
% Now using Q10 of 4.1 based on Conrad 2023, Ho+2025 (cf. Liu+2025, van Hulzen +1999) 
% NOT USED

Ea_CH4 = 100000; % NOMINAL; roughly 100 kJ/mol for methanogenesis (on archaeal population, community, and ecosystem scales with 
% remarkable consistency across scales) per Conrad 2023, consistent broadly 
% with ~4.1 mean high-T Q10 value in Ho et al. 2025 though less extreme in
% slope
Ea_minCH4 = 79000; % see below
Ea_maxCH4 = 122000; % per Conrad 2023 (ref. Yvon-Durocher 2014), range in mean CH4 Ea vals for culture/community/field systems from 79-122 kJ/mol
Ea_fermentation = 48900; % 48.9 kJ/mol avg Ea for marine sedimentary fermentation (Weston and Joye 2005), 
% assuming that fermentation rate T-response can control methanogenesis T-response via substrate limitation (cf. Valentine+1994)

T_ref_CH4 = 14; % referenced to PI GMST
% assuming PIM long-term reference GMST of 14 C (cf. Beerling et al. 2009; Nema et al. 2012)

% REFERENCE EMISSIONS (PI)
PIWetlandRefEmissionCH4 = (200e12./(16.04)); % using a typical natural modern wetland emission flux of around 200 Tg/yr (0.2 Gt/yr)
% (in reasonable ranges for modern natural and PI - see Mitsch + 2013 and
% refs, Beerling+2009, Saunois+2024

PI_otherCH4Emissions = (50e12./(16.04)); % geological, marine, insect, wildfire, and other sources as constant over deep time
% ignore potential GMST and pO2 dependence feedbacks for termites/marine anoxia/clathrates/wildfire frequency, etc.
% cf. Global Carbon Project, Methane Budget 2024, Saunois+2024
% (citation at: https://www.globalcarbonproject.org/methanebudget/24/publications.htm)

%% Compute CH4 Emissions over Phanerozoic at 1 Ma resolution 
% timesteps for 1-Ma resolution over past 375 Ma
timeslices = flip(0:1e6:375e6).'; % flip(3.5e6:1e6:375e6).'; flip(3e6:1e6:375e6).'

Phanero_landA = interp1(4.5e9-In_data.PhanBiomes.LandTime,In_data.PhanBiomes.LandA,timeslices,'pchip').*SA_Earth./1e12; 
% Earth land surface area over Phanerozoic in Mkm^2, Scotese+2021 Fig 21 red line (fraction land)

Phanero_fwet_BC89 = interp1(BC89_fWetTime,BC89_fWet_Area,timeslices,'pchip'); % Phanerozoic history of relative coal wetland extent

% EXTRAPOLATES 5-Ma Pliocene wetland area through PI Holocene, avoids unrealistically low
% wetland areal fraction due to ultra-high Quaternary/Pleistocene sediment
% preservation bias (also Pliocene and PI wetlands may have had similar
% extent, cf. Hopcroft+2020 SM) - comparable to BC89 data with last 5.5 Ma
% assumedly constant in wetland area
Phanero_fwet_PlioPI = interp1(BC89_fWetTime,BC89_fWet_Area,5.5e6,'pchip'); % corresponds to BC89 final mid-age point
Phanero_fwet_BC89(Phanero_fwet_BC89 < 0) = 0; % avoids negative interp vals
Phanero_fwet_BC89(372:end,1) = Phanero_fwet_PlioPI; % Phanero_fwet_BC89(371,1); % CONSTANT since 5.5 Ma Pliocene as above
Phanero_Awet_BC89 = Phanero_fwet_BC89.*Phanero_landA;
% wetland area (Mkm^2), based on coal:continental sed fraction over time (Berner and Canfield, 1989)
% note that period from 5.5 Ma to PI (0.0 Ma) has assumed constant wetland extent (see above)

% Convert to PI-normalized scaling factors
% normalized to PI value for scaling factor - assumed that not all wetlands preserved in coal basins, so assumed ceteris paribus
Phanero_Awet_Rel_BC89 = Phanero_Awet_BC89./Phanero_Awet_BC89(end); % 

Phanero_Awet_Rel_lo = 1 + 0.75.*(Phanero_Awet_Rel_BC89-1); % dampened variability from PI
Phanero_Awet_Rel_lo(Phanero_Awet_Rel_lo < 0) = 0;
Phanero_Awet_Rel_hi = 1 + 1.25.*(Phanero_Awet_Rel_BC89-1); % enhanced variability from PI
Phanero_Awet_Rel_hi(Phanero_Awet_Rel_hi < 0) = 0;

% global mean surface T (GMST, deg C)
GMST_C_hires = interp1(GMST_t,GMST_C,timeslices,'pchip');
GMST_C16 =  interp1(GMST_traw,GMST_C16raw,timeslices,'pchip');
GMST_C84 =  interp1(GMST_traw,GMST_C84raw,timeslices,'pchip');
% Also equatorial T
Teq_C_hires = interp1(PhanTeq_t,PhanTeq_C,timeslices,'pchip');

% 0.21 atm = 0.212730e6 dynes/cm^2, matches rescaling in photochem scripts
pO2_hires = interp1(pO2Phantime,pO2_atm,timeslices,'pchip')./0.21; % PAL O2; 
pO2_min = interp1(pO2Phantime,pO2_minraw,timeslices,'pchip')./0.21; % PAL O2
pO2_max = interp1(pO2Phantime,pO2_maxraw,timeslices,'pchip')./0.21; % PAL O2

timeUVb_hires = 1e-9.*timeslices; % rescaled to Ga from years

% compute emissions at 1 Ma resolution with interpolated forcings

%gamma_T_CH4_hires_min = Q10_CH2O_photo.^((GMST_C_hires - T_ref_CH4)./10); % typical metabolic Q10 formulation for 
% photosynthesis/CH2O fixation by Rubisco/respiration and decay, for substrate-limited CH4genesis with 
% plant T-response controlling substrate availability across climate changes 

% SENSITIVITY TEST gammaT plot only - NOT USED 
gamma_T_CH4_hires_min = Q10_CH4.^((GMST_C_hires - T_ref_CH4)./10); % Q10 formulation per Zhu+2014
gamma_T_CH4_Beerling = exp(0.1678.*GMST_C_hires)./exp(0.1678.*T_ref_CH4); % PI self-normalized function from Beerling

%  ONLY USE Conrad 2023 formulation
% however, Ea values vary between max and min of ranges in Conrad (79-122
% kJ/mol, 100 kJ/mol as nominal)
gamma_T_CH4_hires_std = exp((-Ea_CH4./Constant.R).*((1./(273.15+GMST_C_hires)) - (1./(273.15+T_ref_CH4)))); % Arrhenius function per Conrad 2023 - nominal GMST history

gamma_T_CH4_hires_stdhiT = exp((-Ea_CH4./Constant.R).*((1./(273.15+GMST_C84)) - (1./(273.15+T_ref_CH4)))); % 
gamma_T_CH4_hires_stdloT = exp((-Ea_CH4./Constant.R).*((1./(273.15+GMST_C16)) - (1./(273.15+T_ref_CH4)))); 
gamma_T_CH4_hires_maxEahiT = exp((-Ea_maxCH4./Constant.R).*((1./(273.15+GMST_C84)) - (1./(273.15+T_ref_CH4)))); % Arrhenius function per Conrad 2023, endmember T and Ea
gamma_T_CH4_hires_minEaloT = exp((-Ea_minCH4./Constant.R).*((1./(273.15+GMST_C16)) - (1./(273.15+T_ref_CH4)))); % 
gamma_T_CH4_hires_minEa = exp((-Ea_minCH4./Constant.R).*((1./(273.15+GMST_C_hires)) - (1./(273.15+T_ref_CH4)))); % lowest/highest bound of mean Ea values for 3 ecosystem scales (culture, sediment, and field - Conrad 2023)
gamma_T_CH4_hires_maxEa = exp((-Ea_maxCH4./Constant.R).*((1./(273.15+GMST_C_hires)) - (1./(273.15+T_ref_CH4)))); % but with nominal GMST history

gamma_T_CH4_hires_Ferm = exp((-Ea_fermentation./Constant.R).*((1./(273.15+GMST_C_hires)) - (1./(273.15+T_ref_CH4)))); % Arrhenius function per Weston and Joye 2005 for fermentation
gamma_T_CH4_hires_FermhiT = exp((-Ea_fermentation./Constant.R).*((1./(273.15+GMST_C84)) - (1./(273.15+T_ref_CH4)))); % 
gamma_T_CH4_hires_FermloT = exp((-Ea_fermentation./Constant.R).*((1./(273.15+GMST_C16)) - (1./(273.15+T_ref_CH4)))); 


% CH4 fluxes (mol CH4/yr) - includes geological, insect, marine, and other sources as constant over deep time

% NOMINAL (std) and max/min envelopes
Flux_CH4_emissions_BC89_hires_std = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_std + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_min = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_lo.*gamma_T_CH4_hires_minEaloT + PI_otherCH4Emissions; % Ea_CH4, Gamma_coal, and T forcing all change vs. nominal
Flux_CH4_emissions_BC89_hires_max = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_hi.*gamma_T_CH4_hires_maxEahiT + PI_otherCH4Emissions; 

% boundary conditions sensitivity tests (single param changed in each - T, Ea, Gamma_coal)
Flux_CH4_emissions_BC89_hires_minT = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_stdloT + PI_otherCH4Emissions; % only T forcing change vs. nominal
Flux_CH4_emissions_BC89_hires_maxT = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_stdhiT + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_minEa = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_minEa + PI_otherCH4Emissions; % only Ea change vs. nominal
Flux_CH4_emissions_BC89_hires_maxEa = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_maxEa + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_minCoal = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_lo.*gamma_T_CH4_hires_std + PI_otherCH4Emissions; % only Gamma_coal changes vs. nominal
Flux_CH4_emissions_BC89_hires_maxCoal = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_hi.*gamma_T_CH4_hires_std + PI_otherCH4Emissions; 

% no-Gamma_coal sensitivity test
Flux_CH4_emissions_BC89_hires_stdNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_std + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_minNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_minEaloT + PI_otherCH4Emissions; % Ea_CH4, Gamma_coal, and T forcing all change vs. nominal
Flux_CH4_emissions_BC89_hires_maxNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_maxEahiT + PI_otherCH4Emissions;

% fermentation-limited sensitivity test, with and without Gamma_coal
Flux_CH4_emissions_BC89_hires_Ferm = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_BC89.*gamma_T_CH4_hires_Ferm + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_minFerm = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_lo.*gamma_T_CH4_hires_FermloT + PI_otherCH4Emissions; % Ea_CH4, Gamma_coal, and T forcing all change vs. nominal
Flux_CH4_emissions_BC89_hires_maxFerm = PIWetlandRefEmissionCH4.*Phanero_Awet_Rel_hi.*gamma_T_CH4_hires_FermhiT + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_FermNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_Ferm + PI_otherCH4Emissions; 
Flux_CH4_emissions_BC89_hires_minFermNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_FermloT + PI_otherCH4Emissions; % Ea_CH4, Gamma_coal, and T forcing all change vs. nominal
Flux_CH4_emissions_BC89_hires_maxFermNoCoal = PIWetlandRefEmissionCH4.*1.*gamma_T_CH4_hires_FermhiT + PI_otherCH4Emissions; 


% Normalized to reference fluxes in PHOTOCHEM
FCH4_tot_hires_BC89_std = Flux_CH4_emissions_BC89_hires_std./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_min = Flux_CH4_emissions_BC89_hires_min./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_max = Flux_CH4_emissions_BC89_hires_max./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr

FCH4_tot_hires_BC89_minT = Flux_CH4_emissions_BC89_hires_minT./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxT = Flux_CH4_emissions_BC89_hires_maxT./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_minEa = Flux_CH4_emissions_BC89_hires_minEa./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxEa = Flux_CH4_emissions_BC89_hires_maxEa./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_minCoal = Flux_CH4_emissions_BC89_hires_minCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxCoal = Flux_CH4_emissions_BC89_hires_maxCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr

FCH4_tot_hires_BC89_stdNoCoal = Flux_CH4_emissions_BC89_hires_stdNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_minNoCoal = Flux_CH4_emissions_BC89_hires_minNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxNoCoal = Flux_CH4_emissions_BC89_hires_maxNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr

FCH4_tot_hires_BC89_Ferm = Flux_CH4_emissions_BC89_hires_Ferm./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_minFerm = Flux_CH4_emissions_BC89_hires_minFerm./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxFerm = Flux_CH4_emissions_BC89_hires_maxFerm./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_FermNoCoal = Flux_CH4_emissions_BC89_hires_FermNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_minFermNoCoal = Flux_CH4_emissions_BC89_hires_minFermNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
FCH4_tot_hires_BC89_maxFermNoCoal = Flux_CH4_emissions_BC89_hires_maxFermNoCoal./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr



%% %% Pleistocene Ice-Core pCH4 vs. modeled T-dependence comparison
timeslices2 = flip(0:1e3:8e5).'; % every ka
%timeslices = cat(1,timeslice1,timeslice2); 

% added constant PIM and mid-Pleistocene     (14.6./21.4)
% values (assumedly equal, Pliocene coal deposition rescaled by PIM/Pliocene wetland area estimates from Beerling+2009 
% Table 2) for comparison against ice core pCH4 record over past 800 ka
% MODIFIED to achieve more reasonable pCH4 at PIM (possibly counteracting
% overestimated Tau_CH4 of 13.375 years (could have been more like 9 years, hence overestimated pCH4)
% TUNED!!! Retuned PIM/Pleistocene:Pliocene wetland areal ratio >50% lower
% than suggested by Beerling+2009, but retrieves PIM emissions which are
% not too much lower than lower bound of prior estimates summarized in
% Table 2
%CoalDepoRateHiRes = (0.85.*9./21.4); % now normalized to mid-Pliocene (3.5 Ma) value of 0.020
% Note that PIM is not included for these "anchor point" calculations

%CoalDepoRate_Pleisto_low = 1.65.*(0.85.*9./21.4); % 1.65 high (~715 ppb), 1.56 mid (660 ppb), 1.4 low (560 ppb) tuned to retrieve reasonable PIM pCH4 ~560 ppb for 1.4 or ~660-680 for 1.56 (per Ruddiman et al. 2015) in PI-tuned PHOTOCHEM model
% interp1(CoalDepoTimeHiRes,CoalDepoRateHiRes,timeslices,'pchip');
%CoalDepoRate_Pleisto_high = CoalDepoRate_Pleisto_low; %1.56.*(0.85.*9./21.4); %CoalDepoRate_Pleisto_low; % tuned to retrieve reasonable PIM pCH4 ~560 ppb (per Ruddiman et al. 2015) in PI-tuned PHOTOCHEM model

PI_globalCH4fluxHi = 16e12 - PI_otherCH4Emissions; % 16e12 = ~715 ppb PI modern (near last deglacial peak); 15e12 = ~660 ppb (mid-range); 13.28e12 = ~560 ppb (Ruddiman+2015)
PI_globalCH4fluxLo = PI_globalCH4fluxHi; % 16e12 - PI_otherCH4Emissions; % 13.28e12 16e12 = ~715 ppb PI modern (near last deglacial peak); 15e12 = ~660 ppb (mid-range); 13.28e12 = ~560 ppb (Ruddiman+2015)

GMST_C_Pleisto = interp1(GMST_t,GMST_C,timeslices2,'pchip');
% regression best-fit from Judd+2024 outputs to estimate equatorial/tropical mean SAT from Snyder+2016 GMST
Teq_C_Pleisto = 0.0146.*(GMST_C_Pleisto).^2 + 0.1612.*GMST_C_Pleisto + 19.624; % R^2 = 0.9311, polynomial best-fit to Tequatorial vs GMST over last 375 Ma

pO2_Pleisto = interp1(pO2Phantime,pO2_atm,timeslices2,'pchip')./0.21; % PAL O2

timeUVb_Pleisto = 1e-9.*timeslices2; % rescaled to 

T_ref_CH4PIM = 14; % 14 C PIM reference T, same as for deep time currently

% compute emissions at 1 Ma resolution with interpolated forcings
%gamma_T_CH4_Pleisto_min = Q10_CH2O_photo.^((GMST_C_Pleisto - T_ref_CH4PIM)./10); % typical metabolic Q10 formulation for respiration/photosynthesis/C-fixation
%gamma_T_CH4_Pleisto_mid = Q10_CH4.^((GMST_C_Pleisto - T_ref_CH4PIM)./10); % Q10 formulation per Zhu+2014
gamma_T_CH4_Pleisto_max = exp((-Ea_CH4./Constant.R).*((1./(273.15+GMST_C_Pleisto)) - (1./(273.15+T_ref_CH4PIM))));

% compute fluxes assuming no wetland area change over last 800 ka, constant other CH4 sources as above (attenuates T-variation in fluxes) 
Flux_CH4_emissions_tot_Pleisto_max = PI_globalCH4fluxHi.*gamma_T_CH4_Pleisto_max + PI_otherCH4Emissions; % PlioRefEmissionCH4.*CoalDepoRate_Pleisto_high
%Flux_CH4_emissions_tot_Pleisto_mid = PI_globalCH4fluxLo.*gamma_T_CH4_Pleisto_mid + PI_otherCH4Emissions; % PlioRefEmissionCH4.*CoalDepoRate_Pleisto_low
% Flux_CH4_emissions_tot_Pleisto_min = PI_globalCH4fluxLo.*gamma_T_CH4_Pleisto_min + PI_otherCH4Emissions; % PlioRefEmissionCH4.*CoalDepoRate_Pleisto_low

% rescale fluxes for input into PHOTOCHEM
FCH4_tot_Pleisto_max = Flux_CH4_emissions_tot_Pleisto_max./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
%FCH4_tot_Pleisto_mid = Flux_CH4_emissions_tot_Pleisto_mid./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr
% FCH4_tot_Pleisto_min = Flux_CH4_emissions_tot_Pleisto_min./(0.4.*1e11.*photochem_fluxscaling); % normalized to reference flux of roughly 171.3 Tg CH4/yr


%% %% Save Output to Run in *PHOTOCHEM*

timeUVb_Pleisto5ka = timeUVb_Pleisto(1:5:length(timeslices2));
GMST_C_Pleisto5ka = GMST_C_Pleisto(1:5:length(timeslices2));
pO2_Pleisto5ka = pO2_Pleisto(1:5:length(timeslices2));
%FCH4_tot_Pleisto_min5ka = FCH4_tot_Pleisto_min(1:5:length(timeslices2));
%FCH4_tot_Pleisto_mid5ka = FCH4_tot_Pleisto_mid(1:5:length(timeslices2));
FCH4_tot_Pleisto_max5ka = FCH4_tot_Pleisto_max(1:5:length(timeslices2));
Teq_C_Pleisto5ka = Teq_C_Pleisto(1:5:length(timeslices2));

save("time_GMST_pO2_FCH4_photochemInputHiRes_updated.mat",'timeUVb_hires', 'GMST_C_hires', 'GMST_C16', 'GMST_C84', 'pO2_hires', 'pO2_max', 'pO2_min', 'FCH4_tot_hires_BC89_std',...
    'FCH4_tot_hires_BC89_min', 'FCH4_tot_hires_BC89_max', 'FCH4_tot_hires_BC89_minCoal', 'FCH4_tot_hires_BC89_maxCoal', 'FCH4_tot_hires_BC89_minT', 'FCH4_tot_hires_BC89_maxT',...
    'FCH4_tot_hires_BC89_minEa', 'FCH4_tot_hires_BC89_maxEa', 'time_CH4FluxB09','FluxCH4_B09','FCH4_tot_hires_BC89_stdNoCoal','FCH4_tot_hires_BC89_minNoCoal','FCH4_tot_hires_BC89_maxNoCoal',...
    'FCH4_tot_hires_BC89_Ferm','FCH4_tot_hires_BC89_minFerm','FCH4_tot_hires_BC89_maxFerm','FCH4_tot_hires_BC89_FermNoCoal','FCH4_tot_hires_BC89_minFermNoCoal','FCH4_tot_hires_BC89_maxFermNoCoal','Teq_C_hires'); % 
% OR
save("time_GMST_pO2_FCH4_photochemInputPleisto.mat",'timeUVb_Pleisto5ka', 'GMST_C_Pleisto5ka', 'pO2_Pleisto5ka', 'FCH4_tot_Pleisto_max5ka','Teq_C_Pleisto5ka'); % 

%% Output Figures

polyx = cat(1,timeslices./1e6,flip(timeslices./1e6)); % for polyshape ranges
polyT = cat(1,GMST_C84,flip(GMST_C16));%.';
polyO2 = cat(1,pO2_max,flip(pO2_min));%.';
polyAcoal = cat(1,Phanero_Awet_Rel_hi,flip(Phanero_Awet_Rel_lo));%.';
polygT = cat(1,gamma_T_CH4_hires_maxEa,flip(gamma_T_CH4_hires_minEa));%.';

polyconst = [33.9, 33.9, 149.24, 149.24]; % constant, so can plot as rectangle
polyKPg = cat(1,-5, 48, 48, -5).'; % constant, so can plot as rectangle
polyKPgrf = cat(1,-20, 1000, 1000, -20).'; % constant, so can plot as rectangle
polyKPglog = cat(1,0.00001, 1000, 1000, 0.00001).'; % constant, so can plot as rectangle

figure(101);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(101),'defaultAxesColorOrder',[lc; rc]);
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on
%plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
%hold on 

plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'k','Marker','d','LineWidth',1.25)
%hold on
plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_std./1e12,'m','Marker','^','LineWidth',2,'Linestyle','-','MarkerSize',3) % 'none' for no marker
hold on
%
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':')
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_max./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':') % 'none' for no marker
%hold on
%
scatter(-5,1.31359e13./1e12,'o','MarkerFaceColor','k')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on

xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('CH_4 Emissions (Tmol CH_4/yr)')
%ylim([-80,1400]) % [-35,650]
ylim([-35,650]) % [-35,650]
title('Global CH_4 Emissions')
pbaspect([2 1 1])
fontsize(24,"points") % 14   'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T per Rubisco, GMST per J+24)',  'Weak \gamma_T per Zhu+(2014) & Judd+(2024)',
L = legend('CH_4 Emissions per Beerling+(2009)',...
    'Revised Phanerozoic CH_4 Emissions',...
    'FontSize',24); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
L.AutoUpdate = 'off';


yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(24,"points") % 14
box on
%title("Coal:Continental Clastics Wetland Area Forcing (Berner+Canfield 1989)")



figure(1022);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(1022),'defaultAxesColorOrder',[lc; rc]);
tiledlayout(1,1,"TileSpacing","compact","Padding","compact")

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_min,'k','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_min,'b-','Marker','v','LineWidth',1.25)
%hold on
plot(timeslices./1e6,gamma_T_CH4_hires_std,'m-','Marker','^','LineWidth',2, 'MarkerSize',3) % 'none' for no marker
hold on
plot(timeslices./1e6,gamma_T_CH4_hires_maxEa,'m:','Marker','none','LineWidth',2) % 'none' for no marker
hold on
plot(timeslices./1e6,gamma_T_CH4_hires_minEa,'m:','Marker','none','LineWidth',2) % 'none' for no marker
hold on
%plot(timeslices./1e6,gamma_T_CH4_Beerling,'k-','Marker','o','LineWidth',1.25) % 'none' for no marker
%hold on
scatter(-5,1,'o','MarkerFaceColor','k')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
title('Phanerozoic GMST-Sensitivity (\gamma_T) Function for Coal Wetland CH_4 Emissions')
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('T-Sensitivity Scaling Factor (\gamma_T)')
%ylim([-1.2,21]) % -0.85,22
ylim([-2,40]) 
pbaspect([2 1 1])
fontsize(24,"points") % 14 'Weak Phanerozoic CH_4 Emission T-Sensitivity (Rubisco-limited, Q_{10}= 2)',
%L = legend('Weak CH_4 Emission \gamma_T (Q_{10}= 2.48 per Zhu+2014)','Strong CH_4 Emission \gamma_T (E_a = 100 kJ/mol per Conrad 2023)','Boreal CH_4 Emission \gamma_T (per Beerling+2009)',...
%    'FontSize',24); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';


yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(24,"points") % 14

box on




figure(1033);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(1033),'defaultAxesColorOrder',[lc; rc]);


tiledlayout(3,1,"TileSpacing","compact","Padding","compact")
%subplot(3,1,1)

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(timeslices./1e6,GMST_C_hires,'r','LineStyle','-','LineWidth',1.25)
hold on
plot(timeslices./1e6,GMST_C16,'r','LineStyle',':','LineWidth',1.25)
hold on
plot(timeslices./1e6,GMST_C84,'r','LineStyle',':','LineWidth',1.25)
hold on
scatter(-5,14,'o','MarkerFaceColor','r')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
annotation('textbox',[.28 .85-0.011 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)')
ylim([6,45])
pbaspect([3 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Global Mean Surface Temperature (GMST)')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on


nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices/1e6,pO2_hires,'b','LineStyle','-','LineWidth',1.25)
hold on
plot(timeslices/1e6,pO2_min,'b','LineStyle',':','LineWidth',1.25)
hold on
plot(timeslices/1e6,pO2_max,'b','LineStyle',':','LineWidth',1.25)
hold on
scatter(-5,1,'o','MarkerFaceColor','b')
hold on
annotation('textbox',[.28 .55-0.033 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on

xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('pO_2 (PAL)')
ylim([0.6,2.2])
pbaspect([3 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Atmospheric O_2 Level (pO_2)')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices/1e6,Phanero_Awet_Rel_BC89,'Color',darkgreen,'LineStyle','-','LineWidth',1.25)
hold on
%
plot(timeslices/1e6,Phanero_Awet_Rel_lo,'Color',darkgreen,'LineStyle',':','LineWidth',1.25)
hold on
plot(timeslices/1e6,Phanero_Awet_Rel_hi,'Color',darkgreen,'LineStyle',':','LineWidth',1.25)
hold on
%
scatter(-5,1,'o','MarkerFaceColor',darkgreen)
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
annotation('textbox',[.28 .25-0.033-0.023 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); 
ylabel('\Gamma_{coal}')
ylim([-0.5,7.5]) % [-0.35,6]
pbaspect([3 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Proxy for Relative Global "Coal Wetland" Area (\Gamma_{coal})')

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on



figure(601);clf

size = 45;

scatter(-5,1.31359e13./1e12,'square','MarkerFaceColor','k')
hold on
%plot(timeslices2./1e3,Flux_CH4_emissions_tot_Pleisto_min./1e12,'-k','Marker','none')
%hold on
plot(timeslices2./1e3,Flux_CH4_emissions_tot_Pleisto_max./1e12,'-m','Marker','none') % 'none' for no marker
hold on
%plot(timeslices2./1e3,Flux_CH4_emissions_tot_Pleisto_mid./1e12,'-b','Marker','none') % 'none' for no marker
%hold on

%plot(timeslices./1e3,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','none')
%hold on
errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
hold on

xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,800]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (ka)'); ylabel('Wetland Emissions of Methane (Tmol CH_4/yr)')
%ylim([-60,970])
pbaspect([3 1 1])
fontsize(18,"points") % 14 ,'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T)',
L = legend('Preindustrial Modern CH_4 Emissions (cf. B+09)',...
    'Revised Phanerozoic CH_4 Emissions (Strong \gamma_T)','FontSize',12); % 14,'Modern pN_2O (337 ppb)'
L.AutoUpdate = 'off';

box on









figure(10335);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(10335),'defaultAxesColorOrder',[lc; rc]);


tiledlayout(4,1,"TileSpacing","compact","Padding","compact")
%subplot(3,1,1)

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(polyshape(polyx,polyT),'FaceColor','k','EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices./1e6,GMST_C_hires,'k','LineStyle','-','LineWidth',1.25)
hold on
%plot(timeslices./1e6,GMST_C16,'r','LineStyle',':','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,GMST_C84,'r','LineStyle',':','LineWidth',1.25)
%hold on
scatter(-5,14,'o','MarkerFaceColor','k')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
annotation('textbox',[.28 .85-0.011 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('GMST (^\circC)')
ylim([6,45])
pbaspect([4 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Global Mean Surface Temperature (GMST)')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on


nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(polyshape(polyx,polyO2),'FaceColor','b','EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices/1e6,pO2_hires,'b','LineStyle','-','LineWidth',1.25)
hold on
%plot(timeslices/1e6,pO2_min,'b','LineStyle',':','LineWidth',1.25)
%hold on
%plot(timeslices/1e6,pO2_max,'b','LineStyle',':','LineWidth',1.25)
%hold on
scatter(-5,1,'o','MarkerFaceColor','b')
hold on
annotation('textbox',[.28 .598 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on

xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('pO_2 (PAL)')
ylim([0.6,2.2])
pbaspect([4 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Atmospheric O_2 Level (pO_2)')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(polyshape(polyx,polyAcoal),'FaceColor',darkgreen,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices/1e6,Phanero_Awet_Rel_BC89,'Color',darkgreen,'LineStyle','-','LineWidth',1.25)
hold on
plot(polyshape([48, 48, 56, 56],cat(1,2, 2.5, 2.5, 2).'),'FaceColor','g','EdgeColor','none'); % cf. Wilton+2019 abstract, Ypresian Eocene (56-48 Ma) wetland fraction is 2-2.5x higher than reference modern value 
% (wetland area estimate somewhat outdated cf. Hopcroft+2020, but mainly we care about relative change between
% Eocene and modern, not absolute values - see Table 2 in Wilton+2019
hold on 
%
%plot(timeslices/1e6,Phanero_Awet_Rel_lo,'Color',darkgreen,'LineStyle',':','LineWidth',1.25)
%hold on
%plot(timeslices/1e6,Phanero_Awet_Rel_hi,'Color',darkgreen,'LineStyle',':','LineWidth',1.25)
%hold on
%
scatter(-5,1,'o','MarkerFaceColor',darkgreen)
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
annotation('textbox',[.28 .355 .1 .1],'String','C','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
ylabel('\Gamma_{coal}')
ylim([-0.5,7.5]) % [-0.35,6]
pbaspect([4 1 1])
fontsize(12,"points") % 14
%L = legend('Global Mean Surface Temperature (GMST)',...
%    'FontSize',12); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';
title('Proxy for Relative Global "Coal Wetland" Area (\Gamma_{coal})')

yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on

nexttile

yyaxis left

size = 45;

%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_lo./1e12,size,'o','MarkerFaceColor','b')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_hi./1e12,size,'^','MarkerFaceColor','m')
%hold on
%scatter(AnchorTime./1e6,Flux_CH4_emissions_tot_noT./1e12,size,'d','MarkerFaceColor','g')
%hold on

%plot(CH4fluxtime./1e6,CH4flux./1e12,'k','Marker','o','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','square','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_tot_OG./1e12,'r','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_min,'k','Marker','d','LineWidth',1.25)
%hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_min,'b-','Marker','v','LineWidth',1.25)
%hold on
plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(polyshape(polyx,polygT),'FaceColor','m','EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 
plot(timeslices./1e6,gamma_T_CH4_hires_std,'m-','Marker','none','LineWidth',1.25)%, 'MarkerSize',3) % 'none' for no marker
hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_maxEa,'m:','Marker','none','LineWidth',1.25) % 'none' for no marker
%hold on
%plot(timeslices./1e6,gamma_T_CH4_hires_minEa,'m:','Marker','none','LineWidth',1.25) % 'none' for no marker
%hold on
%plot(timeslices./1e6,gamma_T_CH4_Beerling,'k-','Marker','o','LineWidth',1.25) % 'none' for no marker
%hold on
scatter(-5,1,'o','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on
annotation('textbox',[.28 .114 .1 .1],'String','D','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on
title('GMST-Sensitivity (\gamma_T) for Coal Wetland CH_4 Emissions')
xline([0],'-k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('T-Sensitivity Factor (\gamma_T)')
%ylim([-1.2,21]) % -0.85,22
ylim([-2,40]) 
pbaspect([4 1 1])
fontsize(12,"points") % 14 'Weak Phanerozoic CH_4 Emission T-Sensitivity (Rubisco-limited, Q_{10}= 2)',
%L = legend('Weak CH_4 Emission \gamma_T (Q_{10}= 2.48 per Zhu+2014)','Strong CH_4 Emission \gamma_T (E_a = 100 kJ/mol per Conrad 2023)','Boreal CH_4 Emission \gamma_T (per Beerling+2009)',...
%    'FontSize',24); % 14,'Modern pN_2O (337 ppb)'     
% 'Revised Phanerozoic CH_4 Emissions (no \gamma_T)','Revised Phanerozoic CH_4 Emissions (\gamma_T per B+09, GMST per J+24)',...
%L.AutoUpdate = 'off';


yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(12,"points") % 14

box on





figure(6011);clf

size = 45;

%scatter(-5,1.31359e13./1e12,'square','MarkerFaceColor','k')
%hold on
%plot(timeslices2./1e3,Flux_CH4_emissions_tot_Pleisto_min./1e12,'-k','Marker','none')
%hold on
plot(sort(GMST_C_hires),sort(gamma_T_CH4_hires_std),'-m','Marker','none','LineWidth',3) % 'none' for no marker
hold on
plot(sort(GMST_C_hires),sort(gamma_T_CH4_Beerling),'--k','Marker','none','LineWidth',3) % 'none' for no marker
hold on

plot(sort(GMST_C_hires),sort(gamma_T_CH4_hires_maxEa),'m:','LineWidth',2) % 'none' for no marker
hold on
plot(sort(GMST_C_hires),sort(gamma_T_CH4_hires_minEa),'m:','LineWidth',2) % 'none' for no marker
hold on
%plot(timeslices./1e3,Flux_CH4_emissions_tot_hires_noT./1e12,'g','Marker','none')
%hold on
%errorbar(-5,10.162,NaN,0.92846,'Color','k','CapSize',1)
%hold on

xline([14],'--k')
hold on
yline([1],'--k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
%set(gca,'XDir','reverse');
%set(gca,'xlim',[-10,800]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('GMST (^\circC)'); ylabel('T-Sensitivity Scaling Factor (\gamma_T)')
xlim([12,36])
%pbaspect([3 1 1])
fontsize(18,"points") % 14 ,'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T)',
L = legend('This Study','Beerling+(2009)',...
    'FontSize',18); % 14,'Modern pN_2O (337 ppb)'
L.AutoUpdate = 'off';

box on






figure(6012);clf
lc = [0 0 0];
rc = [0 0 0];
set(figure(6012),'defaultAxesColorOrder',[lc; rc]);


tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
%subplot(3,1,1)

nexttile

yyaxis left
size = 45;

plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
hold on 

plot(timeslices./1e6,Phanero_landA./Phanero_landA(end),'c-','Marker','none','LineWidth',1.25) % 'none' for no marker
hold on
plot(timeslices./1e6,Phanero_fwet_BC89./Phanero_fwet_BC89(end),'g-','Marker','none','LineWidth',1.25) % 'none' for no marker
hold on
plot(timeslices./1e6,Phanero_Awet_Rel_BC89,'color',darkgreen,'Marker','none','LineWidth',2) % 'none' for no marker
hold on
%plot(In_data.PhanBiomes.LiHumidtime, In_data.PhanBiomes.LiHumidPercent./In_data.PhanBiomes.LiHumidPercent(end),'g-.','Marker','none','LineWidth',1.25)
scatter(-5,1,'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','o') % 'none' for no marker
hold on
xline([0],'-k')
hold on
plot(polyshape([48, 48, 56, 56],cat(1,2, 2.5, 2.5, 2).'),'FaceColor',darkgreen,'EdgeColor','none'); % cf. Wilton+2019 abstract, Ypresian Eocene (56-48 Ma) wetland fraction is 2-2.5x higher than reference modern value 
% (wetland area estimate somewhat outdated cf. Hopcroft+2020, but mainly we care about relative change between
% Eocene and modern, not absolute values - see Table 2 in Wilton+2019
hold on 
% lower bounds on wetland area from coal forest area estimates from Cleal and Thomas 2005,
% assumes 4 Mkm^2 wetland area per Wilton+2019, though note much higher estimate in Hopcroft+2020 and refs 
% (would make these lower bounds even low and harder to evaluate)
% ages are taken as midpoints (if "middle" or "-"), beginning points (if
% "early"), or endpoints (if "late") for each stage of the Carbon-Permian
% plot(323.4, 467./4000,'^','Color',darkgreen) 
% hold on
% plot(315.2, 1786./4000,'^','Color',darkgreen) 
% hold on
% plot(311.1, 1721./4000,'^','Color',darkgreen) 
% hold on
% plot(307, 2395./4000,'^','Color',darkgreen) 
% hold on
% plot(305.35, 1131./4000,'^','Color',darkgreen) 
% hold on
% plot(301.3, 1087./4000,'^','Color',darkgreen) 
% hold on
% plot(296.21, 1590./4000,'^','Color',darkgreen) 
% hold on
% plot(290.1, 1690./4000,'^','Color',darkgreen) 
% hold on
% plot(278.85, 105./4000,'^','Color',darkgreen) 
% hold on
% plot(256.825, 395./4000,'^','Color',darkgreen) 
% hold on
% plot(253.021, 140./4000,'^','Color',darkgreen) 
% hold on
%hold on
%yline([1],'--k')
%plot(4.5-(time./1e9),PGC.pH.surface,'color','r','LineStyle','-')
%set(gca,'XDir','reverse');
%set(gca,'xlim',[-10,800]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); ylabel('Relative Units (normalized to PI)')
set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375],'ylim',[-0.45,8]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%xlim([12,36])
%pbaspect([3 1 1])
fontsize(18,"points") % 14 ,'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T)',
L = legend('','A_{land} (Global Land Area)','f_{coal} (Coal Wetland Fraction of Land Area)','\Gamma_{coal} (Global Coal Wetland Area)',...
    'FontSize',18); % 14,'Modern pN_2O (337 ppb)'
L.AutoUpdate = 'off';

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[-10,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(24,"points") % 14

box on



figure(601211);clf
lc = [0 0 0];
rc = [0 0 0];
set(figure(601211),'defaultAxesColorOrder',[lc; rc]);


tiledlayout(1,1,"TileSpacing","compact","Padding","compact")
%subplot(3,1,1)

nexttile

yyaxis left
size = 45;

%plot(polyshape(polyconst,polyKPg),'FaceColor',KPgcolor,'EdgeColor','none'); % PRIOR - Judd et al. 2024 GMST (nominal, 50%ile)
%hold on 

plot(In_data.PhanBiomes.LiHumidtime, In_data.PhanBiomes.LiHumidPercent./In_data.PhanBiomes.LiHumidPercent(end),'g-','Marker','none','LineWidth',2)

xlabel('Age Before Present (Ma)'); ylabel('Relative Humid Land Area (normalized to PI)')
set(gca,'XDir','reverse');
set(gca,'xlim',[0,250]) % ,'xtick',time_ticks,'ylim',[0,1.5e5]
%xlim([12,36])
%pbaspect([3 1 1])
fontsize(18,"points") % 14 ,'Revised Phanerozoic CH_4 Emissions (Weak \gamma_T)',
%L = legend('','A_{land} (Global Land Area)','f_{coal} (Coal Wetland Fraction of Land Area)','\Gamma_{coal} (Global Coal Wetland Area)',...
%    'FontSize',18); % 14,'Modern pN_2O (337 ppb)'
%L.AutoUpdate = 'off';
title('Humid Land Area Fraction (Li et al. 2025)')

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);
yticks([]);

set(gca,'XDir','reverse');
set(gca,'xlim',[0,250]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
fontsize(24,"points") % 14

box on
