%% Plot Comparisons of PHOTOCHEM vs. Beerling+2009 Cambridge 2-D CTM
%
%
% %

Constant.mia = 1.773e20; %5e18./0.029, original EarthN formula;

%% Inputs
% Primary input structure
In_data = load('biome_output.mat');

BeerlingpCH4 = readtable("All_CSV_files/B+09_pCH4_ppb_vs_Ma_constpO2_WPD.csv"); % WPD of original Figure 3a, constant 21% pO2
BeerlingpCH4_ppb_350Ma = cat(1,0,0,0,BeerlingpCH4{:,2}); % roughly 10 Ma intervals from 350 Ma to 0 Ma, adds empty zeroes back to 380 Ma

%BeerlingpOH = readtable("All_CSV_files/B+09_pOH_10^6cm3_vs_Ma_constpO2_WPD.csv"); % WPD of Fig 3b, constant 21% pO2
%BeerlingpOH_380Ma = BeerlingpOH{3:end,2}; % roughly 10 Ma intervals from 380 Ma to 0 Ma

TotalTime = flip(0:10e6:380e6).';

% CH4 lifetime from B+09, estimated
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
FluxCH4_B09 = CH4flux(2:end)./(0.4.*26.7e12); % converts into multiples of base-flux
% replace weird PI output with rescaled value from B+09 Table 2 for PI
FluxCH4_B09(end) = 1.23; % 1.22995414 = 1.23 * (0.4.*26.7e12) mol CH4/yr from wetlands global

TauCH4_B09 = 1e-9.*BeerlingpCH4_ppb_350Ma.*Constant.mia./(FluxCH4_B09.*(0.4.*26.7e12)); 
% computes steady-state emission lifetime of CH4 by dividing estimated CH4 burden by global wetland emissions


% Fig 2a from Holmes 2018
GEOS_blue = readtable("All_CSV_files\Holmes2018_fig2a_GEOS-Chem_blue.csv");
GEOS_blue_burden = GEOS_blue{:,1};
GEOS_blue_Tau = GEOS_blue{:,2};

Oslo_green = readtable("All_CSV_files\Holmes2018_fig2a_Oslo_CTM2_green.csv");
Oslo_green_burden = Oslo_green{:,1};
Oslo_green_Tau = Oslo_green{:,2};

% Holmes 2018 GEOS-Chem anchor points (approximate)
HolmesBurdensCH4 = [1118.36; 2236.93; 2983.88; 4474.58; 5367.2; 6710.22; 8951.4; 17909.35; 35836.68; 71719.3];

Holmes_pCH4 = HolmesBurdensCH4.*1e12.*1.013e6./(16.04.*Constant.mia); % converts to mixing ratio, then to dynes/cm^2 for application in PHOTOCHEM

% photochem at PI conditions, no dry deposition of CH4
PIBC_std = load('CH4_O3_outputs_B09FCH4_NoVdepSoil14Cinterp_PIbcTest_HolmesSS.mat'); % 14 C interp
pOH_modeled14C = PIBC_std.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,3)./3.7e-14;
pCH4_modeled14C = PIBC_std.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1).*1e9; % ppb
Tau_CH4_modeled14C = PIBC_std.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); % already in years

% Photochem_recalc_GEOS = load("CH4_O3_outputs_Holmes_GEOSChem_anchors_setpCH4_ModernAtmoTest_PIbc.mat");
% Photochem_GEOS_pCH4 = Photochem_recalc_GEOS.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound_O2flux_CH4O2ratio(:,1); % mixing ratio
% Photochem_GEOS_TauCH4 = Photochem_recalc_GEOS.pCH4_TauCH4_pOHtropo_pO3tropo_surf_bound_O2flux_CH4O2ratio(:,2); % years

% PI photochem with CH4 dry deposition at nominal Vdep
PI_Photochem_soilCons = load("CH4_O3_outputs_B09FCH4_effVdepSoil00005_14Cinterp_PIbcTest_SSHolmes.mat");
Photochem_soil_pCH4 = PI_Photochem_soilCons.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,1); % mixing ratio
Photochem_soil_TauCH4 = PI_Photochem_soilCons.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,2); % years

%% Figures


figure(1032);clf

%semilogx(16.04.*1e-9.*pCH4_modeled24C(1:end-1).*Constant.mia./1e12,Tau_CH4_modeled24C(1:end-1),'c','Marker','o','LineStyle','none','MarkerEdgeColor','c','MarkerFaceColor','c') % excluding PI value as outlier/altered
%hold on
semilogx(16.04.*1e-9.*pCH4_modeled14C(1:end-1).*Constant.mia./1e12,Tau_CH4_modeled14C(1:end-1),'k','Marker','*','LineStyle','none','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10) % excluding PI value as outlier/altered
hold on
semilogx(16.04.*Photochem_soil_pCH4(1:end-1).*Constant.mia./1e12,Photochem_soil_TauCH4(1:end-1),'m','Marker','*','MarkerFaceColor','m','LineStyle','none','MarkerSize',10) %
hold on
%semilogx(16.04.*1e-9.*pCH4_modeledModern(1:end-1).*Constant.mia./1e12,Tau_CH4_modeledModern(1:end-1),'r','Marker','^','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','r') % excluding PI value as outlier/altered
%hold on
%semilogx(16.04.*1e-9.*pCH4_modeledModern80(1:end-1).*Constant.mia./1e12,Tau_CH4_modeledModern80(1:end-1),'y','Marker','^','LineStyle','none','MarkerEdgeColor','r','MarkerFaceColor','y') % excluding PI value as outlier/altered
%hold on
semilogx(16.04.*1e-9.*BeerlingpCH4_ppb_350Ma(1:end-1).*Constant.mia./1e12,TauCH4_B09(1:end-1),'g','Marker','d','LineStyle','none','MarkerEdgeColor','g','MarkerFaceColor','b','MarkerSize',10) % excluding PI value as outlier/altered
hold on
semilogx(GEOS_blue_burden,GEOS_blue_Tau,'color',[0 0 0.8],'Marker','none');
hold on
semilogx(Oslo_green_burden,Oslo_green_Tau,'color',[0 0.6 0],'Marker','none');
hold on
%semilogx(GEOS_orange_burden,GEOS_orange_Tau,'color',[254,97,0]./255,'Marker','none');
% Photochem_GEOS_pCH4
%hold on
%semilogx(16.04.*Photochem_GEOS_pCH4.*Constant.mia./1e12,Photochem_GEOS_TauCH4,'color',[0.7 0.6 0],'Marker','square','LineStyle','none','MarkerEdgeColor',[0.7 0.6 0],'MarkerFaceColor',[0.7 0.6 0]);
%hold on
xline(1922,'k-');
hold on
xline(4920,'k--');
hold on
xline(10275,'k:');
hold on
%semilogx(16.04.*2.7e-5.*Constant.mia./1e12,21.2,'g','Marker','*','MarkerFaceColor','g','LineStyle','none') %
%hold on
%semilogx(16.04.*1.12e-6.*Constant.mia./1e12,8.02,'m','Marker','*','MarkerFaceColor','m','LineStyle','none') %
%hold on
%semilogx(16.04.*8.918e-6.*Constant.mia./1e12,14.075,'g','Marker','*','MarkerFaceColor','g')
%hold on
%semilogx(16.04.*2.7e-6.*Constant.mia./1e12,8.52,'g','Marker','*','MarkerFaceColor','g')
%hold on
%semilogx(16.04.*9.137e-7.*Constant.mia./1e12,5.77,'g','Marker','*','MarkerFaceColor','g')
%hold on
%semilogx(16.04.*7.735e-7.*Constant.mia./1e12,5.51,'g','Marker','*','MarkerFaceColor','g')
%hold on
xlim([700, 1e5]);
%semilogx(16.04.*9.957e-6.*Constant.mia./1e12,15.8,'m','Marker','*','MarkerFaceColor','m','LineStyle','none') %
%hold on
%semilogx(16.04.*1.31e-6.*Constant.mia./1e12,8.32,'m','Marker','*','MarkerFaceColor','m','LineStyle','none') %
%hold on
%semilogx(16.04.*3.52e-6.*Constant.mia./1e12,11.18,'m','Marker','*','MarkerFaceColor','m','LineStyle','none') %
%hold on
%plot(Naik_OH,Naik_TauCH4,'k','Marker','square','LineStyle','none','MarkerEdgeC
% olor','k','MarkerFaceColor','k')
%hold on
xlabel('Atmospheric CH_4 Burden (Tg CH_4)'); ylabel('CH_4 Emission Lifetime (\tau_{CH_4}, years)')
L = legend('\itphotochem\rm, no soil sink','\itphotochem\rm, with soil sink',...
    'Calculated from Beerling+2009, Cambridge 2-D CTM',...
    'GEOS-Chem (Holmes 2018)','Oslo CTM2 (Holmes 2018, cf. Isaksen+2011)',...
    '1750 CH_4 Burden (cf. Holmes 2018)','2010 CH_4 Burden (cf. Holmes 2018)','2100 CH_4 Burden (RCP8.5, cf. Holmes 2018)',...
    'FontSize',24); % 'Photochem PI w/ modern dT/dz, 21% pO_2, w/ soil Vdep sink', 'Modern PHOTOCHEM per Wogan+2025', ,'PHOTOCHEM PI tuned, modern dT/dz, 21% pO_2, prescribed pCH_4 per GEOS-Chem',
%,'Revised Phanerozoic pCH_4 (Strong \gamma_T)',...
%    'Revised Phanerozoic pCH_4 (No \gamma_T)' 'GEOS-Chem \tau_p (Holmes 2018)',
L.AutoUpdate = 'off';
%title('Atmospheric CH_4 Lifetime vs. Burden')
% ALL photochem runs computed with CH_4 emissions digitized from
% Beerling+2009 to match B+09 data. EXCLUDES PI results, since skewed in
% some of these older output files (since corrected in Emissions calcs)
fontsize(24,"points") % 14
box on



