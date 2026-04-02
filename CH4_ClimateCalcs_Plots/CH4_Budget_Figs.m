% % % % % % % % % % % % % % % % % %
%% Methane Budget Figure
%
%
%
%
%
%
% % % % % % % % % % % % % % % % % %

%% Handy Constants
% plot colors
colorCO2 = [120,94,240]./255;
colorCH4 = [254,97,0]./255;
colorsolar = [255,176,0]./255;
colorOther = [0,176,188]./255;
Emcolor = [255,50,70]./255;

% basic constants
Constant.mia = 1.773e20; %5e18./0.029, original EarthN formula;
% number of moles (molecules of gas) corresponding to 1 atm pressure, from EONS: 
% hence atmospheric Resv values need to be corrected for molecular stoichiometry as needed
% all these values are in atm (partial pressure already)

% Thermodynamic constants
Constant.R = 8.3144; % J./mol*K; ideal gas constant = n*kB
Constant.N_avogadro = 6.022e23; % Avogadro's number (molecules./mol)

% From EONS
v.S_Pref         = 1361;             % W/m2; present day solar constant (Fs)
v.const.bol      = 5.67e-8;          % W/m^2K^4; Stefan-Boltzmann constant
v.ea.alb         = 0.29; % 0.3;      % earth albedo (~0.3 to 0.22) - Judd + 2024 value of 0.29 used here
%mr.CO2 = PGC.CO2; % ALWAYS SAME
%mr.CH4 = Flux.pCH4_max; %PGC.CH4; % 7.15e-7; % PIM constant test-case 
% Flux.pCH4_max; % Flux.pCH4_wetlands; % Flux.pCH4_marine; %
%mr.N2O = Flux.Partial_N2O.atm; % all values account for moles of molecules (not of N, C, or O atoms in molecules)
ClimateSensitivityK = 1.5; %1.9    2; % 1.5; % 3; % 3; % climate sensitivity (deg K/(W/m^2))
% constant for now to reproduce 1 deg C change in T per W/m^2 - Goldblatt, personal comm.
% can make this a T-dependent power-law (or pCO2/RF dependent - cf. He et
% al. 2023)
% could go as high as 8 per Judd et al. 2024
T_ref = 287.15; % 14 C preindustrial modern reference temperature (all GHGs at PIM reference values, modern S(t))
rf.Albedo = v.ea.alb; %PGC.Albedo; % v.ea.alb; % can be variable?
% -0.0000618557.*(4500-time./1e6) + 0.29; % 0.3; 

% spatiotemporal scaling factors!
s_yr = (60.*60.*24.*365.25); % seconds per year
SA_Earth = 5.101e14; % m^2; earth surface area, from EONS (Horne and Goldblatt 2024), essentially invariant over time (cf. Scotese et al. 2021a)
photochem_fluxscaling = s_yr.*1e4.*SA_Earth./Constant.N_avogadro; % s*cm^2*mol CH4/(yr*global area*molecules CH4); converts photochem flux units (molec/cm^2/s) to mol/yr

otherCH4Emissions = (50e12./(16.04))./1e12; % Tmol CH4/yr; geological, marine, insect, wildfire, and other sources as constant over deep time
% ignore potential GMST and pO2 dependence feedbacks for termites/marine anoxia/clathrates/wildfire frequency, etc.
% cf. Global Carbon Project, Methane Budget 2024 (citation at: https://www.globalcarbonproject.org/methanebudget/24/publications.htm)

%% Load Inputs

timeslices = flip(0:1e6:375e6).';  % time interval from 

FluxInputs = load('time_GMST_pO2_FCH4_photochemInputHiRes_updated.mat');
%FluxCH4_min = FluxInputs.FCH4_tot_hires_BC89_min.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
FluxCH4_mid = FluxInputs.FCH4_tot_hires_BC89_std.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
timeFluxB09 = FluxInputs.time_CH4FluxB09.*1e9; % rescaled to yrs from Ga BP
FluxCH4_B09 = FluxInputs.FluxCH4_B09.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
FluxCH4_max = FluxInputs.FCH4_tot_hires_BC89_max.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
FluxCH4_min = FluxInputs.FCH4_tot_hires_BC89_min.*(0.4.*1e11.*photochem_fluxscaling)./1e12; % Tmol CH4/yr
%GMST_C = FluxInputs.GMST_C_hires; % GMST 

%pCH4OutputsMax = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_strong_gammaT.mat'); % _5pt5MaPI
%SoilLossMax = pCH4OutputsMax.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8).*(photochem_fluxscaling)./1e12; % Tmol CH4/yr soil loss
%O2fluxMax = pCH4OutputsMax.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7).*(photochem_fluxscaling)./1e12; % Tmol O2/yr flux
%O2rainoutMax = pCH4OutputsMax.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,9).*(photochem_fluxscaling)./1e12; % Tmol O2/yr flux

pCH4OutputsMid = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_STDREF.mat'); % _5pt5MaPI
SoilLossMid = pCH4OutputsMid.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8).*(photochem_fluxscaling)./1e12; % Tmol CH4/yr soil loss
O2fluxMid = pCH4OutputsMid.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,7).*(photochem_fluxscaling)./1e12; % Tmol O2/yr
O2rainoutMid = pCH4OutputsMid.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,9).*(photochem_fluxscaling)./1e12; % Tmol O2/yr

%pCH4OutputsMin = load('CH4_O3_outputs_PhaneroHiRes_revisedFinal_min_gammaT.mat');
%SoilLossMin = pCH4OutputsMin.pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux(:,8).*(photochem_fluxscaling)./1e12; % Tmol CH4/yr soil loss







figure(333);clf

lc = [0 0 0];
rc = [0 0 0];
set(figure(333),'defaultAxesColorOrder',[lc; rc]);

polyx = cat(1,375,0,flip(timeslices./1e6));
polyx2 = cat(1,timeslices./1e6,flip(timeslices./1e6));
polyconst = [0, 0, 375, 375];


tiledlayout(2,1,"TileSpacing","compact","Padding","compact")
%subplot(3,1,1)

%polyFCH4em_min = cat(1,otherCH4Emissions,otherCH4Emissions,flip(FluxCH4_min));
%polyFCH4loss_min = cat(1,-SoilLossMin,-flip(FluxCH4_min)); % at steady state solution, loss = emission flux
%polyCH4other_min = cat(1,otherCH4Emissions,0,0,otherCH4Emissions).'; % constant, so can plot as rectangle
%polyFCH4soilLoss_min = cat(1,0,0,-flip(SoilLossMin)); %

%nexttile

%yyaxis left

%A = plot(polyshape(polyx,polyFCH4em_min),'FaceColor','g');
%hold on
%B = plot(polyshape(polyx2,polyFCH4loss_min),'FaceColor','m');
%hold on

%C = plot(polyshape(polyconst,polyCH4other_min),'FaceColor','b');
%hold on
%D = plot(polyshape(polyx,polyFCH4soilLoss_min),'FaceColor','y');
%hold on

%yline(0,'k--');
%hold on

%plot(timeslices./1e6,FluxCH4_min,'g-')
%hold on
%plot(timeslices./1e6,-FluxCH4_min,'m-')
%hold on
%yline(otherCH4Emissions,'b-')
%hold on
%plot(timeslices./1e6,-SoilLossMin,'y-')
%hold on

%annotation('textbox',[.34 .85-0.011+0.0045 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on

%set(gca,'XDir','reverse');
%set(gca,'xlim',[0,375]) % ,'ylim',[0.2,8]
% ,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(8,"points") % 14
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
%xlabel('Age Before Present (Ma)'); 
%xticklabels([]);
%yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6]);
%yticklabels({'0.2','','','0.5','','','','','1','2','','','5',''});
%Zeta = get(gca,'ytick');
%set(gca,'YTickLabel',Zeta);
%zoomH = zoom(gcf); 
%set(zoomH,'ActionPostCallback',{@zoom_mypostcallback});
%ylabel('CH_4 Source and Sink Fluxes (Tmol CH_4/yr)')
%ylim([3e-8,2e-5])
%pbaspect([5 2 2])

%L = legend([A],'Modeled pN_2O','FontSize',8,'Position',[0.5 0.7 0.1 0.1]); % 14,'Modern pN_2O (337 ppb)'   ,'PIM pN_2O (~0.270 ppm)'
%L.AutoUpdate = 'off';

%title('Low \gamma_T Emission Scenario')
%set(gca,'yscale','log')
%fontsize(12,"points") % 14

% L = legend([A, C, B, D],'Total CH_4 Emissions','Non-Wetland CH_4 Sources','Total Atmospheric CH_4 Sinks','Soil CH_4 Sink','FontSize',8); % ,'Position',[0.5 0.7 0.1 0.1]
% 14,'Modern pN_2O (337 ppb)'   ,'PIM pN_2O (~0.270 ppm)'
% L.AutoUpdate = 'off';

%yyaxis right

%geotimescale_Mills_JFHmod_375Ma;
%hold on
%PhanTransitions;
%set(gca,'YTickLabel',[]);

%set(gca,'XDir','reverse');
%set(gca,'xlim',[0,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(18,"points") % 14
%fontsize(12,"points") % 14

%figure(101);clf

%lc = [0 0 0];
%rc = [0 0 0];
%set(figure(101),'defaultAxesColorOrder',[lc; rc]);
%tiledlayout(1,1,"TileSpacing","compact","Padding","compact")

polyFCH4em_minmax = cat(1,16.04.*FluxCH4_max./1e3,flip(16.04.*FluxCH4_min./1e3));

nexttile

yyaxis left

size = 45;

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
semilogy(timeFluxB09./1e6,16.04.*FluxCH4_B09./1e3,'k','Marker','o','LineWidth',1.25,'LineStyle','-')
hold on
%[Nx,Ny] = boundary(polyshape(polyx2,polyFCH4em_minmax));
%patch(Nx,Ny,Emcolor) % ,'FaceColor',Emcolor,'EdgeColor',Emcolor
plot(polyshape(polyx2,polyFCH4em_minmax),'FaceColor',Emcolor,'EdgeColor',Emcolor);
hold on
semilogy(timeslices./1e6,16.04.*FluxCH4_mid./1e3,'m','Marker','none','LineWidth',2,'Linestyle','-') %,'MarkerSize',3 'none' for no marker
hold on
%
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_min./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':')
%hold on
%plot(timeslices./1e6,Flux_CH4_emissions_BC89_hires_max./1e12,'m','Marker','none','LineWidth',1.5,'Linestyle',':') % 'none' for no marker
%hold on
%
semilogy(-5,16.04.*1.31359e13./1e15,'o','MarkerFaceColor','k')
hold on
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
L = legend('CH_4 Emissions per Beerling+(2009)','',...
    'Revised CH_4 Emissions',...
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
pbaspect([2 1 1])

% %%%%%%%%%%%%%%%%%%

polyFCH4em_mid = cat(1,16.04.*otherCH4Emissions./1e3,16.04.*otherCH4Emissions./1e3,flip(16.04.*FluxCH4_mid./1e3));
polyFCH4loss_mid = cat(1,-16.04.*SoilLossMid./1e3,-flip(16.04.*FluxCH4_mid./1e3)); % at steady state solution, loss = emission flux
polyCH4other_mid = cat(1,16.04.*otherCH4Emissions./1e3,0,0,16.04.*otherCH4Emissions./1e3).'; % constant, so can plot as rectangle
polyFCH4soilLoss_mid = cat(1,0,0,-flip(16.04.*SoilLossMid./1e3)); %

%nexttile

yyaxis left

A = plot(polyshape(polyx,polyFCH4em_mid),'FaceColor','g');
hold on
B = plot(polyshape(polyx2,polyFCH4loss_mid),'FaceColor','m');
hold on

C = plot(polyshape(polyconst,polyCH4other_mid),'FaceColor','b');
hold on
D = plot(polyshape(polyx,polyFCH4soilLoss_mid),'FaceColor','y');
hold on

yline(0,'k--');
hold on

plot(timeslices./1e6,16.04.*FluxCH4_mid./1e3,'g-')
hold on
plot(timeslices./1e6,-16.04.*FluxCH4_mid./1e3,'m-')
hold on
yline(16.04.*otherCH4Emissions./1e3,'b-')
hold on
plot(timeslices./1e6,-16.04.*SoilLossMid./1e3,'y-')
hold on

%annotation('textbox',[.34 .55-0.022-0.014+0.24 .1 .1],'String','A','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
%hold on

set(gca,'XDir','reverse');
set(gca,'xlim',[0,375]) % ,'ylim',[0.2,8]
% ,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(8,"points") % 14
%set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
%set(gca,'yaxislocation','left')
xlabel('Age Before Present (Ma)'); 
%xticklabels([]);
%yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6]);
%yticklabels({'0.2','','','0.5','','','','','1','2','','','5',''});
%Zeta = get(gca,'ytick');
%set(gca,'YTickLabel',Zeta);
%zoomH = zoom(gcf); 
%set(zoomH,'ActionPostCallback',{@zoom_mypostcallback});
ylabel('CH_4 Fluxes (Pg CH_4/yr)')
annotation('textbox',[.25 .55-0.212 .1 .1],'String','B','EdgeColor','k','FitBoxToText','on','HorizontalAlignment','center')
hold on
ylim([-11,11])
%ylim([-650,650])
%pbaspect([5 2 2])

L = legend([A, C, B, D],'Wetland CH_4 Emissions','Other CH_4 Emissions','Atmospheric CH_4 Sinks','Soil CH_4 Sink','FontSize',16); % ,'Position',[0.5 0.7 0.1 0.1]
% 14,'Modern pN_2O (337 ppb)'   ,'PIM pN_2O (~0.270 ppm)'
L.AutoUpdate = 'off';

%title('Phanerozoic CH_4 Budget')
%set(gca,'yscale','log')
fontsize(16,"points") % 14

yyaxis right

geotimescale_Mills_JFHmod_375Ma;
hold on
PhanTransitions;
set(gca,'YTickLabel',[]);

set(gca,'XDir','reverse');
set(gca,'xlim',[0,375]) % 550,'xtick',time_ticks,'ylim',[0,1.5e5]
%fontsize(20,"points") % 14
fontsize(16,"points") % 14
pbaspect([2 1 1])
