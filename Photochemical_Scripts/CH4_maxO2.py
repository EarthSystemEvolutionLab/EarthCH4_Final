#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:02:23 2026

@author: johnherring
"""

import numpy as np
import scipy.io as sio
#from joblib import Parallel, delayed
from matplotlib import pyplot as plt
import utils
from IPython.display import clear_output


# modified from core code given by Nicholas Wogan (2025, personal comm.) and with additional code from Daniel Garduno-Ruiz (2025, personal comm.)

#InputCode = 2 # 1 = Anchors, 2 = HiRes, 3 = Pleistocene
# IGNORE ANCHORS

#if InputCode == 1:
 #   InputMAT = sio.loadmat('time_GMST_pO2_FCH4_photochemInputAnchors.mat') # Can use anchor points only (20-ish) or hi-res interpolation
 #   Range_timeUV = InputMAT['timeUVb'] # already in Ga BP   # [0.0, 0.6, 2.4] # 3 cases, Ga BP
  #  Range_GMST = InputMAT['AnchorGMST_C'] # in deg C GMST[45] #, 14] #, 20, 25, 30, 35, 40, 45] # 8 cases, GMST in C
  #  Range_pO2 = InputMAT['AnchorpO2'] # in PAL O2.  [0.001, 0.01, 0.1, 0.5, 1, 1.5, 2] # 7 cases, pO2 in multiples of PAL
    #Range_FCH4 = InputMAT['FCH4_tot_lo'] # LOW T-sensitivity reconstruction in xPIM emission flux  [40, 70, 100] #[0.584, 1, 5.84, 11.67] # 4 values spanning Beerling 2009 range, 
   # Range_FCH4 = InputMAT['FCH4_tot_hi'] # HIGH T-sensitivity reconstruction in xPIM emission flux
#elif InputCode == 2:

InputMAT = sio.loadmat('time_GMST_pO2_FCH4_photochemInputHiRes_updated.mat')
Range_timeUV = InputMAT['timeUVb_hires']#[113] # already in Ga BP   # [0.0, 0.6, 2.4] # 3 cases, Ga BP
Range_GMST = InputMAT['GMST_C_hires']#[113] # in deg C GMST[45] #, 14] #, 20, 25, 30, 35, 40, 45] # 8 cases, GMST in C
Range_pO2 = InputMAT['pO2_max']#[113] # in PAL O2.  [0.001, 0.01, 0.1, 0.5, 1, 1.5, 2] # 7 cases, pO2 in multiples of PAL
#Range_FCH4 = InputMAT['FCH4_tot_hires_lo'] # LOW T-sensitivity reconstruction in xPIM emission flux  [40, 70, 100] #[0.584, 1, 5.84, 11.67] # 4 values spanning Beerling 2009 range, 
#Range_FCH4 = InputMAT['FCH4_tot_hires_Macro_max'] # HIGH T-sensitivity reconstruction in xPIM emission flux
Range_FCH4 = InputMAT['FCH4_tot_hires_BC89_std']#[113] # HIGH T-sensitivity reconstruction in xPIM emission flux
#Range_FCH4 = InputMAT['FCH4_tot_hires_Macro_min'] # HIGH T-sensitivity reconstruction in xPIM emission flux
#else:

  #  InputMAT = sio.loadmat('time_GMST_pO2_FCH4_photochemInputPleisto.mat')
  #  Range_timeUV = InputMAT['timeUVb_Pleisto5ka'] # already in Ga BP   # [0.0, 0.6, 2.4] # 3 cases, Ga BP
  #  Range_GMST = InputMAT['GMST_C_Pleisto5ka'] # in deg C GMST[45] #, 14] #, 20, 25, 30, 35, 40, 45] # 8 cases, GMST in C
   # Range_pO2 = InputMAT['pO2_Pleisto5ka'] # in PAL O2.  [0.001, 0.01, 0.1, 0.5, 1, 1.5, 2] # 7 cases, pO2 in multiples of PAL
    #Range_FCH4 = InputMAT['FCH4_tot_Pleisto_lo5ka'] # LOW T-sensitivity reconstruction in xPIM emission flux  [40, 70, 100] #[0.584, 1, 5.84, 11.67] # 4 values spanning Beerling 2009 range, 
 #   Range_FCH4 = InputMAT['FCH4_tot_Pleisto_hi5ka'] # HIGH T-sensitivity reconstruction in xPIM emission flux
# OR
# InputMAT = sio.loadmat('time_GMST_pO2_FCH4_photochemInputPleisto.mat')

#Range_timeUV = InputMAT['AnchorOutputArr'][:,1] # already in Ga BP   # [0.0, 0.6, 2.4] # 3 cases, Ga BP
#Range_GMST = InputMAT['AnchorOutputArr'][:,2] # in deg C GMST[45] #, 14] #, 20, 25, 30, 35, 40, 45] # 8 cases, GMST in C
#Range_pO2 = InputMAT['AnchorOutputArr'][:,3] # in PAL O2.  [0.001, 0.01, 0.1, 0.5, 1, 1.5, 2] # 7 cases, pO2 in multiples of PAL
#Range_FCH4 = InputMAT['AnchorOutputArr'][:,4] # LOW T-sensitivity reconstruction in xPIM emission flux  [40, 70, 100] #[0.584, 1, 5.84, 11.67] # 4 values spanning Beerling 2009 range, 
#Range_FCH4 = InputMAT[:,5] # HIGH T-sensitivity reconstruction in xPIM emission flux


# as multiples of the modern (post-industrial) imposed reference photochem flux (0.4e11 molec/cm^2/s, scales to ~171 Tg CH4/yr)
# here translate to 100, ~171 (PIM), 1000, and 2000 Tg CH4/yr
# Total number of test cases = 7*3*4 = 84 per script
# NOMINAL = [2 5 1 2]

CH4_Array = np.zeros((len(Range_GMST), 9)) # 8 outputs in product array
#1 7 3 4  initializes empty "zero" 4-D array in which to store interpolation matrix
# 4 across to contain pCH4, Tau_CH4, pO3 (DU) and surface tropo pO3 + bug check (last 2 columns must be identical)
# should be accessed as (ii, jj, kk, ll) for 8, 7, 6, 4
counter = 0 # initializes progress counter

for ii in range(len(Range_GMST)):
   # for jj in range(len(Range_pO2)):
      #  for kk in range(len(Range_timeUV)):
         #   for ll in range(len(Range_FCH4)):
                # runs script iteratively
                

                TauSun = Range_timeUV[ii] # sets stellar age

                pc = utils.EvoAtmosphereJohn(
                    'input/zahnle_earth.yaml', # Chemical reactions
                    'input/settings_Earth_old.yaml', # Bunch of settings are in here
                    'input/atmosphere_ModernEarth.txt', # The "initial condition"
                    age_of_sun=TauSun[0] # 0.0 # set age of sun here in billions of years ago.
                )
                pc.var.verbose = 0 # Turn off printing

                # Massie and Hunten 1981 N2O profile data (Wogan 2025 personal comm.)
                # https://github.com/Nicholaswogan/planetary_atmosphere_observations/blob/bc42c45f86735bffee8038c3a49770a61ed29c6c/Earth.yaml
                #data_alt_N2O = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
                #data_mix_N2O = [3.16e-7, 3.16e-7, 3.16e-7, 3.04e-7, 2.51e-7, 1.77e-7, 1.12e-7, 5.1e-8, 1.6e-8, 5.4e-9, 2.1e-9, 1.1e-9, 6.2e-10, 3.5e-10]

                # Massie and Hunten 1981 O3 profile data (Wogan 2025 personal comm.)
                #data_alt_O3 = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
                #data_mix_O3 = [3.0e-8, 3.2e-8, 5.0e-8, 3.3e-7, 2.2e-6, 5.7e-6, 6.5e-6, 6.2e-6, 4.7e-6, 2.9e-6, 1.6e-6, 8.6e-7, 5.6e-7, 4.1e-7]

                # Set the oxygen partial pressure here (dynes/cm^2)
                scalingO2 = 1*Range_pO2[ii] # lo = 0.001 # hi = 2 # nominal = 1
                #print('ScalingO2 =', scalingO2[0])
                pc.set_lower_bc('O2', bc_type='press', press=scalingO2[0]*0.212730e6) # Modern value = 0.212730e6 dynes/cm^2 
                # 1 dyne/cm^2 = 0.1 Pa or 1/1.013e6 atm

                # Set N2O flux into atmosphere
                # Modern Earth flux = 1.0e9 molecules/cm^2/s, increased to 1.5e9 per Schwieterman et al. 2022 (cf. Wogan, personal comm.)
                # vdep = 1e-4 cm/s, a deposition velocity, which seems to work for Modern Earth # halved to give more reasonable surface pN2O given preindustrial fluxes
                # height = -1, means that the flux is emitted only into surface layer'ArithmeticError
                scalingF_N2O = 1 # 1 #0.812 #*0.946 # 1.367 # 1 
                # 1.143 = 0.4 Tmol N2O/yr (preindustrial Schwieterman et al. 2022)
                # 1.367 = case for modern post-industrial atmo sink flux from Tian et al. 2024 = 0.48 Tmol N2O/yr, would be same as emission flux if modern pN2O (337 ppb) were at steady-state
                # 1.367*0.946 = modern 2020 atmo sink of N2O (pseudo-steady-state emission flux equivalent) scaled down to 2000 levels per Table 3, Tian et al. 2024 (stratospheric sink 12.2 vs. 12.9 modern, note not the same as modern 13.4 Tg N/yr flux in Fig. 1)
                # (1.5e9 molec/cm^2/s), renders 310 ppb close to 337 (modern) or 321 (Prather et al. 2015) ppb!
                # 0.875 # 0.812 
                # # lo = 0.25 # hi = 250 # nominal = 1
                # 0.875 * 1.5e9 = 3.5e11 mol N2O/yr global, current EarthN2O output => 275 ppb with 5e-5 vdep! 
                # 0.812 = 3.25e11 mol N2O/yr global, 
                FN2O = scalingF_N2O*1.3125e9 # 0.35 Tmol N2O/yr global, current model output and in Prather et al. 2015 range (near upper limit)
                # 1.5e9 Schwieterman, almost exactly correct rescaling of 0.4 Tmol N2O/yr global emission to atmo
                ScalingVdep = 0.3 # 0.7 # 0.7 #1 # 0.5 #0.5 # 1
                # 0.3 works best with 14 C GMST and moist adiabatic profile (isothermal stratosphere), decent pO3 and pN2O profile and good Tau_N2O (slightly high tropo pO3, but ok)
                # 0.7 works great! 7e-5 cm/s # this slightly deteriorates stratospheric pO3 fit to data, but not very significantly
                # overall, improved pN2O fit and slightly poorer stratospheric pO3 fit seem to balance, reasonable match to data upheld whether using 7e-5 or 1e-4 vdep
                # 0.6 not much better, but lengthens lifetime by 6 years and exaggerates yN2O up to 296 ppb instead of 283 ppb
                # hence 0.6 OK, but 0.7 probably better (may be conservative???)
                vdep_N2O = ScalingVdep*1e-4 #1e-4
                #pc.set_lower_bc('N2O', bc_type='vdep + dist flux', flux=FN2O, vdep=vdep_N2O, height=-1)

                scalingF_CH4 = Range_FCH4[ii] # 0.4-0.5 preindustrial wetland flux of CH4 (Beerling et al. 2009) only 10.66 (midrange) to simulated 13.13 Tmol CH4/yr, vs. 26.7 Tmol/yr in std forcing below
                FCH4 = 0.4*scalingF_CH4[0]*1e11 # standard 1e11 flux from Wogan (personal comm) = 26.7 Tmol CH4/yr global = 428.4 Tg CH4/yr, here scaled down to match 171 Tg/yr CH4 (midrange estimate from prior studies in Table 2 Beerling + 2009)
                # 267.1291 = conversion from molec/cm^2/s to moles/yr global!
                vdep_CH4 = 0.0005 # = 0.002/4; based on simple diffusion-limitation parameterization inspired by Ridgwell+1999 #0 # per Wogan (personal comm) and 
                pc.set_lower_bc('CH4', bc_type='vdep + dist flux', flux=FCH4, vdep=vdep_CH4, height=-1)


                T = Range_GMST[ii] # 14 
                # 20 # 20 # 24 Celsius, reasonable preindustrial modern average (14 rather low)
                surf_temp = 273.15 + T[0] # Celsius to Kelvin
                # FROM Garduno-Ruiz 2025, personal comm.
                # moist adiabat temperature profiles (units: K)
                temp_profiles = np.loadtxt('input/moist_adiabat_temp_profiles.txt')
                # eddy diffusivity profiles (units: cm^2/s)
                eddy_diff_profiles = np.loadtxt('input/moist_adiabat_eddy_profiles.txt')
                # pressure profiles (units: Pa, needs to be converted to bar)
                press_profiles = np.loadtxt('input/moist_adiabat_press_profiles.txt')*10 #/1e5 renders dynes/cm^2 from Pa

                def interp_profiles(profiles, surf_temp):
                    '''Interpolate temperature or eddy diffusivity profile given the surface 
                    temperature   
                    input:
                        profiles (np.array 1d): array of known profiles
                        surf_temp (float): surface temperature of interpolated profile (k)
                    returns:
                        interpolated_profile (np.array)
                    '''
                    temps = np.arange(240, 360.5, 0.5)
                    nz = profiles.shape[1]
                    interpolated_profile = np.zeros(nz)
                    for i in range(0, nz):
                        level = profiles[:, i]
                        interpolated_profile[i] = np.interp(surf_temp, temps, level)
                    return interpolated_profile



                def set_atm_structure(pc, surf_temp):
                    '''
                    Sets temperature, pressure and eddy diffusivity profiles, and tropopause
                    height in photochemical model 
                    input:
                        pc (photochempy object): instance of photochempy class
                        surf_temp (float): surface temperature (k)
                    returns:
                        none
                    '''
                    
                    temp_profile = interp_profiles(temp_profiles, surf_temp)
                    edd_profile = interp_profiles(eddy_diff_profiles, surf_temp)
                    press_profile = interp_profiles(press_profiles, surf_temp)
                    
                    # index of tropopause
                    jtrop = np.argmin(temp_profile) + 1

                    # set temperature and eddy diffusivity profiles and surface pressure
                    #pc.set_temperature(temp_profile, trop_alt = pc.var.z[jtrop])
                    
                    

                #len(pc.var.temperature)
                #pc.set_temperature(temp_profile, trop_alt = pc.var.z[jtrop])
                #set_atm_structure(pc, surf_temp)
                temp_profile_raw = interp_profiles(temp_profiles, surf_temp)
                edd_profile_raw = interp_profiles(eddy_diff_profiles, surf_temp)
                press_profile_raw = interp_profiles(press_profiles, surf_temp)

                temp_profile = temp_profile_raw[::2]
                edd_profile = edd_profile_raw[::2]
                press_profile = press_profile_raw[::2] #units of dynes/cm^2 per requirements of function (see documentation) 

                jtrop = np.argmin(temp_profile) + 1 # locates tropopause "cold trap" elevation

                tropoP = press_profile[jtrop] # converts bar P to dynes/cm^2

                #print(jtrop)

                pc.set_press_temp_edd(press_profile,temp_profile,edd_profile,trop_p=tropoP)
                #print(tropoP)

                
                #
                # Initialize an integrator
                pc.initialize_robust_stepper(pc.wrk.usol) 
                pc.find_steady_state()
                #

                # Integrate to equilibrium and plot along the way
                #while True:
                    
                #    clear_output(wait=True)

                    # Make plot
                    #fig,ax = plt.subplots(1,1,figsize=[6,5])
                    
                    # This function returns the state of the atmosphere in dictionary
                    # sol = pc.mole_fraction_dict()

                    # Plots species
                    #species=['H2O','O2','O3','NO','N2O']
                    #for i,sp in enumerate(species):
                    #    ax.plot(sol[sp],pc.var.z/1e5, label=sp) # variation to plot elevation not pressure profile
                       

                    #ax.plot(data_mix_N2O,data_alt_N2O,label="N2O Profile Data") #
                    #ax.plot(data_mix_O3,data_alt_O3,label="O3 Profile Data") #
                    #ax.plot(pc.var.temperature,pc.var.z/1e5,label="T Profile") #
                    #ax.set_xscale('log')
                    #ax.set_yscale('log')
                    #ax.invert_yaxis()
                    #ax.grid(alpha=0.4)
                    #ax.set_xlim(1e-10,1)
                    # ax.set_ylabel('Pressure (bars)') 
                    # pressure version of plot only goes to 1 bar (~1 atm) = constant P preserved
                    #ax.set_ylabel('Elevation above surface (km)')
                    #ax.set_xlabel('Mixing ratio')
                    #ax.legend(ncol=1,bbox_to_anchor=(1,1.0),loc='upper left')
                    #ax.text(0.02, 1.04, 't = '+'%e s'%pc.wrk.tn, \
                     #   size = 15,ha='left', va='bottom',transform=ax.transAxes)
                    #plt.show()

                    # break
                 #   for i in range(50):
                 #       give_up, converged = pc.robust_step()
                 #       if give_up or converged:
                 #           break
                 #   if give_up or converged:
                 #       break    

                # After converged check N2O mole fraction at surface
                sol = pc.mole_fraction_dict()
                print('Surface yCH4 = ',sol['CH4'][0])
                #print(sol['O2'][0])

                CH4_Array[ii,0] = sol['CH4'][0] # ground pCH4: mixing ratio (assumedly partial pressure at ~ 1 atm surface P) of CH4 at surface, so ground pCH4 
                # (should be near-constant through tropo)

                # This is the flux of N2O needed from the surface to the atmosphere
                # to maintain the current N2O concentration
                surface_flux_CH4 = pc.gas_fluxes()[0]['CH4']
                print('Surface flux CH4 = %.2e'%(surface_flux_CH4))
                CH4_Array[ii,4] = surface_flux_CH4
                # Flux I set into the atmosphere at the lower boundary
                boundary_flux = - sol['CH4'][0]*pc.wrk.density[0]*vdep_CH4 + FCH4
                print('Boundary condition flux of CH4 = %.2e'%boundary_flux)
                CH4_Array[ii,5] = boundary_flux # columns 5 and 6 of output (surface flux and boundary flux) must be equal! Bug/failure check 

                # The fluxes balance because we are in a steady state.

                # atmospheric lifetime script 
                # pl = pc.production_and_loss('N2O',pc.wrk.usol)
                # net = -np.sum(pl.integrated_production) + np.sum(pl.integrated_loss)
                # net = np.sum(pl.integrated_production) - np.sum(pl.integrated_loss) # correct net flux, but signs should be reversed for clarity
                ind = pc.dat.species_names.index('CH4')
                dz = pc.var.z[1] - pc.var.z[0]
                col = np.sum(pc.wrk.usol[ind,:]*dz)
                # col = np.sum(pc.wrk.usol[ind,:]*dz)

                #print('N2O lifetime in seconds =', col/net)

                #print('N2O molecular lifetime in years =', (col/net)/31536000)

                Tau_CH4 = (col/FCH4)/(60*60*24*365.25) # N2O "emission lifetime" (column burden/emission flux) converted to years

                print(ii+1,'CH4 emission lifetime in years =', Tau_CH4) # (col/FN2O)/31536000)
                CH4_Array[ii,1] = Tau_CH4 # accompanies surface pCH4

                
                counter = counter + 1 # updates progress bar
                print('Progress %: ',100*counter/len(Range_GMST)) # 84.  168) #1680) # prints progress bar
                print(pc.check_for_convergence())
                #print(pc.var.temperature[0])
                #print(pc.var.edd[0])
                #print(sol['pressure'][0]/1e6)
                N = pc.gas_fluxes()
                
                #Tau_CH4_Array[ii, jj, kk, ll] = Tau_CH4 # loads most recent output into matching array cell

                    # finds O3 column density and saves it in Dobson Units
                indO3 = pc.dat.species_names.index('O3')
                colO3 = (np.sum(pc.wrk.usol[indO3,:]*dz))/2.69e16 # crucial output for O3 column thickness, 300 DU typical broadly
                # with conversion to Dobson Units (DU) from molecules/cm^2 column (using 2.69e16 molec/cm^2 = 1 DU) - see ozonewatch.gsfc.nasa.gov/facts/dobson_SH.html
                CH4_Array[ii,2] = colO3 # 3rd column contains pO3 column in DU
                CH4_Array[ii,3] = sol['O3'][0] # resolves surface tropo pO3 in 4th column

                # index of strato O3 layer peak
                #jstratoO3 = np.argmax(sol['O3']*press_profile) + 1 # multiplies atmospheric pressure profile by yO3 profile (yO3*P = pO3 per Dalton's Law) and finds maximum pO3
                #CH4_Array[ii,6] = sol['OH'][0] # OH in col 7
                #CH4_Array[ii,7] = sol['CO'][0] # CO in col 7
                CH4_Array[ii,6] = pc.gas_fluxes()[0]['O2'] # O2 fluxes to maintain prescribed pO2
                CH4_Array[ii,7] = sol['CH4'][0]*pc.wrk.density[0]*vdep_CH4 # soil uptake/dry deposition CH4 flux
                #surface_flux_CH4/pc.gas_fluxes()[0]['O2'] # ratio of CH4:O2 fluxes
                pl=pc.production_and_loss('O2',pc.wrk.usol)
                O2rainout = (np.sum(pl.loss[:, pl.loss_rx.index('rainout')]*dz))
                #print('O2 rainout = %.2e'%(O2rainout))
                CH4_Array[ii,8] = O2rainout #

                 # collates output for each timestep
                #Tau_N2O[isp] = (colN2O/FN2O)/31536000
                #Col_O3[isp] = colO3
                #pO3_stratomax[isp] = sol['O3'][jstratoO3]*press_profile[jstratoO3] # to compare maximum pO3 in stratosphere against modern O3 layer peak

# AFTER hyperloop terminates, save final matrix output to .mat file for use in model

Tau_CH4_arr = {'pCH4_TauCH4_colO3_pO3tropo_surf_bound_O2flux_soilDryDepFlux': CH4_Array}

#if InputCode == 1:
  #  sio.savemat('CH4_O3_outputs_PhaneroAnchors.mat',Tau_CH4_arr)
#elif InputCode == 2:

sio.savemat('CH4_O3_outputs_PhaneroHiRes_revisedFinal_maxO2.mat',Tau_CH4_arr) #revisedFinal_strong_gammaT
