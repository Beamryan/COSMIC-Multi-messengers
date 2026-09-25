# -*- coding: utf-8 -*-
# Copyright (C) Katelyn Breivik (2017 - 2019)
# https://zenodo.org/record/3905313#.YJ6QjJNKjyg  -- cite if use this
#
# MW Maker is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MW Maker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# --------------------------------------------------------------------
# REQUIRES:
#    -- GW_calcs.py   (same directory)  DEPRECATED
#    -- Legwork       http://legwork.readthedocs.io
# --------------------------------------------------------------------
# VERSION AUTHOR :  Shane L. Larson (s.larson@northwestern.edu)
#
#  11 Feb 2021 Modifications
#  -- Requirements in Comments
#  -- Modified for my directory setup on EURISKO
#
# 16 Mar 2022 Modifications
#     Made MC_samp as file to import, rather than from COSMIC
#     also updated physical constants to be consistent with MC_samp
#     CORRECTED ERROR IN PARSEC CONVERSION --> x10^16 meters (error was x10^19)
#     though it doesn't look like it is used anywhere...
#
# -----------------------------------------------------------------------
# VERSION AUTHOR : Natalie Gottschlich (natalie.e.gottschlich@gmail.com)
# 
# 10/31/2022 Modifications:
#  -- Changed variable 'model' to 'gx_component' to avoid degeneracy with variable in MC_samp
#
#  5/1/2023 Modifications:
#  -- Re-calculated radii before roche-lobe selection to ensure the correct radii are used
#     (fixing COSMIC issue)
#
# -----------------------------------------------------------------------
# VERSION AUTHOR : Shane L. Larson (slarson@clarkson.edu)
#
# 04/23/2025 Modifications
# -- Updated calles to GW_calcs to us Legwork instead; evolutions for NS and BH
#    were not integrating correctly with the original integrator in GW_calcs
#
# -----------------------------------------------------------------------


'''`MW_Maker`
'''


# =============================
#  Python imports
# =============================

# from cosmic import MC_samp
# Extracted MC_samp from COSMIC 3.3 and included here as standalone file
import MC_samp
import numpy as np
import pandas as pd
import legwork       # replaces legacy:  import GW_calcs
import astropy.units as u

# ==========================================================================
#  Setup CONSTANTS
# ==========================================================================
# sll: these should all be  updated to higher precision; especially the
#      yr_sec, which in principle affects LISA data analysis; don't know
#      how it propagates through the COSMIC production though
#      Follow values I use in my C headers, documented provenance for values

Rsun_au = 1/215.032 # Radius of the sun in AU -- AU/solar radii
day_yr = 1/365.25   # interesting -- what is this and is it consistent with yr_sec? don't see it used
yr_sec = 3.155e7    # seconds in a year  [Sidereal Year = 31558149.7635456 sec]
G = 6.67384e-11     # SI units
c = 299792458.0     # SI units: m/s
Msun = 1.9891e30    # SI units: kg
parsec = 3.0856775814913673e16  # SI units: meters

# ==========================================================================
# FUNCTION: wdRadius
# AUTHOR  : Natalie Gottschlich
# ==========================================================================

# NG: Added this function 5/1 to calculate the radii
def wdRadius(m_wd):
    """Gives radius of WD from mass of WD following Hurley 2000, eq. 91
    
    Parameters
    ----------
    m_wd : series/array
        mass of the white dwarf in solar masses
        
    
    Returns
    -------
    r_wd : float array
        radius of the white dwarf in solar radii
    """
    
    r_ns = 1.4e-5
    m_ch = 1.44
    
    r_wd = 0.0115 * np.sqrt((m_ch / m_wd)**(2/3) - (m_wd / m_ch)**(2/3))
    
    radii = np.where(r_wd < r_ns, r_ns, r_wd)
    
    return radii


# ==========================================================================
# FUNCTION: porb_from_sep
# ==========================================================================

def porb_from_sep(m1, m2, sep):
    """Gives orbital period from m1, m2, sep
    following Kepler III
    
    Parameters
    ----------
    m1 : float/array
        primary mass in msun
        
    m2 : float/array
        secondary mass in msun
        
    sep : float/array
        binary separation in au
        

    Returns
    -------
    porb : float/array
        orbital period in years
    """
    
    porb_2 = sep**3/(m1+m2)  # [sll] sep in AU, m in Msun -- Porb in YR 
    return porb_2**0.5


# ==========================================================================
# FUNCTION: peak_gw_freq
# ==========================================================================

# (sll) -- 23 April 2025
# Copied this function from the old GW_calcs to here, since this function
# is not obviously in Legwork; right now just trying to get the code to
# work -- will evaluate how this is used later to do LISA band cutoffs
# (that's what it looks like) and decide the right way to include/implement
# it.
# KMB says this is from Linqing Wen's paper:

def peak_gw_freq(m1, m2, ecc, porb):
    """Computes the peak gravitational-wave frequency for an
    eccentric binary system. Units are SI

    Parameters
    ----------
    m1 : float or array
        primary mass [kg]
    m2 : float or array    
        secondary mass [kg]
    ecc : float or array
        eccentricity
    porb : float or array
        orbital period [s]

    Returns
    -------
    f_gw_peak : float or array
        peak gravitational-wave frequency [Hz]
    """

    # convert the orbital period into a separation using Kepler III
    sep_m = (G/(4*np.pi**2)*porb**2*(m1+m2))**(1./3.)

    # peak of the peters harmonic sequence; from Linqing Wen's paper... 2003? 2001?
    f_gw_peak = ((G*(m1+m2))**0.5/np.pi) * (1+ecc)**1.1954/(sep_m*(1-ecc)**2)**1.5
    return f_gw_peak



# ==========================================================================
# FUNCTION: SFH
# ==========================================================================

def SFH(n_pop, gx_component):
    if gx_component == 'ThinDisk':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 10000, n_pop)
    elif gx_component == 'Bulge':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 1000, n_pop)
    elif gx_component == 'ThickDisk':
        # times in Myr following COSMIC
        t_birth = np.random.uniform(0, 1000, n_pop)
    return t_birth



# ==========================================================================
# FUNCTION: t_evol_GW    -- returns evolution under GW time in Myr
#                        -- Assumes birth_time and formation_time
#                           are in Myr
# ==========================================================================

def t_evol_GW(birth_time, formation_time, gx_component):

    # import pdb
    # pdb.set_trace()

    if gx_component == 'ThinDisk':
        # times in Myr following COSMIC
        t_evol_GW = 10000 - (birth_time + formation_time)
    elif gx_component == 'Bulge':
        # times in Myr following COSMIC
        t_evol_GW = 10000 - (birth_time + formation_time)
    elif gx_component == 'ThickDisk':
        # times in Myr following COSMIC
        t_evol_GW = 11000 - (birth_time + formation_time)

    # import pdb
    # pdb.set_trace()

    return t_evol_GW



# ==========================================================================
# FUNCTION: R_RL
# ==========================================================================

def R_RL(q, a):
    """ Computes the Eggleton 1983 Roche Lobe radius
    in units of supplied semimajor axis
    
    Parameters
    ----------
    q : float or array
        m_RL/m_companion
    
    a : float or array
        semimajor axis
        
    Returns
    -------
    RL : float or array
        Roche radius
    """
    numerator = 0.49*q**(2./3.)
    denominator = (0.6*q**(2./3.)+ np.log(1+q**(1./3.)))
    RL = a*numerator/denominator
    return RL


# ==========================================================================
# FUNCTION: read_dat
# ==========================================================================
# ### (sll): modded this function to take a filename, so it doesn't have to
#            construct the filename, per my block comment below.
#            Original function header call was:
# def read_dat(dat_path, kstar_types):
def read_dat(dat_path, dat_file, kstar_types):
    if len(kstar_types) == 2:
        kstar_1 = kstar_types[0]
        kstar_2 = kstar_types[1]
        kstar_plot = '{0}_{1}'.format(kstar_1, kstar_2)
        kstar_1 = [kstar_1]
        kstar_2 = [kstar_2]
    else:
        kstar_1 = kstar_types[0]
        kstar_2_lo = kstar_types[1]
        kstar_2_hi = kstar_types[2]
        kstar_plot = '{0}_{1}_{2}'.format(kstar_1, kstar_2_lo, kstar_2_hi)
        kstar_1 = [kstar_1]
        kstar_2 = np.arange(kstar_2_lo, kstar_2_hi)

    
    # NOTE: this is the dataframe so to get the total mass we need to take the final mass
    # ### (sll) Here I've commented out the original lines in kmb's MW_maker.py since they
    #           don't match the filenames currently generated by COSMIC. The "DeltaBurst"
    #           part of the filename seems to have been replaced by the start and duration
    #           of star formation, eg: SFstart_13700.0_SFduration_0.0
    #           I will figure out how to automate this in the future, but for the moment
    #           change these at runtime.  The full name I currently generate is:
    #               dat_kstar1_10_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.017
    # m_sim = pd.read_hdf(dat_path+'dat_DeltaBurst_short_'+kstar_plot+'.h5', key='mass_stars')
    # m_sim_tot = m_sim.max()[0]
    # conv = pd.read_hdf(dat_path+'dat_DeltaBurst_short_'+kstar_plot+'.h5', key='conv')
    # return conv, m_sim, m_sim_tot

    m_sim = pd.read_hdf(dat_path+dat_file, key='mass_stars')
    m_sim_tot = m_sim.max()[0]
    conv = pd.read_hdf(dat_path+dat_file, key='conv')
    return conv, m_sim, m_sim_tot



# ==========================================================================
# FUNCTION: GW_evol
# ==========================================================================

def GW_evol(pop):
    """Computes the final state of binaries in pop
    based on GW evolution over t_evol_GW years
    
    Parameters
    ----------
    pop : pandas.DataFrame
        contains all the binary info for a binary pop that
        has been scaled to the Galaxy
        
    Returns
    -------
    sep_final : array
        final separation after GW evolution
    ecc_final : array
        final eccentricity after GW evolution
    
    """

    sep_finalLW = np.zeros(len(pop)) * u.Rsun  # sll, 25 Apr 2025: added * u.Rsun for Legwork; was sep_final = np.zeros(len(pop))
    sep_final = np.zeros(len(pop))             # sll, 25 Apr 2025: keep for final output without units
    ecc_final = np.zeros(len(pop))

    # units: msun, au, years
    ecc = np.array(pop.ecc)
    m1 = np.array(pop.mass_1) * u.Msun       # sll, 25 Apr 2025: added * u.Msun for Legwork; was m1 = np.array(pop.mass_1)
    m2 = np.array(pop.mass_2) * u.Msun       # sll, 25 Apr 2025: added * u.Msun for Legwork; was m2 = np.array(pop.mass_2)
    sep = np.array(pop.sep) * u.Rsun         # sll, 25 Apr 2025: u.Rsun for Legwork; was: sep = np.array(pop.sep*Rsun_au)
    times = np.array(pop.t_evol_GW) * u.Myr  # sll, 25 Apr 2025: u.Myr for Legwork; was: times = np.array(pop.t_evol_GW*1e6)
    circ_ind, = np.where(ecc <= 0.01)
    ecc_ind, = np.where(ecc > 0.01)
    
# ---------------------------------------------------------------------------------
# ECCENTRINC BINARIES: using Legwork
# ---------------------------------------------------------------------------------
    # sll 30 Apr 2025: added to have Legwork calculate if eccentric binaries merge
    #     so they can be filtered out.
    
    t_merge_ecc = legwork.evol.get_t_merge_ecc(ecc_i=ecc[ecc_ind],
                                               m_1=m1[ecc_ind],
                                               m_2=m2[ecc_ind],
                                               a_i=sep[ecc_ind],
                                               exact=False)

    ecc_ind_merge, = np.where(t_merge_ecc < times[ecc_ind])  # these ecc binaries have merged before today
    ecc_ind_alive, = np.where(t_merge_ecc > times[ecc_ind])  # these ecc binaries have not merged, so are UCBs today
    sep_finalLW[ecc_ind[ecc_ind_merge]] = 0.0

# ---------------------------------------------------------------------------------
# Loop through eccentric binaries, evolve using Legwork
# ---------------------------------------------------------------------------------
    for ind in ecc_ind_alive:
# This is the old call to GW_calcs in MW_makerN2:
#        sep_final[ind], ecc_final[ind] = GW_calcs.peters_sep_ecc(m1=m1[ind], m2=m1[ind],
#                                                                 a_0=sep[ind],
#                                                                 e_0=ecc[ind],
#                                                                 t_evol=times[ind])
# Legwork tutorial example:
#     ecc_evol, f_orb_evol = evol.evol_ecc(ecc_i=ecc_i[inspiral],
#                                          m_1=m_1[inspiral], m_2=m_2[inspiral],
#                                          f_orb_i=f_orb_i[inspiral],
#                                          timesteps=timesteps, t_before=10 * u.yr)

        ecc_final[ind], sep_finalLW[ind] = legwork.evol.evol_ecc(ecc_i=ecc[ind],
                                          m_1=m1[ind], m_2=m2[ind],
                                          a_i=sep[ind],
                                          t_evol=times[ind],
                                          t_before= 10 * u.yr,
                                          n_step=1,
                                          output_vars=["ecc","a"],
                                          avoid_merger=False)
        
    
# ---------------------------------------------------------------------------------
# CIRCULAR BINARIES: using Legwork
# ---------------------------------------------------------------------------------
# This is the old call to GW_calcs in MW_makerN2:
#    t_merge_circ = GW_calcs.peters_t_merge_circ(m1=m1[circ_ind],
#                                       m2=m2[circ_ind], 
#                                       a0=sep[circ_ind])
# Legwork tutorial example:
#   t_merge = evol.get_t_merge_circ(m_1=m_1, m_2=m_2, f_orb_i=f_orb_i).to(u.yr)
    t_merge_circ = legwork.evol.get_t_merge_circ(m_1=m1[circ_ind],
                                                 m_2=m2[circ_ind], 
                                                 a_i=sep[circ_ind])


    ind_merge, = np.where(t_merge_circ < times[circ_ind]) # these have merged before today
    ind_alive, = np.where(t_merge_circ > times[circ_ind]) # these have not merged, so are DWD today
    sep_finalLW[circ_ind[ind_merge]] = 0.0
    
    # import pdb
    # pdb.set_trace()

    # print('Number of Merged Circ Systems = ',np.where(t_merge_circ < times[circ_ind]))
    
# This is the old call to GW_calcs in MW_makerN2:
#    sep_final[circ_ind[ind_alive]] = GW_calcs.peters_a_circ(a_0=sep[circ_ind[ind_alive]], 
#                                                            m1=m1[circ_ind[ind_alive]], 
#                                                            m2=m2[circ_ind[ind_alive]], 
#                                                            t_evol=times[circ_ind[ind_alive]])
# Legwork tutorial example:
#   a_evol = evol.evol_circ(m_1=m_1, m_2=m_2, f_orb_i=f_orb_i, output_vars="a")
# 28 Apr 2025: added no.squeeze() to remove extraneous 1 dimensional array entry
#    sep_finalLW[circ_ind[ind_alive]] = np.squeeze(legwork.evol.evol_circ(a_i=sep[circ_ind[ind_alive]],
#                                                            m_1=m1[circ_ind[ind_alive]],
#                                                            m_2=m2[circ_ind[ind_alive]],
#                                                            t_evol=times[circ_ind[ind_alive]],
#                                                            n_step=1,
#                                                            output_vars="a"))

    sep_out = np.squeeze(legwork.evol.evol_circ(a_i=sep[circ_ind[ind_alive]],
                                                            m_1=m1[circ_ind[ind_alive]],
                                                            m_2=m2[circ_ind[ind_alive]],
                                                            t_evol=times[circ_ind[ind_alive]],
                                                            n_step=2,
                                                            output_vars="a"))
                                                            
    sep_finalLW[circ_ind[ind_alive]] = sep_out[:,1]       # kmb helped here, 12 Jun 2025

    ecc_final[circ_ind] = 0.0

    # convert the separation back to Rsun
    # sep_final = sep_final/Rsun_au
    
    # import pdb
    # pdb.set_trace()

    sep_final = sep_finalLW.value   # this extracts the number without the unit
                                    # Legwork uses units, rest of codes do not
                                    # Legwork units for sep_final are: solRad

    return sep_final, ecc_final



# ==========================================================================
# FUNCTION: gx_sample
# ==========================================================================

def gx_sample(conv, gx_component, n_pop):
    # sample the Gx population and assign birth and GW evol times
    pop = conv.sample(n_pop, replace=True)
    print('Incoming population size is: ',len(pop))

    pop['t_birth'] = SFH(gx_component=gx_component, n_pop=n_pop)
    print('tBirth population size is: ',len(pop))
    pop['t_evol_GW'] = t_evol_GW(birth_time=pop.t_birth,
                                 formation_time=pop.tphys,
                                 gx_component=gx_component) # time need to evolve to today from creation of DWD
    
    # select out systems which haven't evolved yet
    print('Initial population size is: ',len(pop))
    pop = pop.loc[pop.t_evol_GW >= 0]  # takes everythng that makes a DWD
    print('After removing non-evolved systems, pop size is: ',len(pop))
    
    # evolve the systems according to Peters 64 evolution
    sep_final, ecc_final = GW_evol(pop)
    
    # import pdb
    # pdb.set_trace()

    # porb_from_sep() returns orbital period in YRS, hence the  */day_yr
    # to convert result to DAYS (std. cosmic units)
    # The *Rsun_au takes the sep in Sol Radii and converts to AU
    porb_final = porb_from_sep(m1=pop.mass_1, m2=pop.mass_2, sep=sep_final*Rsun_au)/day_yr
    
    pop['sep_final'] = sep_final
    pop['ecc_final'] = ecc_final
    pop['porb_final'] = porb_final
    
    # import pdb
    # pdb.set_trace()

    # select out those which have merged
    # print('The size of the population is: ',len(pop))
    pop = pop.loc[(pop.sep_final > 0) & (pop.porb_final*86400 < 1e7)] # why the second case? looks like a lisa cut?
    print('After removing merged systems, pop size is: ',len(pop))

    #recalculate the radii of the WDs
    #These two lines added 5/1 (natGott)
    # -----------------------------------------
    # select out those which fill roche lobes
    if pop.kstar_1.isin([10,11,12]).all() and pop.kstar_2.isin([10,11,12]).all():
        pop['rad_1'] = wdRadius(pop['mass_1'])
        pop['rad_2'] = wdRadius(pop['mass_2'])

    
    # select out those which fill roche lobes
    if pop.kstar_1.isin([13,14]).all():
        pop = pop
    else:
        pop['RL_1_final'] = R_RL(pop.mass_1/pop.mass_2, pop.sep_final)
        pop = pop.loc[pop.rad_1 <= pop.RL_1_final]
        print('kstar1 number with rstar < RL:',len(pop))
        
    if pop.kstar_2.isin([13,14]).all():
        pop = pop
    else:
        pop['RL_2_final'] = R_RL(pop.mass_2/pop.mass_1, pop.sep_final)
        pop = pop.loc[pop.rad_2 <= pop.RL_2_final]
        print('kstar2 number with rstar < RL:',len(pop))
    
    # old call: pop['f_gw_peak'] = GW_calcs.peak_gw_freq()
    pop['f_gw_peak'] = peak_gw_freq(m1=pop.mass_1*2e30, 
                                                     m2=pop.mass_2*2e30, 
                                                     ecc=pop.ecc_final, 
                                                     porb=pop.porb_final*86400)
    return pop



# ==========================================================================
# FUNCTION: LISA_galaxy
# ==========================================================================
# This version of LISA_Galaxy from KMB, June 2021
# adapted to write to file directly to avoid memory issues
# on desktop machines.
# --------------------------------------------------------------------------
# (sll) -- original kmb definition, where only path was passed in
#       and used to then construct a local filename. Now pass in
#       path and filename already combined
# def LISA_Galaxy(conv, model, m_sim_tot, kstars, path_write):
# --------------------------------------------------------------------------

def LISA_Galaxy(conv, gx_component, m_sim_tot, kstars, dat_write):
    n_pop = MC_samp.mass_weighted_number(component_mass=MC_samp.select_component_mass(gx_component),
                                         dat=conv,
                                         total_sampled_mass=m_sim_tot)
    print('\nThe number of {0} binaries in the {1} is: {2}'.format(kstars, gx_component, n_pop),flush=True)
    if np.all(conv.kstar_1 < 14):
        conv = conv.loc[conv.sep < 100]    # reduces number of binaries in pop to only those with a <= 100 Rsun
        n_pop = MC_samp.mass_weighted_number(component_mass=MC_samp.select_component_mass(gx_component),
                                             dat=conv,
                                             total_sampled_mass=m_sim_tot)
        print('filtered numbers with a < 100 Rsun {}'.format(n_pop)) 

# (sll) -- commented out Katie's filename creation and appended to path.
#       instead I pass in the entire PATH/FILENAME as "dat_write"
# define filename from a path that is passed into the function
# dat_write = '{}/gx_{}.h5'.format(path_write, kstars)

# open an hdf file to hold the data
    dat_store = pd.HDFStore(dat_write)
    
    # debugging
    myCount = 0
    
    if n_pop >= 5e5:
        pop_LISA = []
        # import pdb
        # pdb.set_trace()
        
        for ii in range(0,1000):
            # do this in chunks for memory; ### sll changed 10 to 100 in line above and below; changed 100 to 1000
            n_samp = int(n_pop/1000)
            gx_samp = gx_sample(conv, gx_component, n_samp)
            pop_LISA = gx_samp.loc[gx_samp.f_gw_peak > 1e-5].copy()     # ### Make a cut, only keep f greater than listed value
            
            # since we are going to write to the file, in the loop we need to compute the positions here too
            xGx, yGx, zGx, inc, omega, OMEGA = MC_samp.galactic_positions(gx_component=gx_component, model='McMillan', size=len(pop_LISA))
            pop_LISA['xGx'] = xGx
            pop_LISA['yGx'] = yGx
            pop_LISA['zGx'] = zGx
            pop_LISA['dist'] = ((xGx - 8.0)**2 + (yGx - 0)**2 + (zGx - 0.2)**2)**0.5

            # import pdb
            # pdb.set_trace()
            
            print('Length pop = ',len(pop_LISA),'  ii = ',ii,flush=True)
            
            myCount = myCount + len(pop_LISA)
            
            dat_store.append('LISA_population', pop_LISA)

    else:
        n_samp = n_pop # ***check this***
        gx_samp = gx_sample(conv, gx_component, n_samp)
        
        pop_LISA = gx_samp.loc[gx_samp.f_gw_peak > 1e-5].copy()    # ### Make a cut, only keep f greater than listed value
    
        # since we are writing to the file, we should compute the positions here before we write
        xGx, yGx, zGx, inc, omega, OMEGA = MC_samp.galactic_positions(gx_component=gx_component, model='McMillan', size=len(pop_LISA))
        pop_LISA['xGx'] = xGx
        pop_LISA['yGx'] = yGx
        pop_LISA['zGx'] = zGx
        pop_LISA['dist'] = ((xGx - 8.0)**2 + (yGx - 0)**2 + (zGx - 0.027)**2)**0.5  # z offset correct 0.2 -> 0.027, sll 8 Nov 2023
        myCount = myCount + len(pop_LISA)
        dat_store.append('LISA_population', pop_LISA)

    
    # make sure to close the dat file!
    dat_store.close()
    
    print('Final Count {0} binaries in {1} is: {2} '.format(kstars, gx_component, myCount),flush=True)

    #since we are writing to the file, we won't pass anything back
    return

