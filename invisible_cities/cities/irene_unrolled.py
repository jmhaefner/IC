"""
-----------------------------------------------------------------------
                                 Irene
-----------------------------------------------------------------------

From ancient Greek, Εἰρήνη: Peace.

This city finds the signal pulses within the waveforms produced by the
detector or by diomira in the case of Monte Carlo data.
This includes a number of tasks:
    - Remove the signal-derivative effect of the PMT waveforms.
    - Calibrate PMTs and produced a PMT-summed waveform.
    - Remove the baseline from the SiPM waveforms and calibrate them.
    - Apply a threshold to the PMT-summed waveform.
    - Find pulses in the PMT-summed waveform.
    - Match the time window of the PMT pulse with those in the SiPMs.
    - Build the PMap object.
"""
import tables as tb
import numpy as np
from scipy.optimize import curve_fit

from .. reco                  import tbl_functions        as tbl
from .. core.random_sampling  import NoiseSampler         as SiPMsNoiseSampler
from .. core                  import system_of_units      as units
from .. io  .run_and_event_io import run_and_event_writer
from .. io  .trigger_io       import       trigger_writer

from .. dataflow            import dataflow as fl
from .. dataflow.dataflow   import push
from .. dataflow.dataflow   import pipe
from .. dataflow.dataflow   import sink
from .. database            import load_db

from .  components import city
from .  components import print_every
from .  components import collect
from .  components import copy_mc_info
from .  components import deconv_pmt
from .  components import calibrate_pmts
from .  components import calibrate_sipms
from .  components import zero_suppress_wfs
from .  components import WfType
from .  components import wf_from_files
from .  components import get_number_of_active_pmts
from .  components import compute_and_write_pmaps
from .  components import compute_and_write_pmaps_mod
from .  components import check_nonempty_indices
from .  components import build_pmap_mod

debugging = True

def weighted_avg_and_std(values, weights, debug=debugging):
    """
    Return the weighted average and standard deviation
    values, weights -- Numpy ndarrays with the same shape.
    """
    if debug:
        print('Input values   =', values)
        print('Input weights  =', weights)
        print('Sum of weights =', np.sum(weights))

    if np.sum(weights) > 0:
        average = np.average(values, weights=weights)
        variance = np.average((values-average)**2, weights=weights)
        if debug: print('Avg, var, std =', average, variance, variance**0.5)
        std = np.sqrt(max(0.,variance))
    else:
        if debug:
            for i in range(10): print('! WARNING ! WARNING ! WARNING ! WARNING !')
            print('WARNING - NOT ENOUGH WEIGHTS TO CALCULATE.')
            for i in range(10): print('! WARNING ! WARNING ! WARNING ! WARNING !')
        average = np.mean(values)
        std = 0.

    return average, std 

def exp_dt(sig):
    """
    expected drift distnce  as a function of sigma of the S2 signal (Zrms from kdst 2018 version)
    parameters from the ad-hoc fit
    """
    p = [-10.00667732,  52.1855012,   12.68195726,  58.66322846, -20.11819297]
    dt = p[0] * sig**4 + p[1] * sig**3 + p[2]*sig**2 + p[3]*sig + p[4]
    return dt

def gauss(x, A, m, v):
    if v <= 0:
        return 1e10
    return A*np.exp(-(x-m)**2 / (2*v))

def offset_gauss(x, A, m, v, C):
    if v <= 0:
        return 1e10
    return C+A*np.exp(-(x-m)**2 / (2*v))

@city
def irene_unrolled(files_in, file_out, compression, event_range, print_mod, detector_db, run_number,
          n_baseline, n_mau, thr_mau, thr_sipm, thr_sipm_type,
          s1_lmin, s1_lmax, s1_tmin, s1_tmax, s1_rebin_stride, s1_stride, thr_csum_s1,
          s2_lmin, s2_lmax, s2_tmin, s2_tmax, s2_rebin_stride, s2_stride, thr_csum_s2, thr_sipm_s2,
          pmt_samp_wid=25*units.ns, sipm_samp_wid=1*units.mus):
    if   thr_sipm_type.lower() == "common":
        # In this case, the threshold is a value in pes
        sipm_thr = thr_sipm

    elif thr_sipm_type.lower() == "individual":
        # In this case, the threshold is a percentual value
        noise_sampler = SiPMsNoiseSampler(detector_db, run_number)
        sipm_thr      = noise_sampler.compute_thresholds(thr_sipm)

    else:
        raise ValueError(f"Unrecognized thr type: {thr_sipm_type}. "
                          "Only valid options are 'common' and 'individual'")

    #### Define data transformations

    # Raw WaveForm to Corrected WaveForm
    rwf_to_cwf       = fl.map(deconv_pmt(detector_db, run_number, n_baseline),
                              args = "pmt",
                              out  = "cwf")

    # Corrected WaveForm to Calibrated Corrected WaveForm
    cwf_to_ccwf      = fl.map(calibrate_pmts(detector_db, run_number, n_mau, thr_mau),
                              args = "cwf",
                              out  = ("ccwfs", "ccwfs_mau", "cwf_sum", "cwf_sum_mau"))

    # Find where waveform is above threshold
    zero_suppress    = fl.map(zero_suppress_wfs(thr_csum_s1, thr_csum_s2),
                              args = ("cwf_sum", "cwf_sum_mau"),
                              out  = ("s1_indices", "s2_indices"))

    # Remove baseline and calibrate SiPMs
    sipm_rwf_to_cal  = fl.map(calibrate_sipms(detector_db, run_number, sipm_thr),
                              item = "sipm")

    event_count_in  = fl.spy_count()
    event_count_out = fl.spy_count()

    evtnum_collect  = collect()

    with tb.open_file(file_out, "w", filters = tbl.filters(compression)) as h5out:

        # Define writers...
        write_event_info_   = run_and_event_writer(h5out)
        write_trigger_info_ = trigger_writer      (h5out, get_number_of_active_pmts(detector_db, run_number))

        # ... and make them sinks

        write_event_info   = sink(write_event_info_  , args=(   "run_number",     "event_number", "timestamp"   ))
        write_trigger_info = sink(write_trigger_info_, args=( "trigger_type", "trigger_channels"                ))

        rwfs = wf_from_files(files_in, WfType.rwf)
        rwfs_test = wf_from_files(files_in, WfType.rwf)
        all_pmaps = []
        nrwfs = len(list(rwfs_test))
        print('Total rwfs:', nrwfs)
        import time
        t0 = time.time() 

        window_min = 100 # Window radius 100 samples = 2500 ns radius = 5 us width
        window_max = 500 # Window radius 500 samples = 12500 ns radius = 25 us width
        window_step = 5 # Increase by 5 samples = 125 ns each step
        window_ranges = np.arange(window_min, window_max, window_step)
        nwindow = len(window_ranges)
        event_numbers = [ k for k in range(nrwfs) ]
        k = -1

        max_events = 1
        nevents = 0

        events_prewindow_sums = [] # Sum of all signals in window from -75 to -25 us
        for rwf in rwfs:

            if nevents >= max_events and max_events > 0:
                continue

            nevents += 1

            k += 1
            #if len(all_pmaps) > 0:
            #    continue

            remaining = 999999.9
            if len(all_pmaps) > 0:
                remaining = (nrwfs - len(all_pmaps)) * (time.time() - t0) / len(all_pmaps)
            print(len(all_pmaps), round(time.time() - t0, 1), ', Est Remaining:', round(remaining, 1))

            try:
                cwf = deconv_pmt(detector_db, run_number, n_baseline)(rwf['pmt'])
                ccwfs, ccwfs_mau, cwf_sum, cwf_sum_mau = calibrate_pmts(detector_db, run_number, n_mau, thr_mau)(cwf)
                s1_indices, s2_indices = zero_suppress_wfs(thr_csum_s1, thr_csum_s2)(cwf_sum, cwf_sum_mau)
                indices_pass = check_nonempty_indices(s1_indices, s2_indices)            
                cal_sipms = calibrate_sipms(detector_db, run_number, sipm_thr)
                sipm = cal_sipms(rwf['sipm'])
                pmaps = []
                for window_r in window_ranges:
                    pmap = build_pmap_mod(detector_db, run_number, pmt_samp_wid, sipm_samp_wid,
                                                 s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
                                                 s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin, thr_sipm_s2, window_r = window_r)
                    pmaps.append(pmap( ccwfs, s1_indices, s2_indices, sipm )  )

                ccwf_sum = [ sum(ccwfs[:,i]) for i in range(len(ccwfs[0])) ]
                signal_peak = np.argmax(ccwf_sum)

                start_sample = int(signal_peak - 75000 * (1 / 25.0))
                end_sample   = int(signal_peak - 25000 * (1 / 25.0))
                prewindow_sum    = np.sum(ccwf_sum[start_sample:end_sample])
                events_prewindow_sums.append(prewindow_sum)

                all_pmaps.append(pmaps)
            except:
                print('EVENT FAILED')
                event_numbers.remove(k)
                continue

        #from IPython import embed
        #embed()
        #quit()

        events_window_maxsipm = []
        events_window_Zrms = []
        events_window_Zgauss = []
        events_window_Wrms = []
        events_window_Wgauss = []
        events_window_r2 = []
        events_window_GaussOffset = []

        for pmap in all_pmaps:
            window_maxsipm = []
            window_Zrms = []
            window_Zgauss = []
            window_Wrms = []
            window_Wgauss = []
            window_r2 = []
            window_GaussOffset = []         

            for w in range(nwindow): 

                signals = pmap[w].s2s[0].pmts.sum_over_sensors
                times   = pmap[w].s2s[0].times

                # Calculate max sipm
                try:
                    imax = np.argmax(pmap[w].s2s[0].sipms.sum_over_times)
                    smax = pmap[w].s2s[0].sipms.ids[imax]
                except:
                    smax = -1
                window_maxsipm.append(smax)

                # Calculate the Zrms
                try:
                    if debugging:
                        print('CURRENT WINDOW', w, '/', nwindow)
                    mean, stdev = weighted_avg_and_std(times, signals)
                    Zrms = exp_dt(stdev/1000)
                    Wrms = stdev
                except:
                    Wrms = -1
                    Zrms = -1
                window_Zrms.append(Zrms)
                window_Wrms.append(Wrms)

                # Calculate the Zgauss
                try:
                    xdata, ydata = times, signals
                    mean0, _ = weighted_avg_and_std(xdata, ydata)
                    amp0 = np.max(ydata)

                    # Get the stdev guess from the less edge dependent FWHM
                    above_half_max = ydata > amp0 / 2
                    start = np.argmax(above_half_max)
                    end = len(above_half_max)-np.argmax(above_half_max[::-1])-1
                    FWHM = xdata[end] - xdata[start]
                    sigma0 = FWHM / (2*np.sqrt(2*np.log(2)))
                    var0 = sigma0**2
                    C0 = 0

                    popt,  pcov  = curve_fit(gauss, xdata, ydata, p0 = (amp0, mean0, var0))
                    popt1, pcov1 = curve_fit(offset_gauss, xdata, ydata, p0 = (amp0, mean0, var0, C0))
                    amp, mean, var = popt
                    ampoff, meanoff, varoff, coff = popt1
                    stdev = var**0.5
                    Wgauss = stdev
                    Zgauss = exp_dt(stdev/1000)

                    residuals = ydata - gauss(xdata, *popt)
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((ydata-np.mean(ydata))**2)
                    r2 = 1 - (ss_res / ss_tot)

                except:
                    Wgauss = -1
                    Zgauss = -1
                    coff = 0
                    r2 = 0
                window_Zgauss.append(Zgauss)
                window_Wgauss.append(Wgauss)
                window_GaussOffset.append(coff)
                window_r2.append(r2)

            # print('\nCHECKING PMAP...')
            # print('Getting signals and times')
            # signals = pmap[200].s2s[0].pmts.sum_over_sensors
            # times   = pmap[200].s2s[0].times
            # print('Getting mean, stdev')
            # mean, stdev = weighted_avg_and_std(times, signals, True)
            # print('\tstdev =', stdev)
            # print('Getting Zrms')
            # Zrms = exp_dt(stdev/1000)
            # print('\tZrms =', Zrms)
            # print('\tEnergy =', pmap[200].s2s[0].total_energy)
            events_window_Wrms.append(window_Wrms)
            events_window_Zrms.append(window_Zrms)
            events_window_Zgauss.append(window_Zgauss)
            events_window_Wgauss.append(window_Wgauss)
            events_window_maxsipm.append(window_maxsipm)
            events_window_GaussOffset.append(window_GaussOffset)
            events_window_r2.append(window_r2)

        DataSiPM    = load_db.DataSiPM(detector_db, run_number = run_number) 
        events_window_maxX = [ [ DataSiPM.X[max(i,0)] for i in window_maxsipm ] for window_maxsipm in events_window_maxsipm ]
        events_window_maxY = [ [ DataSiPM.Y[max(i,0)] for i in window_maxsipm ] for window_maxsipm in events_window_maxsipm ]
        events_window_rms = [ [ all_pmaps[i][w].s2s[0].rms for w in range(nwindow) ] for i in range(len(all_pmaps)) ]
        events_window_width = [ [ all_pmaps[i][w].s2s[0].width for w in range(nwindow) ] for i in range(len(all_pmaps)) ]
        events_window_energy = [ [ all_pmaps[i][w].s2s[0].total_energy for w in range(nwindow) ] for i in range(len(all_pmaps)) ] 
        events_window_charge = [ [ all_pmaps[i][w].s2s[0].total_charge for w in range(nwindow) ] for i in range(len(all_pmaps)) ]
        events_sumwf = [ all_pmaps[i][-1].s2s[0].pmts.sum_over_sensors.tolist() for i in range(len(all_pmaps)) ] 
        events_times = [ all_pmaps[i][-1].s2s[0].times.tolist() for i in range(len(all_pmaps)) ] 

        if debugging:
            print('Wrms vectors:')
            print(events_window_Wrms)
            print('EVENT ENERGIES:')
            print([ window_energy[-1] for window_energy in events_window_energy ])

        window_out_name = file_out.replace('h5', 'json')
        outdict  = {'events_window_maxX' : events_window_maxX,
                    'events_window_maxY' : events_window_maxY,
                    'events_window_rms' : events_window_rms,
                    'events_window_width' : events_window_width,
                    'events_window_energy' : events_window_energy,
                    'events_window_charge' : events_window_charge,
                    'events_window_Zrms' : events_window_Zrms,
                    'events_window_Zgauss' : events_window_Zgauss,
                    'events_window_Wrms' : events_window_Wrms,
                    'events_window_Wgauss' : events_window_Wgauss,
                    'events_window_r2' : events_window_r2,
                    'events_window_GaussOffset' : events_window_GaussOffset,
                    'events_sumwf' : events_sumwf,
                    'events_times' : events_times,
                    'event_numbers': event_numbers,
                    'events_prewindow_sums': events_prewindow_sums}

        import json
        with open(window_out_name, 'w') as fileout:
            json.dump(outdict, fileout)

        return True
