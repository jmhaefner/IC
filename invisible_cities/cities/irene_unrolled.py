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

def weighted_avg_and_std(values, weights, debug=False):
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
        if debug: print('Average, variance =', average, variance)
        std = np.sqrt(max(0.,variance))
    else:
        if debug:
            print('WARNING - NOT ENOUGH WEIGHTS TO CALCULATE.')
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
        for rwf in rwfs:
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
                window_ranges = np.arange(100, 500, 1)
                pmaps = []
                for window_r in window_ranges:
                    pmap = build_pmap_mod(detector_db, run_number, pmt_samp_wid, sipm_samp_wid,
                                                 s1_lmax, s1_lmin, s1_rebin_stride, s1_stride, s1_tmax, s1_tmin,
                                                 s2_lmax, s2_lmin, s2_rebin_stride, s2_stride, s2_tmax, s2_tmin, thr_sipm_s2, window_r = window_r)
                    pmaps.append(pmap( ccwfs, s1_indices, s2_indices, sipm )  )

                all_pmaps.append(pmaps)
            except:
                print('EVENT FAILED')
                continue


        events_window_maxsipm = []
        events_window_Zrms = []

        for pmap in all_pmaps:
            window_maxsipm = []
            window_Zrms = []
            for w in range(400): 
                try:
                    # Calculate for max sipm
                    imax = np.argmax(pmap[w].s2s[0].sipms.sum_over_times)
                    smax = pmap[w].s2s[0].sipms.ids[imax]
                    window_maxsipm.append(smax)
                except:
                    window_maxsipm.append(-1)

                try:
                    # Calculate the Zrms
                    signals = pmap[w].s2s[0].pmts.sum_over_sensors
                    times   = pmap[w].s2s[0].times
                    mean, stdev = weighted_avg_and_std(times, signals)
                    Zrms = exp_dt(stdev/1000)
                    window_Zrms.append(Zrms)
                except:
                    Zrms = -1
                    window_Zrms.append(Zrms)

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
            events_window_Zrms.append(window_Zrms)
            events_window_maxsipm.append(window_maxsipm)

        DataSiPM    = load_db.DataSiPM(detector_db, run_number = run_number) 
        events_window_maxX = [ [ DataSiPM.X[max(i,0)] for i in window_maxsipm ] for window_maxsipm in events_window_maxsipm ]
        events_window_maxY = [ [ DataSiPM.Y[max(i,0)] for i in window_maxsipm ] for window_maxsipm in events_window_maxsipm ]
        events_window_rms = [ [ all_pmaps[i][w].s2s[0].rms for w in range(400) ] for i in range(len(all_pmaps)) ]
        events_window_width = [ [ all_pmaps[i][w].s2s[0].width for w in range(400) ] for i in range(len(all_pmaps)) ]
        events_window_energy = [ [ all_pmaps[i][w].s2s[0].total_energy for w in range(400) ] for i in range(len(all_pmaps)) ] 
        events_window_charge = [ [ all_pmaps[i][w].s2s[0].total_charge for w in range(400) ] for i in range(len(all_pmaps)) ]
        events_window_sumwf = [ [ all_pmaps[i][w].s2s[0].pmts.sum_over_sensors for w in range(400) ] for i in range(len(all_pmaps)) ] 

        # print('Entering interactive session...')
        # from IPython import embed
        # embed()
        # quit()

        fileout = open('variable_window_out.txt', 'w')
        fileout.write('events_window_maxX='+str(events_window_maxX)+'\n')
        fileout.write('events_window_maxY='+str(events_window_maxY)+'\n')
        fileout.write('events_window_rms='+str(events_window_rms)+'\n')
        fileout.write('events_window_width='+str(events_window_width)+'\n')
        fileout.write('events_window_energy='+str(events_window_energy)+'\n')
        fileout.write('events_window_charge='+str(events_window_charge)+'\n')
        fileout.write('events_window_Zrms='+str(events_window_Zrms)+'\n')
        fileout.write('events_window_sumwf='+str(events_window_sumwf)+'\n')
        fileout.close()

        return True
