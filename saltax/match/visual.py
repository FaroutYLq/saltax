import matplotlib.pyplot as plt
import numpy as np

def plot_wf(ind, st_salt, st_simu, st_data, runid, matched_simu, 
            event_ext_window_ns=2.4e6, 
            s1_ext_window_samples=25, 
            s2_ext_window_samples=100):
    """
    Plot waveforms for a single event, including full event waveform and zoomed-in S1 and S2 waveforms.
    :param ind: index of the event in the matched_simu dataframe
    :param st_salt: saltax context for salt mode
    :param st_simu: saltax context for simu mode
    :param st_data: saltax context for data mode
    :param runid: runid of the event, example: '066666'
    :param matched_simu: simu event_info already matched to truth. (Not necessarily matched to data!)
    :param event_ext_window_ns: time window in ns to plot around the event, default 2.4e6 ns = 2.4 ms
    :param s1_ext_window_samples: time window in samples to plot around S1, default 25 samples
    :param s2_ext_window_samples: time window in samples to plot around S2, default 100 samples
    """
    print("Loading peaks and lone_hits for run %s event %s"%(runid, ind))

    # Get time ranges in indices for events, S1 and S2
    extended_simu_event_timerange_ns = (matched_simu['s1_time'][ind]-event_ext_window_ns, 
                                        matched_simu['s2_endtime'][ind]+event_ext_window_ns)
    matched_simu_s1_timerange_i = (int((matched_simu['s1_time'][ind] - extended_simu_event_timerange_ns[0])/10), 
                                   int((matched_simu['s1_endtime'][ind]-extended_simu_event_timerange_ns[0])/10))
    matched_simu_s2_timerange_i = (int((matched_simu['s2_time'][ind] - extended_simu_event_timerange_ns[0])/10), 
                                   int((matched_simu['s2_endtime'][ind]-extended_simu_event_timerange_ns[0])/10))
    

    # Get peaks and lone hits for the event
    peaks_salt_selected = st_salt.get_array(runid, "peaks", 
                                            time_range=extended_simu_event_timerange_ns, 
                                            progress_bar=False)
    peaks_simu_selected = st_simu.get_array(runid, "peaks", 
                                            time_range=extended_simu_event_timerange_ns, 
                                            progress_bar=False)
    peaks_data_selected = st_data.get_array(runid, "peaks", 
                                            time_range=extended_simu_event_timerange_ns, 
                                            progress_bar=False)
    lhs_salt_selected = st_salt.get_array(runid, "lone_hits", 
                                          time_range=extended_simu_event_timerange_ns, 
                                          progress_bar=False)
    lhs_simu_selected = st_simu.get_array(runid, "lone_hits", 
                                          time_range=extended_simu_event_timerange_ns, 
                                          progress_bar=False)
    lhs_data_selected = st_data.get_array(runid, "lone_hits", 
                                          time_range=extended_simu_event_timerange_ns, 
                                          progress_bar=False)


    # Get waveforms for the event
    print("Building waveforms...")
    total_length = int((extended_simu_event_timerange_ns[1] - extended_simu_event_timerange_ns[0])/10)
    to_pes = st_data.get_single_plugin(runid, 'peaklets').to_pe
    
    wf_salt_s1 = np.zeros(total_length)
    wf_simu_s1 = np.zeros(total_length)
    wf_salt_s2 = np.zeros(total_length)
    wf_simu_s2 = np.zeros(total_length)
    wf_salt_others = np.zeros(total_length)
    wf_simu_others = np.zeros(total_length)
    wf_data = np.zeros(total_length)
    
    if len(peaks_salt_selected):
        for p in peaks_salt_selected:
            start_i = int((p['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            length = p['length']
            dt = p['dt']
            if p['type'] == 1:
                for i in range(length):
                    wf_salt_s1[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
            elif p['type'] == 2:
                for i in range(length):
                    wf_salt_s2[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
            else:
                for i in range(length):
                    wf_salt_others[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
    if len(lhs_salt_selected):
        for l in lhs_salt_selected:
            time_i = int((l['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            amp = l['area'] * to_pes[l['channel']]
            wf_salt_others[time_i] += amp/10
    
    if len(peaks_simu_selected):
        for p in peaks_simu_selected:
            start_i = int((p['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            length = p['length']
            dt = p['dt']
            if p['type'] == 1:
                for i in range(length):
                    wf_simu_s1[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
            elif p['type'] == 2:
                for i in range(length):
                    wf_simu_s2[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
            else:
                for i in range(length):
                    wf_simu_others[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
    if len(lhs_simu_selected):
        for l in lhs_simu_selected:
            time_i = int((l['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            amp = l['area'] * to_pes[l['channel']]
            wf_simu_others[time_i] += amp/10
    
    if len(peaks_data_selected):
        for p in peaks_data_selected:
            start_i = int((p['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            length = p['length']
            dt = p['dt']
            for i in range(length):
                wf_data[start_i+i*int(dt/10):start_i+(i+1)*int(dt/10)] = p['data'][i] / dt * 10
    if len(lhs_data_selected):
        for l in lhs_data_selected:
            time_i = int((l['time'] - int(extended_simu_event_timerange_ns[0]))/10)
            amp = l['area'] * to_pes[l['channel']]
            wf_data[time_i] += amp/10
    
    # Plot full event waveform
    print("Plotting waveforms for the whole event...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), dpi=150)
    ax1.plot(wf_data, label='Data', color='k', alpha=0.2)
    ax1.plot(wf_simu_s1, label='Simu S1', color='tab:blue')
    ax1.plot(wf_simu_s2, label='Simu S2', color='tab:orange')
    ax1.plot(wf_simu_others, label='Simu Others', color='tab:green')
    ax1.axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    ax1.axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    ax1.legend()
    ax1.set_title("Run %s Event %s: Simu CS1=%sPE, Simu CS2=%sPE"%(runid, ind, 
                                                        int(10*matched_simu['cs1'][ind])/10, 
                                                        int(10*matched_simu['cs2'][ind])/10))
    
    ax2.plot(wf_data, label='Data', color='k', alpha=0.2)
    ax2.plot(wf_salt_s1, label='Salt S1', color='b')
    ax2.plot(wf_salt_s2, label='Salt S2', color='r')
    ax2.plot(wf_salt_others, label='Salt Others', color='g')
    ax2.axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    ax2.axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    ax2.legend()
    
    ax1.set_ylabel("Amplitude [PE/10ns]")
    ax2.set_xlabel("Time [10ns]")
    ax2.set_ylabel("Amplitude [PE/10ns]")
    fig.show()
    
    # Zoom into S1 and S2 waveforms
    print("Zooming into S1 and S2 respectively...")
    fig, axs = plt.subplots(2, 2, figsize=(15, 8), dpi=150)
    axs[0,0].plot(wf_data, label='Data', color='k', alpha=0.2)
    axs[0,0].plot(wf_simu_s1, label='Simu S1', color='tab:blue')
    axs[0,0].plot(wf_simu_s2, label='Simu S2', color='tab:orange')
    axs[0,0].plot(wf_simu_others, label='Simu Others', color='tab:green')
    axs[0,0].axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    axs[0,0].axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    axs[0,0].set_xlim(matched_simu_s1_timerange_i[0]-s1_ext_window_samples, 
                      matched_simu_s1_timerange_i[1]+s1_ext_window_samples)
    axs[0,0].set_ylabel("Amplitude [PE/10ns]")
    axs[0,0].legend()
    
    axs[0,1].plot(wf_data, label='Data', color='k', alpha=0.2)
    axs[0,1].plot(wf_simu_s1, label='Simu S1', color='tab:blue')
    axs[0,1].plot(wf_simu_s2, label='Simu S2', color='tab:orange')
    axs[0,1].plot(wf_simu_others, label='Simu Others', color='tab:green')
    axs[0,1].axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    axs[0,1].axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    axs[0,1].set_xlim(matched_simu_s2_timerange_i[0]-s2_ext_window_samples, 
                      matched_simu_s2_timerange_i[1]+s2_ext_window_samples)
    axs[0,1].legend()
    
    axs[1,0].plot(wf_data, label='Data', color='k', alpha=0.2)
    axs[1,0].plot(wf_salt_s1, label='Salt S1', color='b')
    axs[1,0].plot(wf_salt_s2, label='Salt S2', color='r')
    axs[1,0].plot(wf_salt_others, label='Salt Others', color='g')
    axs[1,0].axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    axs[1,0].axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    axs[1,0].set_xlim(matched_simu_s1_timerange_i[0]-s1_ext_window_samples, 
                      matched_simu_s1_timerange_i[1]+s1_ext_window_samples)
    axs[1,0].set_xlabel("Time [10ns]")
    axs[1,0].set_ylabel("Amplitude [PE/10ns]")
    axs[1,0].legend()
    
    axs[1,1].plot(wf_data, label='Data', color='k', alpha=0.2)
    axs[1,1].plot(wf_salt_s1, label='Salt S1', color='b')
    axs[1,1].plot(wf_salt_s2, label='Salt S2', color='r')
    axs[1,1].plot(wf_salt_others, label='Salt Others', color='g')
    axs[1,1].axvspan(matched_simu_s1_timerange_i[0],matched_simu_s1_timerange_i[1], color='tab:blue', alpha=0.1, label='Simu S1 Range')
    axs[1,1].axvspan(matched_simu_s2_timerange_i[0],matched_simu_s2_timerange_i[1], color='tab:orange', alpha=0.1, label='Simu S2 Range')
    axs[1,1].set_xlim(matched_simu_s2_timerange_i[0]-s2_ext_window_samples, 
                      matched_simu_s2_timerange_i[1]+s2_ext_window_samples)
    axs[1,1].set_xlabel("Time [10ns]")
    axs[1,1].legend()
    fig.show()
    