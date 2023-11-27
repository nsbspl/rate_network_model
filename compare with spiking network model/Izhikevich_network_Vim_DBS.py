import nest
from scipy import io
import scipy.io.matlab
import scipy.signal
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

def setNestOpts():
    nest.ResetKernel()
    nest.resolution = 0.1
    nest.local_num_threads = 7

##
# Stat
##

DBS_freq = 100

def rasterPlot(ev):
    srcs = ev['senders']
    srcs_ = set(ev['senders'])
    m = min(srcs_)
    M = max(srcs_)
    tims = ev['times']
    for s in srcs_:
        c = f'#2757{int((s-m)/(M-m)*255):02X}'
        plt.vlines(tims[srcs == s], s-m-1, s-m, c)


def static_syn(w):
    return {
      'weight': w,
      'delay': 1 #1
    }

def plastic_syn(w):
    return {
        'synapse_model': 'tsodyks2_synapse',
        'U': 0.5, # 0.5; real | Parameter determining the increase in u
        'u': 0.35, # 0.4; The probability of release (U_se) [0,1], default=0.5
        'tau_fac': 1e2, # 1e2; ms | Time constant for facilitation
        'tau_rec': 1e2, # 1e2; ms | Time constant for depression
        'x': 1, #1, real | Initial fraction of syn vesicles ready
        'weight': w*2, #w*2
        'delay': 1,
    }

def spikeStats(ev, time_len, n=1): # compute mean FR of a sample neuron in a neuronal group
    srcs = ev['senders']
    means = np.bincount(srcs)[-n:]/time_len
    return np.mean(means)

def runSynapseTests(freq, n_pulses, static=False,):
    setNestOpts()

    tstart = 200
    isi = 1000/freq
    tdur = isi * n_pulses
    tstop = tstart + tdur + isi
    
    DBS_ts = np.arange(tstart,tstop,isi)
    
    DBS = nest.Create(
        "spike_generator",
        params={
            "precise_times": True,
            "spike_times": DBS_ts
        }
    )

    prt = nest.Create("parrot_neuron_ps",1)
    nrn = nest.Create("izhikevich", 1)

    vm = nest.Create('multimeter', 1, {
        'record_from': ['V_m'],
        'interval': 0.1
    })

    nest.Connect(vm, nrn)

    syn = None
    if static:
        syn = static_syn
    else:
        syn = plastic_syn
    
    nest.Connect(DBS, prt)
    nest.Connect(
        prt,
        nrn,
        syn_spec=syn(1)
    )
    
    nest.Simulate(int(tstop+isi))

    v = np.array(vm.get()['events']['V_m'])
    start_idx = int(tstart / 0.1)

    return v[start_idx:]



def testSynapseAtFreqs(static):
    freqs = [5,10,30,50,100]
    n_f = len(freqs)
    for i,f in enumerate(freqs):
        v = runSynapseTests(f,20, static)
        plt.subplot(n_f*100+10+i+1)
        plt.title(f"{f} Hz")
        plt.plot(v)

def runSim(weights, exc_static=False, inh_static=False, label="NA"):
    ###
    # Basic Definitions
    ###

    tstop = 5

    setNestOpts()

    ##
    # DBS setup
    ##

    t_before_DBS = 2 #s; the initial 2s is DBS-OFF.
                       #DBS start at 2s, may (or may not) last till the end of tstop = 5


    #IDBS = scipy.io.matlab.loadmat(f"./I_DBS_fq/I_DBS_{DBS_freq}Hz.mat")
    IDBS = scipy.io.matlab.loadmat(f"/Users/yupengtian/miniconda3/code_sim_Izhi_23_Jun/I_DBS_fq/I_DBS_{DBS_freq}Hz.mat")
    IDBS = IDBS[list(IDBS.keys())[-1]].T[0]/40
    ts = np.array([(i+1)*0.1 for i, _ in enumerate(IDBS)])+ t_before_DBS*1000 # IDBS time series


    ##
    # Neuronal groups
    ##

    pop_scale = 1

    vim_n = int(20 * pop_scale) # num of Vim neurons
    cer_n = int(100 * pop_scale) # num of Cerebellum neurons
    trn_n = int(40 * pop_scale) # num of TRN neurons

    # Bias currents
    vim_I = 0.45 
    cer_I = 6.5 
    trn_I = 0.45 



    ##
    # Define synapses
    ##

    wee = weights[0]
    wie = weights[1]
    wei = weights[2]
    wii = weights[3]

    def syn(w, excOrInh):
        isStatic = None
        if excOrInh == "exc":
            isStatic = exc_static
        else:
            isStatic = inh_static
        if isStatic:
            return static_syn(w)
        else:
            return plastic_syn(w)

    
    conn_prob_vim = 0.1/pop_scale # sparse conn prob
    conn_prob_cer = 0.1/pop_scale
    conn_prob_trn = 0.25/pop_scale #0.25
    
    
    
    conn_dict_vim = {"rule": "fixed_indegree","indegree": int(conn_prob_vim*vim_n)}
    conn_dict_cer = {"rule": "fixed_indegree","indegree": int(conn_prob_cer*cer_n)}
    conn_dict_trn = {"rule": "fixed_indegree","indegree": int(conn_prob_trn*trn_n)}

    #def sparse_conn(src, dst, conn_dict, w, excOrInh):
        #nest.Connect(src, dst, conn_dict, syn_spec=syn(w, excOrInh))
        
    all_syn_parrots = []
    
    def sparse_conn(src, dst, conn_dict, w, excOrInh):
        n_src = len(src)
        n_nodes = conn_dict["indegree"]
        for d in dst:
            src_i = np.random.choice(np.arange(n_src),n_nodes,replace=False)
            src_i.sort()
            src_ = src[src_i]
            prts = nest.Create("parrot_neuron",n_nodes)
            all_syn_parrots.extend(prts.tolist())
            nest.Connect(src_, prts, 'one_to_one') 
            nest.Connect(prts, d, syn_spec=syn(w, excOrInh))
  


    ##
    # PSTH computing method setup
    ##

    bin_w = 1 # ms; for computing PSTH with constant kernel
    bin_edge = np.arange(0,tstop*1000+bin_w,bin_w) # the bins for computing PSTH with constant kernel

    sigmas = {
        '5': 8.1,
        '10': 5.7,
        '20': 4.4,
        '30': 3.9,
        '50': 3.5,
        '100': 53.2,
        '200': 31.9
    }

    sigma = sigmas[str(DBS_freq)]
    s3 = 3*sigma
    gx = np.arange(-s3, s3, bin_w)
    PSTH_window = np.exp(-((gx/sigma)**2)/2)/((2*np.pi*sigma**2)**0.5) #PSTH Gaussian kernel

    time_vec = np.arange(1,tstop*1000+bin_w,bin_w) # time_vec for plot

    ###
    # Simulation
    ###

    ##
    # Izhi model
    ##

    vim = nest.Create("izhikevich", vim_n, params={
        "V_th": nest.random.normal(30, vim_I/2),
        "a": 0.02,
        "b": 0.25,
        "c": -65,
        "d": 0.05
    })

    vim_bias = nest.Create("noise_generator", params={
        "mean": vim_I,
        "std": vim_I/2 #vim_I/2
    })

    nest.Connect(vim_bias, vim)


    cer = nest.Create("izhikevich", cer_n, params={
        "V_th": nest.random.normal(30, cer_I/2),
        "a": 0.02,
        "b": 0.2, #0.2
        "c": -65,
        "d": 8
    })

    cer_bias = nest.Create("noise_generator", params={
        "mean": cer_I,
        "std": cer_I/2
    })

    nest.Connect(cer_bias, cer)


    trn = nest.Create("izhikevich", trn_n, params={
        "V_th": nest.random.normal(30, abs(trn_I)/2),
        "a": 0.015, 
        "b": 0.25, 
        "c": -65,
        "d": 2.05 
    })

    trn_bias = nest.Create("noise_generator", params={
        "mean": trn_I,
        "std": abs(trn_I)/2
    })

    nest.Connect(trn_bias, trn)

    ##
    # Connection specifications (sparse)
    ##

    # 1. Vim
    # Vim->Vim
    sparse_conn(vim, vim, conn_dict_vim, wee, "exc")
    # Vim->Cer
    sparse_conn(vim, cer, conn_dict_vim, wee, "exc")
    # Vim->TRN
    sparse_conn(vim, trn, conn_dict_vim, wie, "exc")

    # 2. Cer
    # Cer->Cer
    sparse_conn(cer, cer, conn_dict_cer, wee, "exc")
    # Cer->VIM
    sparse_conn(cer, vim, conn_dict_cer, wee, "exc")
    # Cer->TRN
    sparse_conn(cer, trn, conn_dict_cer, wie, "exc")

    # 3. TRN
    # TRN->TRN
    sparse_conn(trn, trn, conn_dict_trn, wii, "inh")
    # TRN->VIM
    sparse_conn(trn, vim, conn_dict_trn, wei, "inh")
    # TRN->Cer
    sparse_conn(trn, cer, conn_dict_trn, wei, "inh")


    ##
    # Induce DBS 
    ##

    # The DBS-induced inputs afferent to Vim are obtained by the current injection (same inputs as the Vim-network rate model)
    DBS_aff = nest.Create("step_current_generator", params={
        "amplitude_times": ts,
        "amplitude_values": IDBS
    })

    nest.Connect(DBS_aff, vim)
    
    # The DBS activation of axons efferent from Vim
    
    DBS_tstart = 2000
    DBS_isi = 1000/DBS_freq
    DBS_tdur = DBS_isi * (3*DBS_freq)
    DBS_tstop = DBS_tstart + DBS_tdur
    
    DBS_ts = np.arange(DBS_tstart,DBS_tstop,DBS_isi)
    
    DBS_act = nest.Create(
        "spike_generator",
        params={
            "precise_times": True,
            "spike_times": DBS_ts
        }
    )

    prt = nest.Create("parrot_neuron_ps",1)
    nest.Connect(DBS_act, prt)    
    
    ###
    all_syn_parrots = nest.NodeCollection(all_syn_parrots)
    #print("all_syn_parrots=", all_syn_parrots)
    
    ## Percent of VIM DBS Activates
    DBS_percent = 0.8

    ## Efficacy of efferent synapses activated per DBS pulse
    DBS_eff = 0.3 

    DBS_tgts = vim[:int(vim_n*DBS_percent)]
    #print("DBS_tgts=", DBS_tgts)
    
    tgt_effs = nest.NodeCollection(
        sorted(set(
            nest.GetConnections(source=DBS_tgts,
                                target=all_syn_parrots).get('target')
        ))
    )
    #print("tgt_effs=", tgt_effs)
    #print(int(len(tgt_effs)*DBS_eff))
    
    tgt_effs = tgt_effs[:int(len(tgt_effs)*DBS_eff)]
    #print("tgt_effs=", tgt_effs)

    nest.Connect(DBS_act, tgt_effs)    
    

    ##
    # Recording
    ##

    vim_spikes = nest.Create("spike_recorder")
    nest.Connect(vim, vim_spikes)

    cer_spikes = nest.Create("spike_recorder")
    nest.Connect(cer, cer_spikes)

    trn_spikes = nest.Create("spike_recorder")
    nest.Connect(trn, trn_spikes)


    ##
    # Simulation part 1 --> compute the baseline firing rates
    #                       from data in the initial 2s (DBS-OFF)
    ##

    nest.Simulate(t_before_DBS*1000)

    # baseline spikes
    vim_ev_b = vim_spikes.get("events")
    cer_ev_b = cer_spikes.get("events")
    trn_ev_b = trn_spikes.get("events")

    print(f"VIM baseline ({label}): {spikeStats(vim_ev_b,t_before_DBS, vim_n)} Hz")
    print(f"CER baseline ({label}): {spikeStats(cer_ev_b,t_before_DBS, cer_n)} Hz")
    print(f"TRN baseline ({label}): {spikeStats(trn_ev_b,t_before_DBS, trn_n)} Hz")


    ##
    # Simulation part 2 --> complete the remaining simulation
    ##

    nest.Simulate((tstop-t_before_DBS)*1000)

    vim_ev = vim_spikes.get("events")
    cer_ev = cer_spikes.get("events")
    trn_ev = trn_spikes.get("events")

    vim_ts = vim_ev['times']
    vim_hist = np.histogram(vim_ts,bins=bin_edge)[0]*1000/vim_n
    kernel_PSTH_vim_ori = scipy.signal.fftconvolve(vim_hist,PSTH_window,mode='same')

    # scale by the mean FR
    mean_FR_vim = (len(vim_ev['times'])/vim_n)/tstop
    scale = mean_FR_vim/np.mean(kernel_PSTH_vim_ori)
    kernel_PSTH_vim = scale*kernel_PSTH_vim_ori

    return time_vec, kernel_PSTH_vim


syn_wt_scale = 12 # 10; 12
# Balanced
weights_BA = np.array([
    0.529,
    5.99e-3,
    -8.95,
    -1.95e-2
])*syn_wt_scale

# Hebbian
weights_H = np.array([
    0.282,
    3.01e-3,
    -0.764,
    -0.147
])*syn_wt_scale

def run_BA_and_H(exc_static, inh_static):
    time_vec, BA = runSim(
        weights_BA, label="BA",
        exc_static=exc_static,
        inh_static=inh_static
    )
    _, H = runSim(
        weights_H, label="H",
        exc_static=exc_static,
        inh_static=inh_static        
    )
    return time_vec, BA, H

###
# Plot and output
###
def plot_BA_and_H(t, BA, H):
    plt.plot(t, BA, color='k', label='BA')
    plt.plot(t, H, color='b', label='Hebbian')
    plt.legend()
    plt.title(f"PSTH ({DBS_freq} Hz) - Vim")
    plt.xlabel("Time (ms)")
    plt.ylabel("Frequency (Hz)")


#plt.figure()
#testSynapseAtFreqs(True)
#testSynapseAtFreqs(False)
    
plt.figure(figsize=(16,12))
t, BA, H = run_BA_and_H(False, False) # TM synapses in all network connections
BA = BA - 10 # consistency with experimental data
plot_BA_and_H(t, BA, H)


io.savemat(f"PSTH_{DBS_freq}Hz_BA.mat", {"data":BA})
io.savemat(f"PSTH_{DBS_freq}Hz_Hebbian.mat", {"data":H})

plt.show()



# raster plots

# plt.figure(figsize=(16,12))
# plt.subplot(221)
# rasterPlot(vim_ev)
# plt.xlabel("Time (ms)")
# plt.subplot(222)
# rasterPlot(cer_ev)
# plt.xlabel("Time (ms)")
# plt.subplot(223)
# rasterPlot(trn_ev)
# plt.xlabel("Time (ms)")
# plt.savefig(f"img/Raster_{freq}.pdf")
# plt.show()
