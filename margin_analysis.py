import csv
from matplotlib import pyplot as plt
import numpy as np
import os
from collections import deque
from itertools import product, islice

from subprocess import call, DEVNULL
from requests.structures import CaseInsensitiveDict

from sim_flatten_new import convert_sim
import multiprocessing



# from pdb import set_trace

'''
* JJ 1 2 ================================================================
.param B1=IC
    B1  1  2 jjmit area=B1
.param LP1=LP
    LP1  2 0 LP1
.param RB1=B0Rs/B1
    RB1 1 101 RB1
.param LRB1=(RB1/Rsheet)*Lsheet
    LRB1 101 0 LRB1

* JJ_NG 4 6 ================================================================
.param B3=IC/1.4
B3  4  6 jjmit area=B3
.param RB3=B0Rs/B3
RB3 4 106 RB3
.param LRB3=(RB3/Rsheet)*Lsheet
LRB3 106 6 LRB3
'''

########################## GLOBAL PARAMETERS ######################################################
Phi0 = 2.067833848E-15
B0 = 1
Ic0 = 0.0001
IcRs = 100e-6*6.859904418
B0Rs = IcRs/Ic0*B0
Rsheet = 2
Lsheet = 1.13e-12
LP = 0.2e-12
IC = 2.5
LB = 2e-12
BiasCoef=0.70
IB1=BiasCoef*Ic0*IC

################################################################################################################

def ramp(tdelay, imax, init = 0):
    return f'pwl(0 0 {tdelay} {imax})'
def pulse(time, pmin = 0, pmax = 5e-4, factor = 1):
    pw = 3e-15 / pmax
    # half_pw = pw/2
    # breakpoint()
    triangles = ['pwl(0 0 '] + [str(q) for t in time for q in (t-pw/2, pmin, t, pmax, t+pw/2, 0)] + [')']
    return ' '.join(triangles)

def jj(f, n1, n2, model='jjmit', area = 1, B0Rs = B0Rs, Rsheet = Rsheet, Lsheet = Lsheet, LP = LP, add_LP=True):
    rb1 = B0Rs/area
    lrb1 = (rb1/Rsheet)*Lsheet
    if add_LP:
        lrb1 += LP
    if str(n2) == '0':
        n_mid_jj = f'{n1}_mid_shunt'
        n_mid_bias = f'{n1}_mid'
        n_jct = f'{n1}_jct'
        return f.write(f'''\nB_{n1}_0   {n1}   {n_mid_jj} jjmit area={area}\nLP_{n1}_0  {n_mid_jj}   0 {LP}\nRB_{n1}_0  {n1} {n_mid_bias} {rb1}\nLRB_{n1}_0 {n_mid_bias} 0 {lrb1}\n''')
    else:
        n_mid_bias = f'{n1}_{n2}_mid'
        n_jct = f'{n1}_{n2}_jct'
        return f.write(f'''\nB_{n1}_{n2}  {n1}  {n2} jjmit area={area}\nRB_{n1}_{n2}  {n1} {n_mid_bias} {rb1}\nLRB_{n1}_{n2} {n_mid_bias} {n2} {lrb1}\n''')

def bias(f, n, IB, LB):
    n_mid = f'{n}_mid_bias'
    return f.write(f'''\nL_bias_{n}  {n} {n_mid} {LB}\nI_bias_{n}  0 {n_mid} {IB}\n''')
def rl(n1, n2, R, L):
    mid_node = f'{n1}_{n2}_rl'
    return f'''\nL_{n1}_{n2} {n1} {mid_node} {L}\nR_{n1}_{n2} {mid_node} {n2} {R}\n'''
def ind(f, n1, n2, L):
    return f.write(f'''\nL_{n1}_{n2} {n1} {n2} {L}\n''')
def res(n1, n2, R):
    return f.write(f'''\nR_{n1}_{n2} {n1} {n2} {R}\n''')

def jtl_diode(n_in, n_out, IC = IC, add_LP = False, LP = LP, BiasCoef = BiasCoef, Ic0 = Ic0, Lout = 1.2e-12):
    bias(f,       n_in, BiasCoef*Ic0*IC, LB)
    ind(f,        n_in,  f'{n_in}1', Phi0/(4*IC*Ic0))
    jj (f,  f'{n_in}1',         '0', area = IC, add_LP = add_LP)
    ind(f,  f'{n_in}1',  f'{n_in}2', Phi0/(2*IC*Ic0))
    jj (f,  f'{n_in}2',         '0', area = IC, add_LP = add_LP)
    ind(f,  f'{n_in}2',  f'{n_in}3', 1.2e-12        )
    # escape jj
    jj (f,  f'{n_in}3',       n_out, area = IC, add_LP = add_LP)
    bias(f,      n_out, BiasCoef*Ic0*IC, LB)

def jtl(f, a, q):

    ind(f, f'{a}0' , f'{a}1', Phi0/(4*IC*Ic0))
    jj( f, f'{a}1' ,    f'0', area = IC, add_LP = True)
    ind(f, f'{a}1' , f'{a}2', Phi0/(4*IC*Ic0))

    bias(f, f'{a}2' , 2*BiasCoef*Ic0*IC, LB)

    ind(f, f'{a}2' , f'{a}3', Phi0/(4*IC*Ic0))
    jj( f, f'{a}3' ,    f'0', area = IC, add_LP = True)
    ind(f, f'{a}3' ,       q, Phi0/(4*IC*Ic0))

def dff_buf(f, a, clk, q, factor = (1,1,1,1)):

    bias(f, f'{a}1' , BiasCoef*Ic0*IC * factor[0], LB)
    bias(f, f'{a}3' ,          Ic0*IC * factor[1], LB)
    bias(f, f'{a}t1', BiasCoef*Ic0*IC * factor[2], LB)
    bias(f, f'{a}5' , BiasCoef*Ic0*IC * factor[3], LB)

    # Input
    ind(f, f'{a}0' , f'{a}1', Phi0/(4*IC*Ic0))
    jj( f, f'{a}1' ,    f'0', area = IC    )
    ind(f, f'{a}1' , f'{a}2', Phi0/(2*IC*Ic0))
    jj( f, f'{a}2' , f'{a}3', area = IC/1.4)
    jj( f, f'{a}3' ,    f'0', area = IC    )
    ind(f, f'{a}3' , f'{a}4', Phi0/(  IC*Ic0))
    jj( f, f'{a}4' ,    f'0', area = IC    )
    #Clock
    ind(f,    clk  ,f'{a}t1', Phi0/(4*IC*Ic0))
    jj( f, f'{a}t1',    f'0', area = IC    )
    ind(f, f'{a}t1',f'{a}t2', Phi0/(2*IC*Ic0))
    jj( f, f'{a}t2',f'{a}4' , area = IC/1.4)
    #Output
    ind(f, f'{a}4' , f'{a}5', Phi0/(2*IC*Ic0))
    jj( f, f'{a}5' ,    f'0', area = IC    )
    jj( f, f'{a}5' , f'{a}6', area = IC    )
    ind(f, f'{a}6' ,       q, Phi0/(4*IC*Ic0))


def dff(f, a, clk, q):
    bias(f, f'{a}1' , BiasCoef*Ic0*IC, LB)
    bias(f, f'{a}3' ,          Ic0*IC, LB)
    bias(f, f'{a}t1', BiasCoef*Ic0*IC, LB)
    bias(f, f'{a}5' , BiasCoef*Ic0*IC, LB)

    # Input
    ind(f,    a    , f'{a}1', Phi0/(4*IC*Ic0))
    jj( f, f'{a}1' ,    f'0', area = IC    )
    ind(f, f'{a}1' , f'{a}2', Phi0/(2*IC*Ic0))
    jj( f, f'{a}2' , f'{a}3', area = IC/1.4)
    jj( f, f'{a}3' ,    f'0', area = IC    )
    ind(f, f'{a}3' , f'{a}4', Phi0/(  IC*Ic0))
    jj( f, f'{a}4' ,    f'0', area = IC    )
    #Clock
    ind(f,    clk  ,f'{a}t1', Phi0/(4*IC*Ic0))
    jj( f, f'{a}t1',    f'0', area = IC    )
    ind(f, f'{a}t1',f'{a}t2', Phi0/(2*IC*Ic0))
    jj( f, f'{a}t2',f'{a}4' , area = IC/1.4)
    #Output
    ind(f, f'{a}4' , f'{a}5', Phi0/(2*IC*Ic0))
    jj( f, f'{a}5' ,    f'0', area = IC    )
    ind(f, f'{a}5' ,       q, Phi0/(4*IC*Ic0))

def xor_branch(f, n_in, jct, factor, Phi0 = Phi0, IC = IC, BiasCoef = BiasCoef, Ic0 = Ic0, LB = LB, add_LP = False, LP = LP):
    bias(f, f'{n_in}1', BiasCoef*Ic0*IC * factor, LB)
    ind(f,     n_in   , f'{n_in}1', Phi0/(4*IC*Ic0))
    jj( f,  f'{n_in}1',        '0', area = IC, add_LP = True)
    ind(f,  f'{n_in}1', f'{n_in}2', Phi0/(2*IC*Ic0))
    jj( f,  f'{n_in}2',        '0', area = IC, add_LP = True)
    ind(f,  f'{n_in}2', f'{n_in}3', 1.2e-12        )
    jj( f,  f'{n_in}3', f'{n_in}4', area = IC * factor, add_LP = True)
    bias(f, f'{n_in}4', BiasCoef*Ic0*IC * factor, LB)
    ind(f,  f'{n_in}4',        jct, Phi0/(2*IC*Ic0))


# def readout(threshold, tvec, ivec, tstep, clk_pulse, twidth = 40e-12): # tstep = tstep,
#     out = []
#     for tstart in clk_pulse:
#         idx = np.logical_and(tvec > tstart, tvec <= tstart + twidth)
#         s = np.trapz(ivec[idx], dx=tstep)
#         out.append(s > threshold)
#     return np.array(out)

# def double_pulse():
def detect_pulse(time, Ivector, clk_pulse, tstep, pulse_area, crossing_th = 250e-6, hold_time = 40e-12, ncycles = None):
# def detect_pulse(tvec, Ivector, clk_pulse, tstep, pulse_area, crossing_th = 300e-6, hold_time = 40e-12):
    ''' Check whether a pulse is produced  '''
    out = []
    t_idx_cycle_start = time.searchsorted(clk_pulse)
    t_idx_hold_end    = time.searchsorted(clk_pulse + hold_time)
    t_idx_cycle_end   = np.append(t_idx_cycle_start[1:], len(time))
        
    ncycles = len(clk_pulse) if ncycles is None else ncycles
    for idx_cycle_start, idx_hold_end, idx_cycle_end in islice(zip(t_idx_cycle_start, t_idx_hold_end, t_idx_cycle_end), ncycles):
        I_during_hold = Ivector[idx_cycle_start:idx_hold_end]
        I_after_hold = Ivector[idx_hold_end:idx_cycle_end]
        if np.any(I_after_hold > crossing_th):
            # out.append('pulse outside initial window')
            out.append('E')
        elif np.any(I_during_hold > crossing_th):
            s = np.trapz(I_during_hold, dx=tstep)
            out.append('1' if s > pulse_area else '0')
        else:
            out.append('0')
    return ''.join(out)
        # elif not np.any(np.diff(Ivec) > threshold):
        #     return False

def readout(threshold, tvec, ivec, tstep, clk_pulse, twidth = 40e-12): # tstep = tstep,
    out = []
    for tstart in clk_pulse:
        idx = np.logical_and(tvec > tstart, tvec <= tstart + twidth)
        s = np.trapz(ivec[idx], dx=tstep)
        out.append(s > threshold)
    return np.array(out)

def is_reset(time, Ivector, clk_pulse, hold_time, window_width, target, tol):
    left = clk_pulse + hold_time
    t_idx_left = time.searchsorted(left)
    t_idx_right = time.searchsorted(left + window_width)
    for idx_l, idx_r in zip(t_idx_left, t_idx_right):
        Ivec = Ivector[idx_l:idx_r+1]
        Imin, Imax = Ivec.min(), Ivec.max()
        if (Imax - target) > tol or (target - Imin) > tol:
            return False
    return True

def is_reset_phase(time, phase, clk_pulse, hold_time, window_width, tol = 0.01, target = None):
    phase = np.remainder(phase, 2*np.pi)
    target = phase[0] if target is None else target
    left = clk_pulse + hold_time
    t_idx_left  = time.searchsorted(left)
    t_idx_right = time.searchsorted(left + window_width)
    # breakpoint()
    for i, (idx_l, idx_r) in enumerate(zip(t_idx_left, t_idx_right)):
        phase_range = phase[idx_l:idx_r+1]
        Imin, Imax = phase_range.min(), phase_range.max()
        # print(Imin, Imax, target, tol, (Imax - target), (target - Imin), (Imax - target) > tol, (target - Imin) > tol)
        if (Imax - target) > tol or (target - Imin) > tol:
            # print(f'Output JJ not reset in cycle {i}')
            return False
    return True

def _branch(f, a, factor = [1,1], r = False):
    # bias(f,f'{a}1', ramp(5e-12, BiasCoef*Ic0*IC), LB)
    # bias(f,f'{a}3', ramp(5e-12, BiasCoef*Ic0*IC), LB)
    bias(f,f'{a}1', BiasCoef*Ic0*IC * factor[0], LB)
    bias(f,f'{a}3', BiasCoef*Ic0*IC * factor[0], LB)

    ind(f, f'{a}0', f'{a}1', Phi0/(4*IC*Ic0) )# * x[ 0])#4
    jj( f, f'{a}1',    f'0', area = IC       )# * x[ 1])
    ind(f, f'{a}1', f'{a}2', Phi0/(2*IC*Ic0) )# * x[ 2])#2
    jj( f, f'{a}2', f'{a}3', area = IC/1.4   )# * x[ 3])
    jj( f, f'{a}3',    f'0', area = IC       )
    ind(f, f'{a}3', f'{a}4', Phi0/(1*IC*Ic0) )#1
    jj( f, f'{a}4',    f'0', area = IC       )
    ind(f, f'{a}4', f'{a}5', Phi0/(2*IC*Ic0) )#2
    if not r:
        jj( f, f'{a}5',   'jct', area = IC/1.4   )
    else:
        ind(f, f'{a}3', f'{a}4e', Phi0/(1*IC*Ic0) )#1
        jj( f, f'{a}4e',    f'0', area = IC       )
        ind(f, f'{a}4e', f'{a}5e', Phi0/(2*IC*Ic0) )#2
        jj( f, f'{a}5e',   'jct', area = IC/1.4   )


    ind(f, f'{a}t',   f't3', 1e-12         * factor[1])
    jj( f, f'{a}4', f'{a}t', area = IC/1.4 * factor[1])

def gen_ckt(factor, file, tstep, tmax):
    with open(file, 'w') as f:
        f.write('*\n\n')

        _branch(f, 'a') #, factor, r = True
        _branch(f, 'b')
        _branch(f, 'c')
        _branch(f, 'd')

        # clock delivery
        ind(f,'clk', 't1', Phi0/(4*IC*Ic0))
        jj( f, 't1',  '0', area = IC      )
        ind(f, 't1', 't2', Phi0/(2*IC*Ic0))
        jj( f, 't2',  '0', area = IC      )
        ind(f, 't2', 't3', 1e-12          )
        bias(f,'t1', BiasCoef*Ic0*IC * factor[0], LB) # has to be between 1.2 and 1.4
        bias(f,'t2', BiasCoef*Ic0*IC * factor[1], LB) # has to be between 1.2 and 1.4
        bias(f,'q1', BiasCoef*Ic0*IC * factor[2], LB) # has to be between 1.2 and 1.4
        # bias(f,'t1', BiasCoef*Ic0*IC * factor[2], LB) # has to be between 1.2 and 1.4
        # bias(f,'t2', BiasCoef*Ic0*IC * factor[2], LB) # has to be between 1.2 and 1.4
        # bias(f,'q1', BiasCoef*Ic0*IC * factor[2], LB) # has to be between 1.2 and 1.4

        # output branch
        ind(f,'jct', 'q1', 1e-12            )
        jj( f, 'q1',  '0', area = IC        )
        ind(f, 'q1', 'q0', Phi0/(4*IC*Ic0)  )

        f.write(f'.tran {tstep} {tmax} 0 {tstep}\n')

        f.write('Ia   0    a0  ' + pulse([p for p,v in zip(in_pulse[0], [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]) if v]) + '\n')
        f.write('Ib   0    b0  ' + pulse([p for p,v in zip(in_pulse[1], [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1]) if v]) + '\n')
        f.write('Ic   0    c0  ' + pulse([p for p,v in zip(in_pulse[2], [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1]) if v]) + '\n')
        f.write('Id   0    d0  ' + pulse([p for p,v in zip(in_pulse[3], [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1]) if v]) + '\n')
        f.write('Iclk 0   clk  ' + pulse(clk_pulse) + '\n')
        f.write('ROUT q0  0   1\n')

        f.write('.model jjmit jj( rtype=1, vg=2.6mV, cap=0.07pF, r0=160, rn=16, icrit=0.1mA)\n')

        f.write(f'.print DEVI ROUT\n')
        f.write(f'.print DEVI LP_q1_0\n')
        f.write(f'.print DEVI LRB_q1_0\n')
        f.write(f'.print PHASE B_q1_0\n')
        f.write('.end')

def axvlines(ax, xvals, color = 'r', linestyle = '--', lw = 0.5):
    for x in xvals:
        ax.axvline(x, color = color, linestyle = linestyle, lw = lw)

def sim_ckt(file, csv_file = None):
    # call(['josim-cli', '-o', csv_file, file, '-m', '1'], stdout=DEVNULL)
    # call(['josim-cli', '-o', csv_file, file, '-m', '1'], stdout=DEVNULL)
    csv_file = convert_sim(file, csv_file)
    with open(csv_file) as f:
        # reading the CSV file
        csvFile = csv.reader(f)
        # displaying the contents of the CSV file
        varnames = next(csvFile)
        vars = []
        for lines in csvFile:
            vars.append([float(q) for q in lines])
    return varnames, np.array(vars)

def sim_ckt_dict(file, csv_file = None, scale_time = 1):
    if csv_file is None:
        base, ext = os.path.splitext(file)
        csv_file = f'{base}.csv'
    call(['josim-cli', '-o', csv_file, file, '-m', '1'], stdout=DEVNULL)
    with open(csv_file) as f:
        # reading the CSV file
        csvFile = csv.reader(f)
        # displaying the contents of the CSV file
        varnames = next(csvFile)
        vars = np.array([[float(q) for q in line] for line in csvFile])
    var_dict = CaseInsensitiveDict()
    for i,vname in enumerate(varnames):
        if vname == 'time':
            var_dict['time'] = vars[:,i] * scale_time
        else:
            # _type = vname[0]
            element = vname[2:-1] #could use regex but seems stable so far
            var_dict[vname] = vars[:,i]
            # var_dict.setdefault(_type,CaseInsensitiveDict())
            # var_dict[_type][element] = vars[:,i]
    # print(csv_file)
    # print(list(var_dict.keys()))
    return var_dict

def from_template(params: dict, template: str):
    '''
    params is a dictionary describing the minimum and maximum factors that can be taken
    '''
    vars = []
    variants = []
    for var, values in params.items():
        vars.append(var)
        variants.append(values)
    
    with open(template) as f:
        lines = f.readlines()
    
    base, ext = os.path.splitext(template)
    names = {}
    for i, values in enumerate(product(*variants)):
        suffix = '_'.join(str(v) for v in values)
        filename = f'{base}_{suffix}{ext}'
        names[tuple(values)] = filename
        with open(filename, 'w') as f: 
            f.write('*\n')
            for var, v in zip(vars, values):
                f.write(f'.param {var}={v:g}\n')
            for line in lines:
                f.write(line)
    return names

if __name__ == '__main__':

    # TODO: replace current sources with generated ones
    # TODO: vary bias parameters 
    # TODO: Ask Tahereh how to do it right
    
    T = 250e-12
    n_inputs = 4
    nperiod = 2 ** n_inputs * 2
    clk_pulse = T * np.arange(1, nperiod + 1)# [240e-12,480e-12,720e-12,960e-12]
    _, *offsets, _ = np.linspace(0, T, n_inputs+2)
    in_pulse = [(clk_pulse - T + off) for off in offsets]
    # tstep = 0.25e-12
    tstep = 1e-12
    tmax  = T*(nperiod+1)
    pulse_area  = 1.25e-15
    hold_time = 45e-12
    window_width = 5e-12
    crossing_th = 250e-6
    itol = 10e-6
    N = len([0])

    params = {  'bf_xa'     :(0.8, 1, 1.2),
                'bf_xb'     :(0.8, 1, 1.2),
                'bf_xab'    :(0.8, 1, 1.2),
                'bf_xc'     :(0.8, 1, 1.2),
                'bf_xd'     :(0.8, 1, 1.2),
                'bf_xcd'    :(0.8, 1, 1.2),
                'bf_xabcd'  :(0.8, 1, 1.2),}
    template = 'or4_template.cir'
    
    names = from_template(params, template)
    pool = multiprocessing.Pool(40)
    results = zip(*pool.map(sim_ckt, list(names.values())))
    
    # result = {}
    # for key, file in names.items():
    #     file = 'or4.cir'
    #     varnames, out = sim_ckt(file)
    #     time = out[:,0]
    #     iout = out[:,varnames.index('I(ROUT_ABCD)')]
    #     bit_seq = detect_pulse(time, iout, clk_pulse, tstep, pulse_area, crossing_th, hold_time, ncycles = 16)
    #     result[key] = bit_seq
    #     print(key, bit_seq)
    
    
    '''
    factor1 = np.arange(0.7, 1.8, 0.3)
    factor2 = np.arange(0.7, 1.8, 0.3)
    factor3 = np.arange(0.7, 1.8, 0.3)
    factor4 = np.arange(0.7, 1.8, 0.3)
    f = [q.ravel() for q in np.meshgrid(factor1, factor2, factor3, indexing = 'ij')]#factor4,
    for i, factor in enumerate(zip(*f)):
        # factor = np.array([factor])
        print(factor)
        gen_ckt(factor, file, tstep, tmax)
        vars = sim_ckt(file)

        time     = vars[:,0]
        iout     = vars[:,1]
        ijct     = vars[:,2] + vars[:,3]
        phase_JJ = vars[:,4]
        cycles = phase_JJ/2/np.pi
        # cycles -= cycles[0]
        target = ijct[0]
        # bit_seq = readout(pulse_area, time, vars[:,1], tstep, clk_pulse)
        bit_seq = detect_pulse(time, iout, clk_pulse, tstep, pulse_area, crossing_th, hold_time)
                       # is_reset_phase(time, phase   , clk_pulse, hold_time, window_width, tol = 0.01, target = None)
        reset_status = is_reset_phase(time, phase_JJ, clk_pulse, hold_time, window_width)
        # if not reset_status:
        #     continue
    '''
    '''
        fig, axs = plt.subplots(2,1, num = ' '.join(f'{q:g}' for q in factor))
        # fig.canvas.set_window_title(' '.join(factor))

        plt.tight_layout()
        axs[0].set_ylabel(' '.join(f'{q:g}' for q in factor))

        axs[0].plot(time, iout, lw=1) # time, ijct,
        axs[0].vlines( clk_pulse, ymin = 0, ymax = 0.5e-3, colors = 'r', linestyles='dashed', lw=0.5)
        for i,t in enumerate(clk_pulse):
            axs[0].text(t, 5e-4, f'{i:04b}',)
        axs[0].grid()
        # p = axs[1].plot(time, cycles, lw=1)
        # axs[1].vlines( clk_pulse, ymin = 0, ymax = cycles.max(), colors = 'r', linestyles='dashed', lw=0.5)
        axs[1].plot(time, np.remainder(phase_JJ,2*np.pi), lw=1)
        axs[1].vlines( clk_pulse, ymin = 0, ymax = 2*np.pi, colors = 'r', linestyles='dashed', lw=0.5)

        axs[1].grid()

        # axs[2].plot(time, np.remainder(phase_JJ,2*np.pi), lw=1)
        # axs[2].vlines( clk_pulse, ymin = 0, ymax = 2*np.pi, colors = 'r', linestyles='dashed', lw=0.5)
        # axs[2].grid()
        plt.show()

        print(bit_seq)
        print(reset_status)
        # break
    '''


    # x  = np.ones(3) * 1.1
    # lb = x * 0.7
    # ub = x * 1.3
    # corners(x, lb, ub, 0.25e-12)
