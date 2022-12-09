
'''
Workaround for JoSIM SFQ simulator.

A collection of functions to place common parameterized SFQ elements within a
netlist.


'''

import numpy as np
from subprocess import call

########################## GLOBAL PARAMETERS ##########################

# 2.067833848E-15/(4*2.5*0.0001)
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


def ramp(tdelay, imax, init = 0):
    'Create a piecewise-linear definition of ramp source'
    return f'pwl(0 {init} {tdelay} {imax})'

def pulse(time, pmax0 = 5e-4, pulse_width = 6e-12, scale = 1):
    '''Create a piecewise-linear sequence with triangular pulses at specific time.
    Inputs:
        time        : list of pulse peaks (in seconds)
        pmax0       : height of a single pulse. Default is 5e-4
        pulse_width : width of a pulse. Default is 6e-12
        scale       : multiplier for the height (and, consequently, area) of the
                    pulse. Can be used to feed a simultaneous clock signal to
                    multiple elements.
    WARNING: this function cannot combine closely spaced pulses. ValueError is
    raised if the peaks are placed too closely.
    '''
    differences = np.diff(time)
    if np.any(differences < pulse_width):
        raise ValueError('The pulse peaks are spaced too closely')
    pmax = pmax0 * scale
    hw = pulse_width / 2 #half width
    triangles = ['pwl(0 0 '] + [str(q) for t in time for q in (t-hw, 0, t, pmax, t+hw, 0)] + [')']
    return ' '.join(triangles)

def jj(f, n1, n2, model='jjmit', A = None, AF = 1, B0Rs = B0Rs, Rsheet = Rsheet, Lsheet = Lsheet, LP = LP, add_LP=True):
    '''Create a JJ with parasitics based on MITLL technology.
    INPUTS:
        f : target file handle
        n1, n2 : input/output nodes
        model : name of the jj model specified in .cir file
        A : area of the JJ
        AF : convenience parameter for area factor
        add_LP : flag indicating whether the LP inductance is added to the shunt
        inductance
    '''
    if AF is not None:
        A = AF * IC
    rb1 = B0Rs/A
    lrb1 = (rb1/Rsheet)*Lsheet
    if add_LP:
        lrb1 += LP
    if str(n2) == '0':
        n_mid_jj = f'{n1}_mid_shunt'
        n_mid_bias = f'{n1}_mid'
        n_jct = f'{n1}_jct'
        return f.write(f'''
************************ JJ from {n1} to GND ************************
B_{n1}_0   {n1}   {n_mid_jj} {model} area={A}
LP_{n1}_0  {n_mid_jj}   {n_jct} {LP}
RB_{n1}_0  {n1} {n_mid_bias} {rb1}
LRB_{n1}_0 {n_mid_bias} {n_jct} {lrb1}
LJCT_{n1}_0 {n_jct} {n2} 1e-18
''')
    else:
        n_mid_bias = f'{n1}_{n2}_mid'
        n_jct = f'{n1}_{n2}_jct'
        return f.write(f'''
************************ JJ from {n1} to {n2} ************************
B_{n1}_{n2}  {n1}  {n_jct} {model} area={A}
RB_{n1}_{n2}  {n1} {n_mid_bias} {rb1}
LRB_{n1}_{n2} {n_mid_bias} {n_jct} {lrb1}
LJCT_{n1}_{n2} {n_jct} {n2} 1e-18
''')

'''
* BIAS 8 ++++++++++++++++++++++++++++++++++++++++++++++++++
.param LB2=LB
LB2  8 10 LB2
.param IB2=IB1
IB2  0 10 pwl(0 0 5p IB2)
'''
def bias(f, n, IB = None, IF = None, LB=LB):
    ''' IB is ignored when IF is given '''
    if IF:
        IB = IF*Ic0*IC
    n_mid = f'{n}_mid_bias'
    return f.write(f'''
************************ Bias line for {n} ************************
L_bias_{n}  {n} {n_mid} {LB}
I_bias_{n}  0 {n_mid} {IB}
''')

def rl(n1, n2, R, L):
    mid_node = f'{n1}_{n2}_rl'
    return f'''
L_{n1}_{n2} {n1} {mid_node} {L}
R_{n1}_{n2} {mid_node} {n2} {R}
'''

def ind(f, n1, n2, L = None, LF = None):
    #LF is inductance scaled by Phi0/(4*IC*Ic0)
    if LF is None:
        return f.write(f'\nL_{n1}_{n2} {n1} {n2} {L}\n')
    elif L is None:
        return f.write(f'\nL_{n1}_{n2} {n1} {n2} {LF * Phi0/(4*IC*Ic0)}\n')
    else:
        print(f'L = {L}')
        print(f'L = {LF}')
        raise ValueError('Either L or LF should be None')

# def LF(f, n1, n2, L):
#     return f.write(f'''
# L_{n1}_{n2} {n1} {n2} {L * Phi0/(4*IC*Ic0)}
# ''')

def res(f, n1, n2, R):
    return f.write(f'''
R_{n1}_{n2} {n1} {n2} {R}
''')

def make_in_clk(n_inputs, clk, nclk, tstep = 1e-12, T = 240e-12, cycles = 1):
    ''' Utility to add clock and input pulses to the circuit
    INPUTS:
        - n_inputs - number of input waveforms to create
        - clk - clock node
        - nclk - number of clock sinks. Needed to scale the clock signal, since the current is split between the sinks
        - T - clock period
        - cycles - int > 2, number of times the sequence is run. Default is 1, i.e. one time without repetition    '''
    n_period = 2 ** n_inputs * cycles
    clk_pulse = T * np.arange(1, n_period + 1)
    _, *offsets, _ = np.linspace(0, T, n_inputs+2)
    in_pulses = [(clk_pulse - T + off) for off in offsets] # could be made with an outer product insteadcal
    tmax=  T*(n_period+1)
    if nclk:
        clk_pulse = f'Iclk 0   {clk}  ' + pulse(clk_pulse, scale = nclk) + '\n'
        clk_meas = f'.print DEVI Iclk\n'
    else:
        clk_pulse = ''
        clk_meas = ''
    return in_pulses, clk_pulse, clk_meas, n_period, tmax

def boiler_plate(f, in_nodes, clk, nclk, out_nodes, T = 240e-12, tstep = 1e-12, analysis = {}, cycles = 1):

    n_inputs = len(in_nodes)
    in_pulses, clk_pulse, clk_meas, n_period, tmax = make_in_clk(n_inputs, clk, nclk, T = T, cycles = cycles)

    for i, (n, in_pulse) in enumerate(zip(in_nodes, in_pulses)):
        idx = np.arange(n_period) % (2**(i+1)) >= 2**i
        f.write(f'I{n}   0    {n}  ' + pulse(in_pulse[idx]) + '\n')
        f.write(f'.print DEVI I{n}\n')

    for n in out_nodes :
        f.write(f'ROUT_{n} {n}  0   1\n')
        f.write(f'.print DEVI ROUT_{n}\n')

    for dev in analysis.get('dev',[]) :
        f.write(f'.print DEVI {dev}\n')

    for n in analysis.get('node',[]) :
        f.write(f'.print V {n}\n')

    for n in analysis.get('phase',[]) :
        f.write(f'.print PHASE {n}\n')

    f.write(clk_pulse)
    f.write(clk_meas)

    f.write(f'.tran {tstep} {tmax}\n')
    f.write('.model jjmit jj( rtype=1, vg=2.6mV, cap=0.07pF, r0=160, rn=16, icrit=0.1mA)\n')

def merge(f, name, a, b, q, mbias = 1,):
    f.write(f'\n*** Merger {name} connecting {a} and {b} to {q} *** \n')

    jj (f, f'{name}_a_1',          f'0', AF = 1)
    jj (f, f'{name}_a_2',          f'0', AF = 1)
    jj (f, f'{name}_a_2', f'{name}_a_3', AF = 1/1.4)

    jj (f, f'{name}_b_1',          f'0', AF = 1)
    jj (f, f'{name}_b_2',          f'0', AF = 1)
    jj (f, f'{name}_b_2', f'{name}_b_3', AF = 1/1.4)

    jj (f, f'{name}_q_2',          f'0', AF = 1)
    jj (f, f'{name}_q_1',          f'0', AF = 1)

    bias(f, f'{name}_a_1', IF = BiasCoef)
    bias(f, f'{name}_b_1', IF = BiasCoef)
    bias(f, f'{name}_q_3', IF = mbias)
    bias(f, f'{name}_q_2', IF = BiasCoef)
    bias(f, f'{name}_q_1', IF = BiasCoef)

    ind(f,    a         , f'{name}_a_1', LF = 1)
    ind(f, f'{name}_a_1', f'{name}_a_2', LF = 2)
    ind(f, f'{name}_a_3', f'{name}_q_3', L  = 1.2e-12)
    ind(f,    b         , f'{name}_b_1', LF = 1)
    ind(f, f'{name}_b_1', f'{name}_b_2', LF = 2)
    ind(f, f'{name}_b_3', f'{name}_q_3', L  = 1.2e-12)

    ind(f, f'{name}_q_3', f'{name}_q_2', LF = 2)
    ind(f, f'{name}_q_2', f'{name}_q_1', LF = 2)
    ind(f, f'{name}_q_1',          q   , LF = 1)

def dff(f, name, a, t, q,):
    f.write(f'\n*** DFF {name} connecting {a} to {q} *** \n')
    jj (f, f'{name}_1',        f'0', AF = 1)
    jj (f, f'{name}_2', f'{name}_3', AF = 1/1.4)
    jj (f, f'{name}_3',        f'0', AF = 1)
    jj (f, f'{name}_4',        f'0', AF = 1)
    jj (f, f'{name}_5', f'{name}_4', AF = 1/1.4)
    jj (f, f'{name}_t',        f'0', AF = 1)
    jj (f, f'{name}_6',        f'0', AF = 1)

    bias(f, f'{name}_1', IF = BiasCoef)
    bias(f, f'{name}_3', IF = 1) # .4
    bias(f, f'{name}_t', IF = BiasCoef)
    bias(f, f'{name}_6', IF = BiasCoef)

    ind(f,        a   , f'{name}_1', LF = 1)
    ind(f, f'{name}_1', f'{name}_2', LF = 2)
    ind(f, f'{name}_3', f'{name}_4', LF = 4)
    ind(f, f'{name}_5', f'{name}_t', LF = 2)
    ind(f,        t   , f'{name}_t', LF = 1)
    ind(f, f'{name}_4', f'{name}_6', LF = 2)
    ind(f, f'{name}_6',        q   , LF = 1)

def buffer_branch(f, name, a, q, bias_pts = None, sizes = {}):
    f.write(f'* BUFFER BRANCH FROM {a} TO {q}\n')
    if isinstance(bias_pts, int):
        bias_pts = [bias_pts]
    elif bias_pts is None:
        bias_pts = []
    for pt in bias_pts:
        if pt == 0:
            bias(f, a, IF = BiasCoef)
        elif pt == 3:
            bias(f, q, IF = BiasCoef)
        elif pt in (1,2):
            bias(f, f'{name}_{pt}', IF = BiasCoef)
        else:
            raise ValueError(f'Invalid bias point {bias_pts}')

    ind(f,           a , f'{name}_1', LF = sizes.get('ind1',1    ))
    jj (f,  f'{name}_1',        f'0', AF = sizes.get('jj1', 1    ))
    ind(f,  f'{name}_1', f'{name}_2', LF = sizes.get('ind2',2    ))
    jj (f,  f'{name}_2',          q , AF = sizes.get('ind2',1/1.4))

def comparator(f, a, t, q, size_top = 0.8, size_bot = 1):
    f.write(f'* COMPARATOR BETWEEN {t} AND {a}\n')
    jj( f, t,   a, AF = size_top)
    jj( f, a, '0', AF = size_bot)
    ind(f, a,   q, LF = 2)


def jtl(f, name, a, q, LF = 2, AF = 1, in_ind = True, out_ind = True, njj = 1):
    # TODO: add bias lines
    if in_ind:
        if out_ind:
            nodes = [a] + [f'{name}_{i}' for i in range(1, njj+1)] + [q]
            jjs =         [f'{name}_{i}' for i in range(1, njj+1)]
        else:
            nodes = [a] + [f'{name}_{i}' for i in range(1, njj  )] + [q]
            jjs =         [f'{name}_{i}' for i in range(1, njj  )] + [q]
    else:
        if out_ind:
            nodes = [a] + [f'{name}_{i}' for i in range(1, njj  )] + [q]
            jjs =   [a] + [f'{name}_{i}' for i in range(1, njj  )]
        else:
            nodes = [a] + [f'{name}_{i}' for i in range(1, njj-1)] + [q]
            jjs =   [a] + [f'{name}_{i}' for i in range(1, njj-1)] + [q]

    for i,j in zip(nodes[:-1], nodes[1:]):
        ind(f, i, j, LF = LF)
    for i in jjs:
        jj( f, i, '0', AF = 1)

def rdff(f, name, a, r, t, q):
    f.write(f'***//*** RDFF {a}, {r} -> {q} ***//***\n')

    # branch A
    ind(f,            a, f'{name}_a1', LF = 1)
    bias(f,f'{name}_a1', IF = BiasCoef)
    jj( f, f'{name}_a1',          '0', AF = 1)
    ind(f, f'{name}_a1', f'{name}_ar', LF = 1)

    # branch R
    ind(f,            r, f'{name}_r1', LF = 1)
    bias(f,f'{name}_r1', IF = BiasCoef)
    jj( f, f'{name}_r1', f'{name}_r2', AF = 1/1.4)
    jj( f, f'{name}_r1',          '0', AF = 1)
    ind(f, f'{name}_r2', f'{name}_r3', LF = 1)

    jj( f, f'{name}_r3', f'{name}_r4', AF = 1/1.4)
    jj( f, f'{name}_r3',          '0', AF = 1)
    ind(f, f'{name}_r4', f'{name}_ar', LF = 1)

    # branch T
    ind(f,            t, f'{name}_t1', LF = 1)
    bias(f,f'{name}_t1', IF = BiasCoef)
    jj( f, f'{name}_t1', f'{name}_t2', AF = 1.4)
    jj( f, f'{name}_t1',          '0', AF = 1)
    ind(f, f'{name}_t2', f'{name}_to', LF = 1)

    # Bridge
    jj( f, f'{name}_to', f'{name}_b1', AF = 1)
    jj( f, f'{name}_to',          '0', AF = 1/1.4)
    ind(f, f'{name}_b1', f'{name}_ar', LF = 1)

    # branch Out
    ind(f,           q , f'{name}_q1', LF = 1)
    jj( f, f'{name}_q1',          '0', AF = 1)
    ind(f, f'{name}_q1', f'{name}_q2', LF = 1)
    bias(f,f'{name}_q2', IF = BiasCoef)
    ind(f, f'{name}_q2', f'{name}_q3', LF = 1)
    jj( f, f'{name}_q3',          '0', AF = 1)
    ind(f, f'{name}_q3', f'{name}_to', LF = 1)


    # jj( f, t,   a, AF = size_top)
    # jj( f, a, '0', AF = size_bot)
    # ind(f, a,   q, LF = 2)
