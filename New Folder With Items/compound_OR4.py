import csv
from matplotlib import pyplot as plt
import numpy as np
import os, sys
from collections import deque

from subprocess import call
from itertools import product

from elements import jj, ind, res, bias, pulse, boiler_plate, BiasCoef, dff, merge

################################################################################################################

registry = {}

def branch(f, name, a, q, registry = registry):
    f.write(f'\n*** Branch {name} from {a} to {q} *** \n')

    bias(f,f'{name}_1', IF = BiasCoef)
    bias(f,f'{name}_3', IF = BiasCoef)

    ind(f,     a      , f'{name}_1', LF = 1)#4
    jj (f, f'{name}_1',        f'0', AF = 1)
    ind(f, f'{name}_1', f'{name}_2', LF = 2)#2
    jj (f, f'{name}_2', f'{name}_3', AF = 1/1.4)
    jj (f, f'{name}_3',        f'0', AF = 1)
    ind(f, f'{name}_3', f'{name}_4', LF = 4)#1
    jj (f, f'{name}_4',        f'0', AF = 1)
    ind(f, f'{name}_4', f'{name}_5', LF = 2)#2
    jj (f, f'{name}_5',          q , AF = 1/1.4)

    ind(f, f'{a}t',   f't3', L  = 1e-12)
    jj (f, f'{a}4', f'{a}t', AF = 1/1.4)



def or4(f, name, a, b, c, d, t, q, registry = registry):

    dff(f, 'dff_A', a, t, f'{a}_jct', registry = registry)
    dff(f, 'dff_B', b, t, f'{b}_jct', registry = registry)
    dff(f, 'dff_C', c, t, f'{c}_jct', registry = registry)
    dff(f, 'dff_D', d, t, f'{d}_jct', registry = registry)

    merge(f, 'A_or_B', f'{a}_jct', f'{b}_jct', 'ab', registry = registry)
    merge(f, 'C_or_D', f'{c}_jct', f'{d}_jct', 'cd', registry = registry)
    merge(f, 'A_or_B_or_C_or_D', 'ab', 'cd', q, registry = registry)

def inv(f, name, a, t, q, registry = registry):
    jj (f, f'{name}_1',        f'0', AF = 1)
    jj (f, f'{name}_2',        f'0', AF = 1)
    jj (f, f'{name}_3', f'{name}_4', AF = 1/1.4)
    jj (f, f'{name}_5',        f'0', AF = 1)
    jj (f, f'{name}_6', f'{name}_7', AF = 1/1.4)
    jj (f, f'{name}_8',        f'0', AF = 1)
    jj (f, f'{name}_9', f'{name}_7', AF = 1/1.4)
    jj (f, f'{name}_7',        f'0', AF = 1)
    jj (f, f'{name}_10',       f'0', AF = 1)

    bias(f, f'{name}_1', IF = BiasCoef)
    bias(f, f'{name}_2', IF = 0.5)
    bias(f, f'{name}_4', IF = BiasCoef)
    bias(f, f'{name}_5', IF = BiasCoef)
    bias(f, f'{name}_10',IF = BiasCoef)

    ind(f,        a   , f'{name}_1' , LF = 1)
    ind(f, f'{name}_1', f'{name}_2' , LF = 2)
    ind(f, f'{name}_2', f'{name}_3' , LF = 2)
    ind(f,        t   , f'{name}_5' , LF = 1)
    ind(f, f'{name}_5', f'{name}_11', L  = 1e-12)
    ind(f, f'{name}_11',f'{name}_6' , LF = 2)
    ind(f, f'{name}_11',f'{name}_12', L  = 2e-12)
    ind(f, f'{name}_8' ,f'{name}_12', L  = 1e-12)
    ind(f, f'{name}_8' ,f'{name}_4' , LF = 4)
    ind(f, f'{name}_4' ,f'{name}_9' , L  = 1e-12)
    ind(f, f'{name}_7' ,f'{name}_10', LF = 2)
    ind(f, f'{name}_10',       q    , LF = 1)

    ind(f, f'{name}_12',f'{name}_112', L = 2e-12)
    res(f, f'{name}_112', 0, 4)

if __name__ == '__main__':
    file = 'or4_test.cir'
    # file = 'dff_test.cir'
    # file = 'or2_test.cir'
    base, ext = os.path.splitext(file)
    csv_file = f'{base}.csv'
    T = 240e-12
    in_nodes = list('abcd')
    out = 'q'
    clk = 't'


    analysis = {'dev' :[f'LJCT_dff_A_1_0',f'LJCT_dff_B_1_0',f'LJCT_dff_C_1_0',f'LJCT_dff_D_1_0',
                        f'LJCT_dff_A_3_0',f'LJCT_dff_B_3_0',f'LJCT_dff_C_3_0',f'LJCT_dff_D_3_0',
                        f'LJCT_dff_A_4_0',f'LJCT_dff_B_4_0',f'LJCT_dff_C_4_0',f'LJCT_dff_D_4_0',
                        f'LJCT_dff_A_6_0',f'LJCT_dff_B_6_0',f'LJCT_dff_C_6_0',f'LJCT_dff_D_6_0'],
                'node':[], #['ab', 'cd',],
                }

    with open(file,'w') as f:
        f.write('*\n\n')
        # or4(f, 'qqq', *in_nodes, clk, out)
        # boiler_plate(f, in_nodes, clk, [out],  T = 240e-12, analysis = analysis)
        # dff(f, 'qqq', *in_nodes, clk, out)

        dff(f, 'dff_A', 'a', 't', 'a_jct')
        dff(f, 'dff_B', 'b', 't', 'b_jct')
        merge(f, 'A_or_B', 'a_jct', 'b_jct', 'ab', mbias = 0.5)

        dff(f, 'dff_C', 'c', 't', 'c_jct')
        dff(f, 'dff_D', 'd', 't', 'd_jct')
        merge(f, 'C_or_D', 'c_jct', 'd_jct', 'cd', mbias = 0.5)

        # dff(f, 'dff_C', 'ab', 't', 'ab2')
        # dff(f, 'dff_C', 'c', 't', 'c_jct')
        merge(f, 'A_or_B_or_C_or_D', 'ab', 'cd', 'abcd', mbias = 1)
        ts = 4
        # analysis =
        boiler_plate(f, ['a','b','c','d'], 't', ts, ['abcd'],  T = 250e-12,
            analysis = {'dev':['LJCT_A_or_B_or_C_or_D_q_1_0'],
                        'phase': ['B_A_or_B_or_C_or_D_q_2_0',
                                  'B_A_or_B_or_C_or_D_q_1_0']},
            cycles = 2)
        f.write('.end')

    print(f'Simulating {file}')
    call(['josim-cli', '-o', csv_file, file, '-V', '1'])

# time = []
# opening the CSV file
# josim-cli -o ./merge_test.csv ./merge_test.cir -V 1

    print(f'Reading {csv_file}')
    with open(csv_file) as f:
        # reading the CSV file
        csvFile = csv.reader(f)
        # displaying the contents of the CSV file
        _, *varnames = next(csvFile)
        # vars = {v:[] for v in varnames}
        vars = []
        for lines in csvFile:
            vars.append([float(q) for q in lines])
        vars = np.array(vars)

    sys.exit()

















T = 240e-12
n_inputs = 4
nperiod = 2 ** n_inputs
clk_pulse = T * np.arange(1, nperiod + 1)# [240e-12,480e-12,720e-12,960e-12]
_, *offsets, _ = np.linspace(0, T, n_inputs+2)
in_pulse = [(clk_pulse - T + off) for off in offsets]
tmax=  T*(nperiod+1)
tstep = 1e-12

if __name__ == '__main__':
    file = 'or4_test.cir'
    base, ext = os.path.splitext(file)
    csv_file = f'{base}.csv'
    # options = np.array([1.2])
    options = list(product(np.linspace(0.7, 1.9, 4), repeat = 4))
    for factor in options:
        fig, (ax, axr) = plt.subplots(2, 1, sharex=True) # sharey=True,
        plt.tight_layout()
        # ax.set_ylabel(f'{factor:g}')
        ax.set_ylabel('_'.join(f'{q:g}' for q in factor))
        print(f'Creating {file} for {factor}')
        with open(file,'w') as f:
            f.write('*\n\n')


            branch('aa', 'jct', factor)
            branch('ab', 'jct', factor)
            branch('ba', 'jct', factor)

            # clock delivery
            ind(f,'clk', 't1', Phi0/(4*IC*Ic0))
            jj( f, 't1',  '0', area = IC    , add_LP = False)
            ind(f, 't1', 't2', Phi0/(2*IC*Ic0))
            jj( f, 't2',  '0', area = IC    , add_LP = False)
            ind(f, 't2', 't3', 1e-12          )

            bias(f,'t1', BiasCoef*Ic0*IC, LB)# *factor
            bias(f,'t2', BiasCoef*Ic0*IC * factor[1], LB)# *factor

            # output branch
            bias(f,'q1', BiasCoef*Ic0*IC * factor[2], LB)# *factor
            ind(f,'jct', 'q1', 1e-12          )
            jj( f, 'q1',  '0', area = IC * factor[3], add_LP = False)
            ind(f, 'q1', 'q0', Phi0/(4*IC*Ic0))

            f.write(f'.tran {tstep} {tmax} 0 0.25p\n')

            f.write('Ia   0    a0  ' + pulse([p for p,v in zip(in_pulse[0], [0,1,0,1,0,1,0,1]) if v]) + '\n')
            f.write('Ib   0    b0  ' + pulse([p for p,v in zip(in_pulse[1], [0,0,1,1,0,0,1,1]) if v]) + '\n')
            f.write('Ic   0    c0  ' + pulse([p for p,v in zip(in_pulse[2], [0,0,0,0,1,1,1,1]) if v]) + '\n')
            f.write('Iclk 0   clk  ' + pulse(clk_pulse) + '\n')
            f.write('ROUT q0  0   1\n')

            f.write('.model jjmit jj( rtype=1, vg=2.6mV, cap=0.07pF, r0=160, rn=16, icrit=0.1mA)\n')

            # f.write('.print DEVI L_a_1\n.print DEVI LP_1_0\n.print DEVI LRB_1_0\n.print DEVI L_1_4\n.print DEVI L_5_8\n.end')
            # f.write('.print DEVV ROUT\n.end')
            f.write('.end')
        print(f'Simulating {file}')
        call(['josim-cli', '-o', csv_file, file, '-V', '1'])

    # time = []
    # opening the CSV file
    # josim-cli -o ./merge_test.csv ./merge_test.cir -V 1

        print(f'Reading {csv_file}')
        with open(csv_file) as f:
            # reading the CSV file
            csvFile = csv.reader(f)
            # displaying the contents of the CSV file
            _, *varnames = next(csvFile)
            # vars = {v:[] for v in varnames}
            vars = []
            for lines in csvFile:
                vars.append([float(q) for q in lines])
            vars = np.array(vars)

        vectors = {}
        time = vars[:,0]
        # fig, ax = plt.subplots()
        plotitems = [q for i,varname in enumerate(varnames, 1) for q in (np.array(time), vars[:,i])]
        # for i in range(9,len(plotitems),2):
        #     plotitems[i] /= 2*np.pi
        #     # plotitems[i] = np.remainder(plotitems[i] - plotitems[i][20], 2*np.pi)
        #     plotitems[i] = np.diff(plotitems[i])
        #     plotitems[i] = np.append(plotitems[i], 0)
        print(f'Plotting')
        p = ax.plot(*plotitems[:8], lw=1)
        pr = axr.plot(*plotitems[8:], lw=1)#, linestyle='dotted'
        ax.vlines( clk_pulse, ymin = 0, ymax = 0.5e-3, colors = 'r', linestyles='dashed', lw=0.5)
        axr.vlines(clk_pulse, ymin = 0, ymax = 0.5e-3, colors = 'r', linestyles='dashed', lw=0.5)
        # axr.vlines(clk_pulse, ymin = 0, ymax = 0.2*np.pi, colors = 'r', linestyles='dashed', lw=0.5)
        # for i,varname in enumerate(varnames, 1):
        #     vectors[varname] = vars[:,i]
        ax.grid()
        axr.grid()
        #     ax.plot(np.array(time), np.array(vectors[varname]), lw=1)
    # fig.legend([f'{_type}({elem})' for _type, elem in legends])
    plt.show()
