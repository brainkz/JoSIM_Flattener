'''
MUX schematic
'''
import csv
from matplotlib import pyplot as plt
import numpy as np
import os, sys
from collections import deque

from subprocess import call
from itertools import product

from elements import jj, ind, res, bias, pulse, boiler_plate, BiasCoef, dff, rdff, merge, buffer_branch, comparator, jtl


if __name__ == '__main__':

    file = 'buf_test.cir'
    base, ext = os.path.splitext(file)
    csv_file = f'{base}.csv'
    T = 240e-12
    tstep = 1e-12
    in_nodes = ['A','R',]
    out_nodes = ['q']
    clk = 't'
    nclk = 1

    analysis = {'dev' :['L_rdff_t2_rdff_to','L_rdff_q3_rdff_to'],
                'node':[], #['ab', 'cd',],
                'phase': ['B_rdff_t1_rdff_t2','B_rdff_to_0', 'B_rdff_to_rdff_b1'],
                }

    with open(file,'w') as f:
        f.write('*\n\n')
        rdff(f, 'rdff', 'A', 'R', clk, 'q')
        boiler_plate(f, in_nodes, clk, nclk, out_nodes, T = T, tstep = tstep, analysis = analysis, cycles = 2)
        f.write('.end')

    print(f'Simulating {file}')
    call(['josim-cli', '-o', csv_file, file, '-V', '1'])

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

def mux(f, name, set, reset, a, b, q,):

    # SET branch
    buffer_branch(f, f'{name}_set', set, f'{name}_s1', bias_pts = [0,3])
    bias(f, f'{name}_s2'              , IF = BiasCoef)
    jj(f,   f'{name}_s2', f'0'        , AF = 1)
    ind(f,  f'{name}_s2', f'{name}_s1', LF = 4)

    #RESET branch
    buffer_branch(f, f'{name}_res', reset, f'{name}_r1', bias_pts = 0)
    bias(f, f'{name}_r2'              , IF = BiasCoef)
    jj (f,  f'{name}_r2', f'0'        , AF = 1)
    ind(f,  f'{name}_r2', f'{name}_r1', LF = 2)

    ind(f,  f'{name}_s1', f'{name}_r1', LF = 4)

    sizes = 1.0, 1.2
    # branch A
    buffer_branch(f, f'{name}_a' ,            a, f'{name}_a1', bias_pts = 0, sizes  = {'jj2': 2})
    comparator(f, f'{name}_a1', f'{name}_s1', f'{name}_a2', *sizes)
    jj (f,  f'{name}_a2', f'0'        , AF = 1)
    jj (f,  f'{name}_a2', f'{name}_q4', AF = 1/1.4)

    # branch B
    buffer_branch(f, f'{name}_b' ,            b, f'{name}_b1', bias_pts = 0, sizes  = {'jj2': 2})
    comparator(f, f'{name}_b1', f'{name}_r1', f'{name}_b2', *sizes)
    jj (f,  f'{name}_b2', f'0'        , AF = 1)
    jj (f,  f'{name}_b2', f'{name}_q2', AF = 1/1.4)

    #output Branch
    ind(f,  f'{name}_q2', f'{name}_q1', LF = 1)
    bias(f, f'{name}_q1'              , IF = BiasCoef)

    jtl(f, f'{name}_qjtl', f'{name}_q1', q, LF = 2, AF = 1, njj = 2)


if __name__ == '__main__':
    file = 'mux_test.cir'
    base, ext = os.path.splitext(file)
    csv_file = f'{base}.csv'
    T = 240e-12
    tstep = 1e-12
    in_nodes = ['S','R','A','B']
    out_nodes = ['q']
    clk = 't'

    analysis = {'dev' :['L_main_s2_main_s1', 'L_main_r2_main_r1', 'L_main_s1_main_r1',
                        'LJCT_main_a1_0', 'LJCT_main_b1_0','LJCT_main_a2_0', 'LJCT_main_b2_0',
                        'LJCT_main_qjtl_1_0', 'L_main_q2_main_q1'],
                        # 'LJCT_main_a_2_main_a1', 'LJCT_main_a_1_0'],
                'node':[], #['ab', 'cd',],
                'phase': ['B_main_s1_main_a1', 'B_main_a_2_main_a1']
                }

    with open(file,'w') as f:
        f.write('*\n\n')
        mux(f, 'main', 'S', 'R', 'A', 'B', 'q',)
        # analysis =
        # boiler_plate(f, in_nodes, clk, nclk, out_nodes, T = 240e-12, analysis = {}, cycles = 1)
        boiler_plate(f, in_nodes, clk, 0, out_nodes, T = T, tstep = tstep, analysis = analysis)

        f.write('.end')

    print(f'Simulating {file}')
    call(['josim-cli', '-o', csv_file, file, '-V', '1'])

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
