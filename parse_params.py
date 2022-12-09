


import re
WS = re.compile(f'\s+')
EQ = re.compile(f'\=')
SUFFIXED_NUM = re.compile(r'(\b(\d+)[FPNUMKXGT]{1}\b)|(\b(\d+)MEG\b)')
ALPHA = re.compile(f'[A-Z]+')

ALLOWED_NODE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890|_'
DEVICES = tuple('RCLKEFHG')
ADVANCED_DEVICES = tuple('BT')
SOURCES = tuple('VIP')
SUFFIX = {  'F':'E-15','P'  :'E-12','N':'E-9','U':'E-6','M':'E-3' ,
            'K':'E+3' ,'MEG':'E+6' ,'X':'E+6','G':'E+9','T':'E+12',}

import os
import ast
import sys
from collections import defaultdict
from pprint import pprint

# Simple expression parse base on SimpleCalc.py from the
# pyparsing module.
# Supports out-of-order definitions.

from pyparsing import ParseException, Word, alphas, alphanums
from collections import deque
from fourFn import BNF, exprStack, evaluate_stack

def parse_device(line):
    name, *inst_io, value = WS.split(line)
    _type = name[0]
    label = name[1:]
    return _type, label, inst_io, value

def parse_subckt_inst(line, subckts, global_params):
    inst_name, subc_name, *io_and_params = WS.split(line)
    label = inst_name[1:]
    subc_io = subckts[subc_name]['io']
    n_io = len(subc_io)
    inst_io = io_and_params[:n_io]
    expressions  = io_and_params[n_io:]
    local_params = {}
    for expr in expressions:
        param, value = EQ.split(expr)
        local_params[param] = value
    return label, subc_name, inst_io, local_params

def parse_jj(line):
    # B5       14  3   jjmitll100 area=1.5
    inst_name, n1, n2, model, *expressions = WS.split(line)
    inst_io = (n1, n2)
    label = inst_name[1:]
    local_params = {}
    for expr in expressions:
        param, value = EQ.split(expr)
        local_params[param] = value
    return label, inst_io, model, local_params

def parse_source(line):
    name, n1, n2, param_str = WS.split(line, 3)
    _type = name[0]
    label = name[1:]
    return _type, label, (n1, n2), '', param_str


def parse_expr(input_string, variables, debug_flag = False):
    # Reset to an empty exprStack
    del exprStack[:]

    arithExpr = BNF()
    ident = Word(alphas, alphanums).setName("identifier")
    assignment = ident("varname") + "=" + arithExpr
    pattern = assignment | arithExpr

    # print(input_string)
    while item := SUFFIXED_NUM.search(input_string):
        suf = ALPHA.search(item[0])[0]
        sub_str = item[0].replace(suf, SUFFIX[suf])
        input_string = input_string[:item.start()] + sub_str + input_string[item.end():]
        print(f'SUBSTITUTING THE SUFFIXES : {input_string}')
        # if s[:item.start()] + '_____' + s[item.end():]
        # item = SUFFIXED_NUM.search('B1341P*MEG242MEG*23.149P/3242U*23.55MEG')
    # try parsing the input string
    try:
        L = pattern.parseString(input_string, parseAll=True)
    except ParseException as err:
        L = ["Parse Failure", input_string, (str(err), err.line, err.column)]
        print(L)
        return False, '', ''

    # show result of parsing the input string
    if debug_flag:
        print(input_string, "->", L)
    if len(L) == 0 or L[0] != "Parse Failure":
        if debug_flag:
            print("exprStack=", exprStack)

        for i, ob in enumerate(exprStack):
            if isinstance(ob, str) and ob in variables:
                exprStack[i] = str(variables[ob])

        # calculate result , display the result to user
        try:
            result = evaluate_stack(exprStack)
        except Exception as e:
            return False, '', ''
        else:
            return True, L.varname, result
            # Assign result to a variable if required
            # if L.varname:
            #     variables[L.varname] = result
            # if debug_flag:
            #     print("variables=", variables)
    else:
        print("Parse Failure")
        err_str, err_line, err_col = L[-1]
        print(err_line)
        print(" " * (err_col - 1) + "^")
        print(err_str)

def strip_uncomment_upper_join(file):
    ''' Iterator over lines in file f.
    Strips trailing whitespace
    Skips fully commented lines
    Converts line to uppercase
    Joins separated lines
    '''
    with open(file, 'r') as f:
        for prev_line in f:
            prev_line = prev_line.strip().upper()
            if not prev_line.startswith('*') or not prev_line:
                break
        if prev_line.startswith('+'):
            raise ValueError('First uncommented line starts with "+"')
        for next_line in f:
            next_line = next_line.strip().upper()
            if next_line.startswith('*') or not next_line:
                continue
            elif next_line == '.END':
                yield prev_line
                yield next_line
                break
            elif next_line.startswith('+'):
                prev_line += ' ' + next_line[1:]
            else:
                yield prev_line
                prev_line = next_line

def parse_params(file):
    params = {}
    expressions = deque()
    for line in strip_uncomment_upper_join(file):
        if line.startswith('.PARAM'):
            _, expr = WS.split(line, 1)
            expressions.appendleft(expr)

    while expressions:
        input_string = expressions.pop()
        print(input_string)
        success, varname, value = parse_expr(input_string, params)
        if success:
            params[varname] = value
        else:
            # Some dependency exists. Parse the string later
            expressions.appendleft(input_string)
    return params

def topsort_subckts(subckts):
    order = []
    dependencies = {sub: d['dependencies'] for sub,d in subckts.items()}
    while dependencies:
        print(dependencies)
        for sub, dep in dependencies.items():
            if not dep and sub not in order:
                order.append(sub)
        dependencies = {sub: dep for sub,dep in dependencies.items() if sub not in order}
        for sub, dep in dependencies.items():
            if sub not in order:
                dependencies[sub] = dependencies[sub].difference(order)
    return order

def parse_subckt_params(file, global_params):
    models = {}
    # First pass
    other_lines = []
    subckts = {'main':{      'name': 'main',
                               'io': [],
                     'local_params': {},
                     'devices' : [],
                     'dependencies': set(),
                 'subckt_inst': [],
                        'raw_lines': [] }}
    line_iter = strip_uncomment_upper_join(file)
    for line in line_iter:
        if line.startswith('.SUBCKT'):
            _, sub_name, *io_and_params = WS.split(line)
            sub_dict = {'name' : sub_name,
                        'io': [],
                        'local_params': {},
                        'dependencies': set(),
                        'devices' : [],
                'subckt_inst': [],
                        'raw_lines': []}
            for item in io_and_params:
                if all(c in ALLOWED_NODE_CHARS for c in item):
                    sub_dict['io'].append(item)
                elif '=' in item:
                    local_param, default = EQ.split(item)
                    sub_dict['local_params'][local_param] = default
                else:
                    raise ValueError(f'Unsupported syntax at the following line:\n{line}\nitem:{item}')
            for line in line_iter:
                if line.startswith('.ENDS'):
                    break
                else:
                    sub_dict['raw_lines'].append(line)
                if line.startswith('X'):
                    label, parent_sub, *_ = WS.split(line)
                    sub_dict['dependencies'].add(parent_sub)
                # elif line.startswith(DEVICES):
                # elif line.startswith(DEVICES + ADVANCED_DEVICES):
                # elif line.startswith(SOURCES):
            subckts[sub_name] = sub_dict

        elif line.startswith('X'):
            subckts['main']['raw_lines'].append(line)
            label, parent_sub, *_ = WS.split(line)
            subckts['main']['dependencies'].add(parent_sub)
        elif line.startswith(DEVICES + ADVANCED_DEVICES):
            subckts['main']['raw_lines'].append(line)
        else:
            other_lines.append(line)

    order = topsort_subckts(subckts)

    for subc in order[::-1]: # propagate parameters top down
        # print(subc)
        variables = global_params.copy()
        variables.update(subckts[subc]['local_params'])
        for line in subckts[subc]['raw_lines']:
            print(line)
            if line.startswith(DEVICES):
                _type, label, inst_io, params = parse_device(line)
                print(line)
                success, _, value = parse_expr(params, variables)
                if not success:
                    breakpoint()
                # assert(success)
                subckts[subc]['devices'].append((_type, label, inst_io, '', value))
                # print(value)
            elif line.startswith('X'):
                label, subc_name, inst_io, custom_params = parse_subckt_inst(line, subckts, global_params)
                print(line)
                resolved_params = subckts[subc_name]['local_params'].copy()
                resolved_params.update(custom_params)
                for param, val_str in resolved_params.items():
                    success, _, val = parse_expr(val_str, variables)
                    if success:
                        resolved_params[param] = val
                    else:
                        print(f'Failed to propagate parameter {val_str}')
                subckts[subc]['subckt_inst'].append((label, subc_name, inst_io, resolved_params))
            elif line.startswith('B'):
                label, inst_io, model, params = parse_jj(line)
                # print(line)
                for param, val_str in params.items():
                    success, _, val = parse_expr(val_str, variables)
                    if success:
                        params[param] = val
                    else:
                        print(f'Failed to propagate parameter {val_str}')
                subckts[subc]['devices'].append(('B', label, inst_io, model, params))
            elif line.startswith(SOURCES):
                _type, label, inst_io, _, param_str = parse_source(line)
                subckts[subc]['devices'].append((_type, label, inst_io, '', param_str))

            # elif line.startswith('T'):
            #     label, inst_io, params = parse_jj(line)
            #     print(line)
            #     for param, val_str in params.items():
            #         success, _, val = parse_expr(val_str, variables)
            #         if success:
            #             params[param] = val
            #         else:
            #             print(f'Failed to propagate parameter {val_str}')
            #     subckts[subc]['devices'].append(('B', label, inst_io, params))

    for parent in order: # replace subckt definitions bottom up
        for x_label, child, inst_io, resolved_params in subckts[parent]['subckt_inst']:
            model_io = subckts[child]['io']
            io_map = dict(zip(model_io, inst_io)) # map child_io to parent io
            io_map['0'] = '0'
            for _type, model_label, model_nodes, model, value in subckts[child]['devices']:
                inst_nodes = [io_map.setdefault(n, f'{n}|{x_label}') for n in model_nodes]
                inst_label = f'{model_label}|{x_label}'
                subckts[parent]['devices'].append((_type, inst_label, inst_nodes, model, value))

    return subckts, order, other_lines

def write_temp_file(subckts, other_lines, temp_file):
    with open(temp_file, 'w') as f:
        f.write('*\n')
        for _type, label, nodes, model, value in subckts['main']['devices']:
            if isinstance(value, dict):
                value_str = ' '.join(f'{k}={v}' for k,v in value.items())
            else:
                value_str = value
            f.write(f'{_type}{label} {" ".join(nodes)} {model} {value_str}\n')
        for line in other_lines:
            f.write(line + '\n')


{'FUNDAMENTAL_R': {'dependencies': set(),
                   'devices': [],
                   'io': ['1', '2'],
                   'local_params': {'RES': '1'},
                   'name': 'FUNDAMENTAL_R',
                   'raw_lines': ['R1     1 2 RES'],
                   'subckt_inst': []},
 'JJ_GND': {'dependencies': {'FUNDAMENTAL_R'},
            'devices': [],
            'io': ['1', '2'],
            'local_params': {'AF': '1', 'LLP': '2E-13', 'LP_ADD': '0'},
            'name': 'JJ_GND',
            'raw_lines': ['LP  MID1  JCT LLP',
                          'X1  FUNDAMENTAL_R   1 MID2 RES=6.859904418/(2.5*AF)',
                          'LRB MID2  JCT (3.875845996E-12/(2.5*AF)+LP_ADD*0.2E-12)',
                          'LJCT JCT    2 1E-18'],
            'subckt_inst': []},}



if __name__ == '__main__':
    # file = 'ignore/josim_test.cir'
    # temp_file = 'ignore/josim_test_temp.cir'

    path = '/Users/brainkz/Documents/GitHub/JoSIM/test/ex_jtl_string.cir'
    # path = 'jtl_test.cir'
    folder, name = os.path.split(path)
    base, ext = os.path.splitext(name)

    temp_file = f'ignore/{base}_temp.cir'
    csv_file = f'ignore/{base}.csv'
    csv_file_temp = f'ignore/{base}_temp.csv'

    global_params = parse_params(path)
    subckts, order, other_lines = parse_subckt_params(path, global_params)
    write_temp_file(subckts, other_lines, temp_file)
    from subprocess import call
    call(['open', temp_file])
    call(['josim-cli', '-o', csv_file, path, '-V', '1'])
    call(['josim-cli', '-o', csv_file_temp, temp_file, '-V', '1'])
    # call(['/Applications/BeSpiceWave.app/Contents/MacOS/BeSpiceWave', csv_file_temp])
    # pprint(subckts)
    #
    # for child in order:
    #     for parent, cfg in subckts.items():
    #         if child in cfg['dependencies']:




    # next_line = next(f, '').strip().upper()
    # while next_line.startswith('+'):
    #     prev_line += next_line[1:]
    #     next_line = next(f, '').strip().upper()
    # prev_line
