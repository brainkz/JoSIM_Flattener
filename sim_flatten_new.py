#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Netlist preprocessing tool to support subcircuit variables in JoSIM.
The main function of this module parses the input netlist, replaces
variables top down, and replaces subcircuit instances flattening the netlist.

The pyparsing module is used to correctly and safely replace the parameters
within the netlist

Example:
    The following netlist is flattened using this tool.

Todo:
    * Add control of verbosity
    * Run additional tests
    * Support local jj models within the subcircuits
"""

import os
import re
import sys
from argparse import ArgumentParser

from itertools import chain
from collections import deque, defaultdict
from subprocess import call
from pprint import pprint

from pyparsing import (
    alphanums,
    alphas,
    ParseException,
    Word,
)

from fourFn import (
    BNF,
    evaluate_stack,
    exprStack,
)


WS = re.compile(r'\s+')
EQ = re.compile(r'\=')
SUFFIXED_NUM = re.compile(r'(\b(\d+)[FPNUMKXGT]{1}\b)|(\b(\d+)MEG\b)')
ALPHA = re.compile(r'[A-Z]+')
BRACKET = re.compile(r'(?<=\().+(?=\))')

ALLOWED_NODE_CHARS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890|_'
TWONODE_DEVICES = tuple('RCLKEFHG')
FOURNODE_DEVICES = tuple('EFHGT')
SOURCES = tuple('VIP')
SUFFIX = {  'F':'E-15','P'  :'E-12','N':'E-9','U':'E-6','M':'E-3' ,
            'K':'E+3' ,'MEG':'E+6' ,'X':'E+6','G':'E+9','T':'E+12',}


CONTROLS = ('.TRAN', '.TEMP', '.NEB', '.SPREAD',
            '.IV', '.PRINT', '.PLOT', '.SAVE',
            '.FILE', '.MODEL')

CONSTANTS = {   'PI'        : 3.141592653589793238463,
                'PHI_ZERO'  : 2.067833831170082E-15,
                'BOLTZMANN' : 1.38064852E-23,
                'EV'        : 1.6021766208E-19,
                'HBAR'      : 1.0545718001391127E-34,
                'C'         : 299792458,
                'MU0'       : 12.566370614E-7,
                'EPS0'      : 8.854187817E-12,
                'SIGMA'     : 3.291059757E-16}

def parse_expr(input_string: str, variables: dict, debug_flag: bool = False):
    """ Safely parse mathematical expression. Based on example fron pyparsing
    module.

    Args:
        input_string: string containing the expression.
        variables: dictionary of parameters visible to the device consiting of
        global parameters and local parameters defined within the subckt
        definition.
        debug_flag: flag to enable printing of additional information

    Returns:
        A 3-tuple in format (status, name, value) where
            * status -- indicates whether parsing is successful
            * name -- variable name if the expression is in form {name}={value}
            * value -- evaluation result

    Example:
        the input parameters
            input_string = "A = 2*B"
            variables = {'B':10}
        produce
            (True, 'A', 20)
    """
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
        # print(f'SUBSTITUTING THE SUFFIXES : {input_string}')
        # if s[:item.start()] + '_____' + s[item.end():]
        # item = SUFFIXED_NUM.search('B1341P*MEG242MEG*23.149P/3242U*23.55MEG')
    # try parsing the input string
    try:
        L = pattern.parseString(input_string, parseAll=True)
    except ParseException as err:
        L = ["Parse Failure", input_string, (str(err), err.line, err.column)]
        # print(L)
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

        try:
            result = evaluate_stack(exprStack)
        except ValueError as e:
            print(e)
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

def strip_uncomment_upper_join(file : str):
    ''' Iterator over lines in {file}. Strips trailing whitespace. Skips fully
    commented lines. Converts line to uppercase. Joins lines separated with (+)

    Args:
        file: path to a file for generate the lines

    Yields:
        Uncomented lines from {file} in uppercase, free from whitespace and
        line breaks.
    '''
    with open(file, 'r', encoding = 'utf-8') as f:
        it = iter(f)
        lines = []
        for line in it:
            line = line.strip().upper()
            if not line or line.startswith('*'):
                continue
            elif line.startswith('+'):
                lines[-1] += ' ' + line[1:]
            else:
                lines.append(line)
        
        # breakpoint()
        for line in lines:
            line = line.strip().upper()
            if line.startswith('*') or not line:
                continue
            elif line.startswith('.INCLUDE'):
                _, include_path = WS.split(line, 1)
                yield from strip_uncomment_upper_join(include_path)
            elif line == '.END':
                break
            elif line.startswith('+'):
                line += ' ' + line[1:]
            else:
                yield line

def subckt_line(line: str):
    _, subc_name, *io_params = line.split()
    io = []
    kwargs = {}
    for i, item in enumerate(io_params):
        if '=' in item:
            kwargs.update(EQ.split(expr) for expr in io_params[i:])
            break
        else:
            io.append(item)
    return subc_name, io, kwargs

def expr2dict(expr):
    if isinstance(expr, str):
        name, val = EQ.split(expr)
        return {name: val}
    if isinstance(expr, list):
        return dict(EQ.split(ex) for ex in expr)

def parse_twonode(line: str):
    device_name, *inst_io, val_str = line.split()
    _type, label = device_name[0], device_name[1:]
    return _type, label, inst_io, val_str

def parse_fournode(line: str):
    device_name, n1, n2, n3, n4, *expressions = line.split()
    assert(all('=' in item for item in expressions))
    _type, label = device_name[0], device_name[1:]
    return _type, label, (n1, n2, n3, n4), expr2dict(expressions)

def parse_jj(line: str):
    inst_name, n1, n2, *items = line.split()

    # Resolve whether the third node is provided
    expressions = []
    if   '=' in items[0]:
        model = ''
        expressions.extend(items)
    elif '=' in items[1]:
        model = items[0]
        expressions.extend(items[1:])
    else:
        model = items[1]
        expressions.extend(items[2:])
    assert(all('=' in item for item in expressions))
    inst_io = (n1, n2)
    label = inst_name[1:]
    return label, inst_io, model, expr2dict(expressions)

def parse_source(line: str):
    device_name, n1, n2, param_str = WS.split(line, 3)
    _type = device_name[0]
    label = device_name[1:]
    return _type, label, (n1, n2), param_str

def parse_subckt_inst(line: str, subckts: dict):
    inst_name, *subc_name_and_io_and_expr = line.split()
    label = inst_name[1:]

    # First, extract the parameter expressions
    subc_name_and_io = []
    expressions = []
    for i, item in enumerate(subc_name_and_io_and_expr):
        if '=' not in item:
            subc_name_and_io.append(item)
        else:
            expressions.extend(subc_name_and_io_and_expr[i:])
            break

    # Resolve whether the model name is supplied first or last
    if subc_name_and_io[0] in subckts: #JSIM mode
        subc_name = subc_name_and_io[0]
        inst_io = subc_name_and_io[1:]
    elif subc_name_and_io[-1] in subckts:
        subc_name = subc_name_and_io[-1]
        inst_io = subc_name_and_io[:-1]
    else:
        raise NameError('SUBCKT declaration syntax not recognized')
    n_io = len(subckts[subc_name]['io'])
    assert(len(inst_io) == n_io)
    return label, subc_name, inst_io, expr2dict(expressions)

def default_kw(line_iter: str = None, params: dict = {}, subckt_name: str = 'main'):
    expr_dict = {}
    if isinstance(line_iter, str):
        line_iter = strip_uncomment_upper_join(line_iter)
    for line in line_iter:
        if line.startswith('.PARAM'):
            _, expr = line.split()
            name, val_str = EQ.split(expr)
            expr_dict[name] = val_str
        elif line.startswith('.SUBCKT'):
            new_subc_name, io, kwargs = subckt_line(line)
            params[new_subc_name] = {'kw': kwargs, 'io': io}
            params = default_kw(line_iter, params, new_subc_name)
        elif line.startswith('.END'):
            break
    params.setdefault(subckt_name,{'kw': [], 'io': []}).update({'expr_dict':expr_dict})
    return params


def get_devices(line_iter: str, params: dict, subckt_name: str = 'main'):
    
    def cprint(*args, **kwargs):
        if subckt_name == 'main':
            print(*args, **kwargs)
    params[subckt_name]['devices'] = []
    params[subckt_name]['instances'] = []
    if isinstance(line_iter, str):
        line_iter = strip_uncomment_upper_join(line_iter)
    
    for line in line_iter:
        cprint(line[0])
        if line.startswith(TWONODE_DEVICES):
            _type, label, io, val_str = parse_twonode(line)
            cprint(f"TWO_NODE: {line}, {_type, label, io, val_str}")
            params[subckt_name]['devices'].append((_type, label, io, '', val_str))
        elif line.startswith(FOURNODE_DEVICES):
            if line.startswith('T'):
                continue
            _type, label, io, expr = parse_fournode(line)
            cprint(f"FOUR_NODE: {line}, {_type, label, io, expr}")
            params[subckt_name]['devices'].append((_type, label, io, '', expr))
        elif line.startswith('B'):
            label, io, model, expr = parse_jj(line)
            cprint(f"JJ: {line}, {label, io, model, expr}")
            params[subckt_name]['devices'].append(('B', label, io, model, expr))
        elif line.startswith(SOURCES):
            _type, label, io, expr = parse_source(line)
            cprint(f"SOURCE: {line}, {_type, label, io, expr}")
            params[subckt_name]['devices'].append((_type, label, io, '', expr))
        elif line.startswith('X'):
            label, inst_subc_name, io, expr = parse_subckt_inst(line, params)
            cprint(f"SUBCKT: {line}, {label, inst_subc_name, io, expr}")
            params[subckt_name]['instances'].append((label, inst_subc_name, io, expr))
        elif line.startswith('.SUBCKT'):
            _, new_subckt_name, *_ = line.split()
            cprint(f"SUBCKT DEF: {line}, {new_subckt_name}")
            params = get_devices(line_iter, params, new_subckt_name)
        # elif subckt_name == 'main':
        #     print(f'COULD NOT MATCH:\n{line}')
        #     breakpoint()
        if line.startswith(('.END', '.ENDS')):
            break
        
    return params

def parse_controls(file):
    line_iter = strip_uncomment_upper_join(file)
    control_lines = []
    for line in line_iter:
        if line.startswith(CONTROLS):
            control_lines.append(line)
        elif line.startswith('.SUBCKT'):
            for line in line_iter:
                if line.startswith('.ENDS'):
                    break
        elif line.startswith('.CONTROL'):
            for line in line_iter:
                if line.startswith('.ENDC'):
                    break
        elif line.startswith('T'):
            control_lines.append(line)
    return control_lines

def evaluate(expr, known = {}, copy = False):
    if copy:
        known = known.copy()
    if isinstance(expr, dict):
        for name, val_str in expr.items():
            if isinstance(val_str, str):
                success, _, val = parse_expr(val_str, known)
                if success:
                    known[name] = val
            else:
                known[name] = val_str
        return known
    elif isinstance(expr, str):
        success, _, val = parse_expr(expr, known)
        return val if success else expr
    else:
        return expr

def write_flat_file(params: dict, commands: list, temp_file: str, global_params: dict,
                    measured_devices: int = 999, measured_nodes: int = 999, measured_phases: str = 999):
    ''' Write flattened netlist {flat_file} based on the heirarchical
    information stored in {subckts}. Append {commands} and write {params}.

    Args:
        subckts : dictionary containing subcircuit parameters, including I/O,
        local keyword parameters, devices and subcircuits instantiated within
        the subcircuit, dependencies w.r.t. other subcircuits, and raw lines
        within the subcircuit.
        commands : list of commands to write within the netlist
        temp_file : temporary file for storing the flattened netlist
        params : dictionary containing the parameters keyed by the subcircuit
        name. Global parameters are placed in params['main']
    '''
    with open(temp_file, 'w') as fobj:
        fobj.write('*\n')
        for name, value in global_params.items():
            # print(f'.PARAM {name}={value}')
            fobj.write(f'.PARAM {name}={value}\n')
        for line in commands:
            # print(line)
            fobj.write(line + '\n')

        for _type, label, nodes, model, value in params['main']['devices']:
            if isinstance(value, dict):
                value_str = ' '.join(f'{k}={v}' for k,v in value.items())
            else:
                value_str = value
            # print(f'{_type}{label} {" ".join(nodes)} {model} {value_str}')
            fobj.write(f'{_type}{label} {" ".join(nodes)} {model} {value_str}\n')
        
        for _type, name, *_ in params['main']['devices']:
            if name.count('|') <= measured_devices:
                fobj.write(f'.print DEVI {_type}{name}\n')
        nodes = set().union(*[q[2] for q in params['main']['devices']])
        nodes.remove('0')
        for node in nodes:
            if node.count('|') <= measured_nodes:
                fobj.write(f'.print V {node}\n')
        for _type, name, *_ in params['main']['devices']:
            if _type == 'B' and name.count('|') <= measured_phases:
                fobj.write(f'.print DEVP B{name}\n')
                

def convert_sim(file, temp_file = None, csv_path = None):
    folder, name = os.path.split(file)
    base, ext = os.path.splitext(name)
    if temp_file is None:
        temp_file = f'{base}_temp{ext}'
    if csv_path is None:
        csv_path = f'{base}.csv'

    # initially, just get the parameters
    params = default_kw(file)
    # then save the control statements
    control_lines = parse_controls(file)
    # next, parse all subckt structures
    params = get_devices(file, params)
    # swap each subckt instance with its equivalent representation
    global_params = evaluate(params['main']['expr_dict'], {})
    while params['main']['instances']:
        new_instances = []
        for label, parent, inst_io, inst_expr in params['main']['instances']:
            model_io = params[parent]['io']
            io_map = dict(zip(model_io, inst_io)) # map child_io to parent io
            io_map['0'] = '0'

            known = evaluate(params[parent]['kw'] | params[parent]['expr_dict'] | inst_expr, global_params, copy=True)
            for _type, dev_label, dev_io, model, expr in params[parent]['devices']:
                dev_label = f'{label}|{dev_label}'
                dev_io = [io_map.get(n,f'{label}|{n}') for n in dev_io]
                if isinstance(expr, dict):
                    necessary_keys = list(expr.keys())
                    expr = evaluate(expr, known, copy=True)
                    expr = {k:v for k,v in expr.items() if k in necessary_keys}
                elif expr.startswith('PWL'):
                    # breakpoint()
                    expr_str = re.search(r'\(.+\)', expr)[0][1:-1] # get what's inside the brackets
                    expr = ' '.join(evaluate(e, known, copy=True) for e in expr_str.strip().split())
                    expr = 'PWL(' + expr + ')'
                else:
                    expr = evaluate(expr, known, copy=True)
                params['main']['devices'].append((_type, dev_label, dev_io, model, expr))
            for child_label, child, child_io, child_expr in params[parent]['instances']:
                child_label = f'{label}|{child_label}'
                child_io = [io_map.get(n,f'{label}|{n}') for n in child_io]
                child_expr = evaluate(params[child]['expr_dict']|params[child]['kw']|child_expr, known, copy=True)
                new_instances.append((child_label, child, child_io, child_expr))
                # break
        params['main']['instances'] = new_instances

    write_flat_file(params, control_lines, temp_file, global_params)
    status = call(['josim-cli', '-o', csv_path, temp_file, '-V', '1'])
    return csv_path


if __name__ == '__main__':
    # file = 'async_gate_test.cir'
    # file = 'A_and_B_xor_C.cir'
    # file = 'A_and_B_xor_C_and_D.cir'
    # file = 'two_gates_seq.cir'
    # file = 'cb_spurious.cir'
    # file = 'desync.cir'
    # file = 'NIMPLY.cir'
    # file = 'desync_2.cir'
    # file = 'bka_1.cir'
    # file = 'bka_2.cir'
    # file = 'bka_3.cir'
    # file = 'bka_component_test.ckt'
    # file = 'Margins_0001.ckt'
    # file = 'bka_3_ptl_testvectors_2.cir'
    # file = 'bka_3_ptl_MITLL_cells.cir'
    # file = 'bka_8bit.cir'
    # file = 'TFF.cir'
    # file = 'tff_fa.cir'
    # file = 'tff_fa_working.cir'
    # file = 'tff_fa_fix.cir'
    file = 'tff_fa_RSFQlib.cir'
    # file = 'pulse_code.cir'
    
    lvl = 1
    
    base, ext = os.path.splitext(file)
    temp_file = f'{base}_temp{ext}'
    csv_path = f'{base}.csv'

    print('PARAMETERS')
    line_iter = strip_uncomment_upper_join(file)
    params = default_kw(line_iter)
    control_lines = parse_controls(file)
    print('DEVICES')
    line_iter = strip_uncomment_upper_join(file)
    params = get_devices(line_iter, params)
    global_params = evaluate(params['main']['expr_dict'], {})
    while params['main']['instances']:
        new_instances = []
        for label, parent, inst_io, inst_expr in params['main']['instances']:
            model_io = params[parent]['io']
            io_map = dict(zip(model_io, inst_io)) # map child_io to parent io
            io_map['0'] = '0'

            known = evaluate(params[parent]['kw'] | params[parent]['expr_dict'] | inst_expr, global_params, copy=True)
            for _type, dev_label, dev_io, model, expr in params[parent]['devices']:
                dev_label = f'{label}|{dev_label}'
                dev_io = [io_map.get(n,f'{label}|{n}') for n in dev_io]
                if isinstance(expr, dict):
                    necessary_keys = list(expr.keys())
                    expr = evaluate(expr, known, copy=True)
                    expr = {k:v for k,v in expr.items() if k in necessary_keys}
                elif expr.startswith('PWL'):
                    expr_str = re.search(r'\(.+\)', expr)[0][1:-1] # get what's inside the brackets
                    expr = 'PWL(' + ' '.join(f'{evaluate(e, known, copy=True):g}' for e in expr_str.strip().split()) + ')'
                else:
                    expr = evaluate(expr, known, copy=True)
                params['main']['devices'].append((_type, dev_label, dev_io, model, expr))
            for child_label, child, child_io, child_expr in params[parent]['instances']:
                child_label = f'{label}|{child_label}'
                child_io = [io_map.get(n,f'{label}|{n}') for n in child_io]
                child_expr = evaluate(params[child]['expr_dict']|params[child]['kw']|child_expr, known, copy=True)
                new_instances.append((child_label, child, child_io, child_expr))
                # break
        params['main']['instances'] = new_instances

    write_flat_file(params, control_lines, temp_file, global_params, lvl, lvl, lvl+2)
    
    print(f"FILE {temp_file} WRITTEN")
    
    status = call(['josim-cli', '-o', csv_path, temp_file, '-V', '1'])
    # status = call(['josim-cli', '-o', csv_path, temp_file, '-V', '1'])

            # inst_expr = inst_expr.extend(params[inst_subc_name]['expr_deque'])
            # inst_expr = inst_expr.extend(params[inst_subc_name]['kw'])

            # breakpoint()
