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

from collections import deque, defaultdict
from subprocess import call

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


def parse_subckt_def(line: str, line_iter, subckts: dict):
    ''' Parse subcircuit definition and update the subckts dictionary

    Args:
        line : string containing the subckt definition. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        line_iter : line generator used during parsing. The generator will
        advance until a line starting with ".ENDS" is encountered.
        subckts : dictionary containing subcircuit parameters. Updated during
        this function execution
    '''
    _, sub_name, *io_and_params = line.split()
    for i, item in enumerate(io_and_params):
        if '=' in item:
            subckts[sub_name]['default_params'] = dict(EQ.split(param) for param in io_and_params[i:])
            break
        else:
            subckts[sub_name]['io'].append(item)
    # assert(all(c in ALLOWED_NODE_CHARS for c in node for node in subckts[sub_name]['io']))
    for line in line_iter:
        if line.startswith('.ENDS'):
            break
        else:
            subckts[sub_name]['raw_lines'].append(line)
        if line.startswith('X'):
            subc_name = parse_subckt_inst(line, subckts, just_model = True)
            subckts[sub_name]['dependencies'].add(subc_name)
    return subckts

def parse_twonode(line: str, variables: dict):
    ''' Parse the two terminal device line in format
    {type}{label} {node+} {node-} {value}

    Args:
        line: string containing the device definition. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        variables: dictionary of parameters visible to the device consiting of
        global parameters and local parameters defined within the subckt
        definition.

    Returns:
        A 4-tuple in format (_type, label, inst_io, value) where
            * _type -- a single character representing type of the device
            * label -- the label of the device
            * inst_io -- a tuple of device nodes
            * value -- string representing the primary parameter of the device,
            i.e., resistance of the resistor, or capacitance of the capacitor.
            Stored in string format to support expressions.

    Example:
        a string
            "R1 in out 3*scale/k"
        is parsed as
            "R", "1", ("in", "out"), "3*scale/k"
        The expression "3*scale/k" is evaluated using the variables stored
        in "variables". If evaluation fails, the expression is returned as is.
    '''
    device_name, *inst_io, val_str = line.split()
    _type, label = device_name[0], device_name[1:]
    success, _, val = parse_expr(val_str, variables)
    if success:
        return _type, label, inst_io, val
    else:
        print(f'Failed to propagate parameter {val_str}')
        return _type, label, inst_io, val_str

def parse_subckt_inst(line: str, subckts: dict, just_model: bool = False):
    ''' Parse the subckt instance declaration line in format
    X{label} {node_1} ... {node_N} {param_1}={value_1} ... {param_M}={value_M}

    Args:
        line: string containing the subcircuit instance. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        subckts:  dictionary containing parameters of each known subcircuit
    Returns:
        A 4-tuple in format (label, subc_name, inst_io, kw_params) where
            * label -- the label of the subckt
            * subc_name -- the name of the subckt
            * inst_io -- a tuple of subckt nodes
            * kw_params -- dictionary mapping the parameters to declared
            values. Values are stored as strings to support parameters

    Example:
        a string
            "X02 JTLSTRING1000 in out size=1.3 factor=10 scale=2"
        is parsed as
            "02", "JTLSTRING1000", ("in", "out"), {"size": "1.3", "factor": "10", "scale": "2"}
    '''
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
    if (subc_name:= subc_name_and_io[0]) in subckts: #JSIM mode
        inst_io = subc_name_and_io[1:]
    elif (subc_name:= subc_name_and_io[-1]) in subckts:
        inst_io = subc_name_and_io[:-1]
    else:
        raise NameError('SUBCKT declaration syntax not recognized')

    if just_model:
        return subc_name
    else:
        subc_io = subckts[subc_name]['io']
        n_io = len(subc_io)
        assert(len(inst_io) == n_io)
        kw_params = {}
        for expr in expressions:
            param, value = EQ.split(expr)
            kw_params[param] = value
        return label, subc_name, inst_io, kw_params

def parse_jj(line: str, variables: dict):
    ''' Parse the jj instance declaration line in format
    B{label} {node+} {node-} {model} {param_1}={value_1} ... {param_M}={value_M}

    Args:
        line: string containing the subcircuit instance. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        variables: dictionary of parameters visible to the device consiting of
        global parameters and local parameters defined within the subckt
        definition.

    Returns:
        A 4-tuple in format (label, inst_io, model, local_params) where
            * label -- the label of the subckt
            * inst_io -- a tuple of subckt nodes
            * subc_name -- the name of the subckt model
            * local_params -- dictionary mapping the parameters to declared
            values. Values are stored as strings to support parameters

    Example:
        a string
            "B02 JTLSTRING1000 in out area=1.3 temp=10 freq=2G"
        is parsed as
            "02", "JTLSTRING1000", ("in", "out"), {"area": "1.3", "temp": "10", "freq": "2G"}
        The expressions are evaluated using the variables stored in "variables".
        If evaluation fails, the expression is returned as is.
    '''
    inst_name, n1, n2, *phase_node_and_model_and_expressions = line.split()

    # Resolve whether the third node is provided
    expressions = []
    if '=' in phase_node_and_model_and_expressions[0]:
        model = ''
        expressions.extend(phase_node_and_model_and_expressions)
    else:
        if '=' in phase_node_and_model_and_expressions[1]:
            model = phase_node_and_model_and_expressions[0]
            expressions.extend(phase_node_and_model_and_expressions[1:])
        else:
            model = phase_node_and_model_and_expressions[1]
            expressions.extend(phase_node_and_model_and_expressions[2:])
    assert(all('=' in item for item in expressions))

    inst_io = (n1, n2)
    label = inst_name[1:]
    local_params = {}
    for expr in expressions:
        name, val_str = EQ.split(expr)
        success, _, val = parse_expr(val_str, variables)
        if success:
            local_params[name] = val
        else:
            print(f'Failed to propagate parameter {val_str}')
            local_params[name] = val_str
    return label, inst_io, model, local_params

def parse_fournode(line: str, variables: dict):
    ''' Parse the two terminal device line in format
    {type}{label} {node1} {node2} {node3} {node4} {parameters}

    Args:
        line: string containing the device definition. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        variables: dictionary of parameters visible to the device consiting of
        global parameters and local parameters defined within the subckt
        definition.

    Returns:
        A 4-tuple in format (_type, label, inst_io, value) where
            * _type -- a single character representing type of the device
            * label -- the label of the device
            * inst_io -- a tuple of device nodes
            * value -- string representing the primary parameter of the device,
            i.e., resistance of the resistor, or capacitance of the capacitor.
            Stored in string format to support expressions.

    Example:
        a string
            "T1 2 0 3 0 Z0=5 TD=10P"
        is parsed as
            "T", "1", ("2", "0", "3", "0"), {"Z0": 5, "TD": 10e-12}
        The expressions are evaluated using the variables stored in "variables".
        If evaluation fails, the expression is returned as is.
    '''
    inst_name, n1, n2, n3, n4, *expressions = line.split()

    assert(all('=' in item for item in expressions))

    inst_io = (n1, n2, n3, n4)
    _type, label = inst_name[0], inst_name[1:]
    local_params = {}
    for expr in expressions:
        name, val_str = EQ.split(expr)

        success, _, val = parse_expr(val_str, variables)
        if success:
            local_params[name] = val
        else:
            print(f'Failed to propagate parameter {val_str}')
            local_params[name] = val_str

    return _type, label, inst_io, local_params

def parse_source(line: str, variables: dict):
    ''' Parse the current, voltage, or phase source in format
    {type}{label} {node1} {node2} {source_parameters}

    Args:
        line: string containing the device definition. Assumed uppercase and
        without trailing whitespace, as produced by strip_uncomment_upper_join()
        variables: dictionary of parameters visible to the device consiting of
        global parameters and local parameters defined within the subckt
        definition.

    Returns:
        A 4-tuple in format (_type, label, inst_io, value) where
            * _type -- a single character representing type of the device
            * label -- the label of the device
            * inst_io -- a tuple of device nodes
            * value -- string representing the primary parameter of the device,
            i.e., resistance of the resistor, or capacitance of the capacitor.
            Stored in string format to support expressions.

    Example:
        a string
            "I5 0 11 PWL( 0 0 TS 0 2*TS 1000U)"
        is parsed as
            "I", "5", ("0", "11"), "PWL( 0 0 TS 0 2*TS 1000E-6)"
        The expressions are evaluated using the variables stored in "variables".
        If evaluation fails, the expression is returned as is.
    '''
    device_name, n1, n2, param_str = WS.split(line, 3)
    _type = device_name[0]
    label = device_name[1:]
    if param_str.startswith(('PWL', 'PULSE', 'SIN', 'CUS', 'NOISE')):
        bracket_match = BRACKET.search(param_str)
        print(bracket_match)
        src_type = param_str[:bracket_match.start() - 1] # determine the type of the source
        if src_type == 'CUS':
            tmp = WS.split(bracket_match[0].strip())
            values = tmp[:1]
            values_str_list = tmp[1:]
        else:
            values = []
            values_str_list = WS.split(bracket_match[0].strip())
        for item in values_str_list:
            success, _, val = parse_expr(item, variables)
            if success:
                values.append(str(val))
            else:
                print(f'Failed to evaluate parameter {item}')
                values.append(item)
            param_str = f"{src_type}({' '.join(values)})"
    elif param_str.startswith('DC'):
        _, val_str = WS.split(param_str)
        success, _, val = parse_expr(val_str, variables)
        if success:
            param_str = f"DC {val}"
        else:
            print(f'Failed to evaluate parameter {val_str}')
    else:
        success, _, val = parse_expr(param_str, variables)
        if success:
            param_str = f"DC {val}"
        else:
            print(f'Failed to evaluate parameter {param_str}')
        # raise ValueError('Unsupported source type')

    return _type, label, (n1, n2), param_str

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

        for line in f:
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

def parse_initial(file : str):
    ''' Initial parsing of a file to determine local and global parameters and
    existence of the subcircuits.

    Args:
        file: path to a file

    Returns:
        params -- dictionary containing the parameters keyed by the subcircuit
        name. Global parameters are placed in params['main']
    '''
    expressions = {'main': deque()}
    params = {'main': {}}
    line_iter = strip_uncomment_upper_join(file)
    for line in line_iter:
        if line.startswith('.SUBCKT'):
            _, subc_name, *_ = line.split()
            expressions[subc_name] = deque()
            params[subc_name] = {}
            for line in line_iter:
                if line.startswith('.ENDS'):
                    break
                if line.startswith('.PARAM'):
                    _, expr = WS.split(line, 1)
                    expressions[subc_name].appendleft(expr)
        if line.startswith('.PARAM'):
            _, expr = WS.split(line, 1)
            expressions['main'].appendleft(expr)

    for subc, expr_deque in expressions.items():
        params[subc] = {}
        combined_params = params['main'].copy()
        while expr_deque:
            input_string = expr_deque.pop()
            print(input_string)
            success, varname, value = parse_expr(input_string, combined_params)
            if success:
                combined_params[varname] = value
                params[subc][varname] = value
            else:
                # Some dependency exists. Parse the string later
                expr_deque.appendleft(input_string)
    return params

def resolve_params(variables):
    order = [v for v in _topsort(variables.items()) if v in variables]
    for var in order:
        value = variables[var]
        if isinstance(value, str):
            success, varname, value = parse_expr(value, variables)
        variables[var] = value
    return variables

def _topsort(dependency_pairs):
    'Sort values subject to dependency constraints'
    num_heads = defaultdict(int)   # num arrows pointing in
    tails = defaultdict(list)      # list of arrows going out
    heads = []                     # unique list of heads in order first seen
    for h, t in dependency_pairs:
        num_heads[t] += 1
        if h in tails:
            tails[h].append(t)
        else:
            tails[h] = [t]
            heads.append(h)

    ordered = [h for h in heads if h not in num_heads]
    for h in ordered:
        for t in tails[h]:
            num_heads[t] -= 1
            if not num_heads[t]:
                ordered.append(t)
    return ordered

def topsort_subckts(subckts: dict):
    ''' Topological sorting of the subcircuit dependencies.

        Args:
            subckts:  dictionary containing parameters of each known subcircuit

        Returns:
            order -- list containing the names of the subcircuits from the
            bottommost (i.e. containing no dependencies), to the topmost (main)
    '''
    order = []
    dependencies = {sub: d['dependencies'] for sub,d in subckts.items()}
    while dependencies:
        print('Detected following dependencies:')
        print(dependencies)
        for sub, dep in dependencies.items():
            if not dep and sub not in order:
                order.append(sub)
        dependencies = {sub: dep for sub,dep in dependencies.items() if sub not in order}
        for sub, dep in dependencies.items():
            if sub not in order:
                dependencies[sub] = dependencies[sub].difference(order)
    return order

def flatten_netlist(file: str, temp_file: str):
    ''' Parse netlist to determine the hierarchy of the subcircuits. The
    netlist is flattened based on the hierarchical information and is written
    to the target file.

    Args:
        file: path to a file
        temp_file : path to a file storing the flattened netlist

    Returns:
        subckts : dictionary containing subcircuit parameters, including I/O,
        local keyword parameters, devices and subcircuits instantiated within
        the subcircuit, dependencies w.r.t. other subcircuits, and raw lines
        within the subcircuit.
        params : dictionary containing the parameters keyed by the subcircuit
        name. Global parameters are placed in params['main']
        order : list containing the names of the subcircuits from the
        bottommost (i.e. containing no dependencies), to the topmost (main)
        commands : commands within the netlist, such as .print, .plot etc.
    '''

    params = parse_initial(file)

    subckts = {q:{       'io': [],
                     'default_params': {},
                     'devices' : [],
                     'dependencies': set(),
                 'subckt_inst': [],
                        'raw_lines': [] }
                for q in ['main'] + list(params)}
    # Record subckt structures and dependencies
    commands = []
    line_iter = strip_uncomment_upper_join(file)
    for line in line_iter:
        if line.startswith('.SUBCKT'):
            subckts = parse_subckt_def(line, line_iter, subckts)
        elif line.startswith('X'): # Top level subcircuit instance
            subckts['main']['raw_lines'].append(line)
            subc_name = parse_subckt_inst(line, subckts, just_model = True)
            subckts['main']['dependencies'].add(subc_name)
        elif line.startswith(TWONODE_DEVICES + SOURCES + FOURNODE_DEVICES):
            subckts['main']['raw_lines'].append(line)
        elif line.startswith('.') and not line.startswith('.PARAM'):
            commands.append(line)
            print('COMMAND ', line)
        else:
            print('OTHER LINE', line)

    order = topsort_subckts(subckts)
    for subc in order[::-1]: # propagate kw arguments top down
        variables = params['main'] | params[subc]
        assert(not any((arg in subckts[subc]['default_params']) for arg in params[subc]))

        for line in subckts[subc]['raw_lines']:
            print(line)
            if line.startswith(TWONODE_DEVICES):
                _type, label, inst_io, value = parse_twonode(line, variables)
                subckts[subc]['devices'].append((_type, label, inst_io, '', value))
            elif line.startswith(FOURNODE_DEVICES):
                _type, label, inst_io, local_params = parse_fournode(line, variables)
                subckts[subc]['devices'].append((_type, label, inst_io, '', local_params))
            elif line.startswith('X'):
                label, subc_name, inst_io, custom_params = parse_subckt_inst(line, subckts)
                resolved_params = subckts[subc_name]['default_params'] | custom_params
                for param, val_str in resolved_params.items():
                    success, _, val = parse_expr(val_str, variables | params[subc_name])
                    if success:
                        resolved_params[param] = val
                    else:
                        print(f'Failed to propagate parameter {val_str}')
                subckts[subc]['subckt_inst'].append((label, subc_name, inst_io, resolved_params))
            elif line.startswith('B'):
                label, inst_io, model, args = parse_jj(line, variables)
                subckts[subc]['devices'].append(('B', label, inst_io, model, args))
            elif line.startswith(SOURCES):
                _type, label, inst_io, param_str = parse_source(line, variables)
                subckts[subc]['devices'].append((_type, label, inst_io, '', param_str))

    for parent in order: # replace subckt definitions bottom up
        for x_label, child, inst_io, custom_params in subckts[parent]['subckt_inst']:
            variables = resolve_params(subckts[parent]['default_params'] | subckts[child]['default_params'] | params['main'] | params[parent] | params[child] | custom_params)
            model_io = subckts[child]['io']
            io_map = dict(zip(model_io, inst_io)) # map child_io to parent io
            io_map['0'] = '0'
            if x_label == 'I_Q3':
                breakpoint()
            for _type, model_label, model_nodes, model, default in subckts[child]['devices']:
                inst_nodes = [io_map.setdefault(n, f'{n}|{x_label}') for n in model_nodes]
                inst_label = f'{model_label}|{x_label}'
                if isinstance(default, dict):
                    value = default.copy()
                    for name, val_str in value.items():
                        if not isinstance(val_str, str):
                            continue
                        success, _, val = parse_expr(val_str, variables)
                        if success:
                            value[name] = val
                        else:
                            print(f'Failed to propagate parameter {val_str}')
                            value[name] = val_str
                elif isinstance(default, str):
                    success, _, val = parse_expr(default, variables)
                    if success:
                        value = val
                    else:
                        print(f'Failed to propagate parameter {default}')
                        value = default
                elif isinstance(default, float) or isinstance(default, int):
                    value = default
                else:
                    raise TypeError(f'Unsupported parameter {default} type ({type(default)})')
                subckts[parent]['devices'].append((_type, inst_label, inst_nodes, model, value))

    write_flat_file(subckts, commands, temp_file, params)
    return subckts, params, order, commands

def write_flat_file(subckts: dict, commands: list, temp_file: str, params: dict):
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
        for name, value in params['main'].items():
            fobj.write(f'.PARAM {name}={value}\n')

        for _type, label, nodes, model, value in subckts['main']['devices']:
            if isinstance(value, dict):
                value_str = ' '.join(f'{k}={v}' for k,v in value.items())
            else:
                value_str = value
            fobj.write(f'{_type}{label} {" ".join(nodes)} {model} {value_str}\n')
        for line in commands:
            fobj.write(line + '\n')

def run(cir_path, csv_path, temp_file, args):
    subckts, params, order, commands = flatten_netlist(cir_path, temp_file)
    status = call(['josim-cli', '-o', csv_path, temp_file, '-V', '1'])
    if status == 0:
        print('Circuit successfully simulated')
        print(f'The result is written to {csv_path}')
    print(f'Temporary file is saved at {temp_file}')

if __name__ == '__main__':
    parser = ArgumentParser(
        prog = 'NetlistFlattener',
        description = 'The utility to parse netlists while supporting keyword parameters',)

    parser.add_argument('netlist_path')
    parser.add_argument('-t', '--tempfile', default = '')
    parser.add_argument('-d', '--delete', action='store_true')
    parser.add_argument('-n', '--no_sim', action='store_true')

    args, josim_args = parser.parse_known_args()
    print('PARSED:', args.netlist_path, '-t', args.tempfile, '-n', args.no_sim)
    print('SKIPPED:',josim_args)
    folder, cir_name = os.path.split(args.netlist_path)
    base, ext = os.path.splitext(cir_name)
    temp_file = f'{base}_temp.cir' if not args.tempfile else args.tempfile

    subckts, params, order, commands = flatten_netlist(args.netlist_path, temp_file)
    if not args.no_sim:
        print(['josim-cli', temp_file] + josim_args)
        status = call(['josim-cli', temp_file] + josim_args)
        if status == 0:
            print('Circuit successfully simulated')
    if args.delete:
        os.remove(temp_file)
        print(f'Temporary file is deleted')
    else:
        print(f'Temporary file is saved at {temp_file}')
