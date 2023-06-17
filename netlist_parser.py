
import os
import re
import sys
import sympy 
import networkx as nx

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
        
class Subckt:
    # def __init__(self, name, terminals, parameters, lines, children):
    #     self.name = name
    #     self.terminals = terminals
    #     self.locals = parameters
    #     self.lines = lines
    #     self.children = children
    
    def __init__(self, name, io_defn, var_params = {}, const_params = {}, raw_lines = [], children = []):
        self.name = name
        self.io = io_defn
        self.var_params = var_params
        self.const_params = const_params
        self.raw_lines = raw_lines
        
    def __str__(self):
        return f"Name: {self.name}\nTerminals: {self.terminals}\nParameters: {self.parameters}\nDefinition: {self.definition}"
    
class Inst:
    def __init__(self, inst_name: str, inst_io: list, subc_name: str, var_params: dict = {} ):
        self.name = inst_name
        self.io = inst_io
        self.subc_name = subc_name
        self.var_params = var_params
    
    

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

def to_sympy(expr_str, subs = {}):
    while item := SUFFIXED_NUM.search(expr_str):
        suf = ALPHA.search(item[0])[0]
        sub_str = item[0].replace(suf, SUFFIX[suf])
        expr_str = expr_str[:item.start()] + sub_str + expr_str[item.end():]
        print(expr_str)
    expr = sympy.parse_expr(expr_str)
    return expr.evalf(subs=subs)
    
def expr_dict(tokens : list, subs: dict = {}):
    expressions = {}
    for expr in tokens:
        key, val = EQ.split(expr)
        expressions[key] = to_sympy(val, subs)
    return expressions
        
def elem_defn(line: str, subs: dict = {}):
    if line[0] in TWONODE_DEVICES:
        tokens = WS.split(line)
        return tokens[0][0], tokens[0][1:], tokens[1:3], expr_dict(tokens[3:], subs)
    elif line[0] in FOURNODE_DEVICES:
        tokens = WS.split(line)
        return tokens[0][0], tokens[0][1:], tokens[1:5], expr_dict(tokens[5:], subs)
    elif line[0] == 'X':
        tokens = WS.split(line)
        for i, t in enumerate(tokens):
            if '=' in t:
                return tokens[0][0], tokens[0][1:], tokens[1:i-1], tokens[i-1], expr_dict(tokens[i:], subs)
        else:
            return tokens[0][0], tokens[0][1:], tokens[1:-1], tokens[-1], {}

def subc_first_line(line):
    _, name, *items = WS.split(line)
    io_defn = []
    local_params = {}
    for item in items:
        if '=' not in item:
            io_defn.append(item)
        else:
            key, val = EQ.split(item)
            local_params[key] = to_sympy(val)
    return name, io_defn, local_params

def create_subckt_sublists(lines) -> tuple[list[list[str]], list[str]]:
    main_lines = []
    subckt_sublists = []
    subckt_lines = []
    is_subckt = False

    for line in lines:
        line = line.strip()

        if line.startswith(".SUBCKT"):
            is_subckt = True
            subckt_lines.append(line)
        elif line.startswith(".ENDS"):
            if is_subckt:
                subckt_lines.append(line)
                subckt_sublists.append(subckt_lines)
                subckt_lines = []
                is_subckt = False
        elif is_subckt:
            subckt_lines.append(line)
        else:
            main_lines.append(line)
    return subckt_sublists, main_lines

def parse_twonode(line: str):
    device_name, *inst_io, val_str = line.split()
    _type, label = device_name[0], device_name[1:]
    return _type, label, inst_io, val_str

def parse_fournode(line: str):
    device_name, n1, n2, n3, n4, *expressions = line.split()
    assert(all('=' in item for item in expressions))
    _type, label = device_name[0], device_name[1:]
    return _type, label, (n1, n2, n3, n4), expr2dict(expressions)

def sub_instances(prefix, input_params, lines, subckts):
    new_lines = []
    params = input_params.copy()
    for line in lines:
        if line.startswith('.PARAM'):
            _, expr = WS.split(line)
            key, val = EQ.split(expr)
            val = to_sympy(val, subs=params)
            params[key] = val
    for line in lines:
        breakpoint()
        if line.startswith(".PARAM"):
            continue
        
        if line.startswith("X"):
            _, inst_name, inst_io, subc_name, inst_params = elem_defn(line, params)
        if line[0] in TWONODE_DEVICES:
            _type, label, inst_io, val_str = parse_twonode(line)
            new_io
            new_line = f"{_type}_{prefix}_{label} "
            
            

def subc_replace(lines : list, subckts : dict):
    subckts = {}
    # tree T tracks subckt instances
    # Nodes represent subckt instances
    # Edges contain parameters passed to the instance
    T = nx.DiGraph() 
    main = Inst('main_inst', [], 'main_subc', {})
    subckt_line_lists, main_lines = create_subckt_sublists(lines)
    
    main_params = CONSTANTS.copy()
    sub_instances('main_inst', CONSTANTS, main_lines, subckts)
    
    for subckt_lines in subckt_line_lists:
        subc_name, io_defn, local_params = subc_first_line(subckt_lines[0])
        # parse_instances(T, 'main_inst', CONSTANTS, main_lines)
    
    return T, main_params
    
    # T.add_node('main_inst') #intentionally lowercase, in case a subckt 'MAIN' is defined
    # T.nodes['main_inst']['lines'] = main_lines
    # T.nodes['main_inst']['params'] = main_params

def extract_subckts(lines):
    
    main = Subckt('main', [], [], const_params={}, raw_lines=[])
    subckts = {}
    DAG = nx.MultiDiGraph()
    DAG.add_node('main') #intentionally lowercase, in case a subcircuit 'MAIN' is defined
    
    is_subcircuit = False

    for i,line in enumerate(lines):
        line = line.strip()
        if line.startswith(".SUBCKT"):
            print(i,line)
            assert (not is_subcircuit), "Nested subcircuit definitions are not supported"
            is_subcircuit = True
            name, io_defn, local_params = subc_first_line(line)
            subc = Subckt(name, io_defn, local_params, const_params={}, raw_lines=[])
        elif line.startswith(".ENDS"):
            print(i,line)
            assert (is_subcircuit), f".ENDS statement outside subcircuit"
            is_subcircuit = False
            subckts[subc.name] = subc
            DAG.add_node(subc.name)
        elif line.startswith(".PARAM"):
            _, expr = WS.split(line)
            key, val = EQ.split(expr)
            if is_subcircuit:
                subc.raw_lines.append(line)
                subc.const_params[key] = to_sympy(val)
            else:
                main.const_params[key] = to_sympy(val)
                
        elif line.startswith("X"):
            _, inst_name, inst_io, subc_name, inst_params = elem_defn(line)
            parent = subc.name if is_subcircuit else 'main'
            DAG.add_edge(parent, subc_name,
                         **{
                             'inst_name':inst_name,
                             'inst_io':inst_io,
                             'inst_params':inst_params
                         })
            if is_subcircuit:
                subc.raw_lines.append(line)
            else:
                main.raw_lines.append(line)
        elif is_subcircuit:
            subc.raw_lines.append(line)
        else:
            main.raw_lines.append(line)
    
    subckts['main'] = main
    print(DAG.edges)
    live_subckts = nx.descendants(DAG, 'main')
    live_subckts.add('main')
    DAG.remove_nodes_from([n for n in DAG if n not in live_subckts])
    subckts = {k:v for k,v in subckts.items() if k in live_subckts}
    return DAG, subckts

def topological_edges(DAG, reverse = False):
    order = list(nx.topological_sort(DAG))
    if reverse:
        order = order[::-1]
    for parent in order:
        for child in DAG[parent]:
            yield parent, child

def propagate_params(DAG, subckts):
    
    it = nx.topological_sort(DAG)
    main_n = next(it)
    main_subc = subckts[main_n]
    DAG.nodes[main_n]['params'] = CONSTANTS | main_subc.const_params
    for n in it:
        subc = subckts[n]
        params = {}
        for u in DAG.pred[n]: 
            params |= DAG.nodes[u]['params'] #gather predecessor params
            # Note: conflicts resolution is not controlled
        params |= subc.var_params
        for e in DAG.in_edges(n): 
            params |= DAG.edges[u,n]['inst_params'] #overwrite with instance params
        params |= subc.const_params #overwrite with const params (defined with .param statements within .subckt body)
        DAG.nodes[n]['params'] = params
    
    for n in reversed(nx.topological_sort(DAG)):
        breakpoint()
        

            
    
    

if __name__ == '__main__':
    FILE = '/Users/brainkz/Documents/GitHub/SFQ/pulse_code.cir'
    lines = list(strip_uncomment_upper_join(FILE))
    # DAG, subckts = extract_subckts(lines)
    # propagate_params(DAG, subckts)
    DAG, subckts = extract_subckts(lines)
    T, main_params = subc_replace(lines, subckts)
    