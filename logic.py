from subprocess import call, check_output
from itertools import islice
import numpy as np

def _resolve_args(arg, num_vars):
    # if num_vars is None:
    #     num_vars = int(np.ceil(np.log2(len(arg))))
    if isinstance(arg, (int, np.int, np.int8, np.int8, np.int16, np.int32, np.int64)):
        arg = f'{arg:0{1 << num_vars}b}'
    elif isinstance(arg, str):
        assert(''.join(set(arg)) in ('0','1','01','10'))
        arg = arg.zfill(1 << num_vars)
    elif isinstance(arg, bytes):
        assert(''.join(set(arg)) in (b'0',b'1',b'01',b'10'))
        arg = arg.zfill(1 << num_vars)
    else:
        raise ValueError('Incorrect Argument')
    return arg

def npn_repr(arg, num_vars):
    arg = _resolve_args(arg, num_vars)
    return call(['/Users/brainkz/Documents/GitHub/qca_exploration/npn_repr', arg])

def p_repr(arg, num_vars):
    arg = _resolve_args(arg, num_vars)
    return call(['/Users/brainkz/Documents/GitHub/qca_exploration/p_repr', arg])

def n_repr(arg, num_vars):
    arg = _resolve_args(arg, num_vars)
    return call(['/Users/brainkz/Documents/GitHub/qca_exploration/n_repr', arg])

def is_degenerate(arg, num_vars):
    arg = _resolve_args(arg, num_vars)
    return call(['/Users/brainkz/Documents/GitHub/qca_exploration/is_degenerate', arg])

def tile_check(arg, num_vars):
    arg = _resolve_args(arg, num_vars)
    out = check_output(['/Users/brainkz/Documents/GitHub/qca_exploration/tile_check_3', arg])
    # print(out)
    is_degen = int(out[ 0: 1],2)
    npn_repr = int(out[ 1: 9],2)
    n_repr   = int(out[ 9:17],2)
    p_repr   = int(out[17:25],2)
    # print(arg)
    # print(out)
    # print(bool(is_degen), npn_repr, n_repr, p_repr)
    return bool(is_degen), npn_repr, n_repr, p_repr

def batched(iterable, n):
    "Batch data into tuples of length n. The last batch may be shorter."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while (batch := tuple(islice(it, n))):
        yield batch

def batch_check(out_vectors, tt_check, num_vars, chunk_size = 1000):
    arguments = {}
    for arg in out_vectors:
        if arg not in tt_check:
            arguments[arg] = _resolve_args(arg, num_vars)
    if not arguments:
        return None
    outs = []
    for args in batched(arguments.values(), chunk_size):
        out = check_output([f'/Users/brainkz/Documents/GitHub/qca_exploration/batch_check_{num_vars}', *args]).decode().strip().split('\n')
        outs.extend(out)
    for out, vec in zip(outs, arguments):
        it = iter(out)
        is_degen = int(next(it), 2)
        args = [it] * (1 << num_vars)
        npn_class, n_class, p_class = list(map(lambda q:int(''.join(q),2), zip(*args)))
        # npn_class_str, n_class_str, p_class_str = list(map(b''.join, zip(*args)))
        # npn_class = int(out[ 1: 9],2)
        # n_class   = int(out[ 9:17],2)
        # p_class   = int(out[17:25],2)
        tt_check[vec] = is_degen, npn_class, n_class, p_class
    return None

def np_check(vectors, tt_check, num_vars):
    arguments = {}
    for arg in vectors:
        if arg not in tt_check:
            arguments[arg] = _resolve_args(arg, num_vars)
    if not arguments:
        return None
    outs = check_output([f'/Users/brainkz/Documents/GitHub/kitty/build/examples/batch_check_np_{num_vars}', *arguments.values()]).decode()
    for out, vec in zip(outs.split('\n'), arguments):
        is_degen = out[0]
        np_class = out[1:]
        tt_check[vec] = (is_degen == '1'), int(np_class,2)
    return None

if __name__ == '__main__':
    nvars = 3
    vectors = list(range(2**2**nvars))
    tt_check = {}
    np_check(vectors, tt_check, nvars)
    

# if __name__ == '__main__':
#     from json import dump, load
#     with open('/Users/brainkz/Documents/GitHub/qca_exploration/npn_classes.json') as f:
#         npn_classes = load(f)
    
#     filtered_classes = {}
#     for nvars, classes in npn_classes.items():
#         nvars = int(nvars)
#         filtered_classes[nvars] = {}
#         for repr, family in classes.items():
#             repr = int(repr)
#             arg = _resolve_args(repr, nvars)
#             is_degen = call(['/Users/brainkz/Documents/GitHub/qca_exploration/is_degenerate', arg])
#             if not is_degen:
#                 filtered_classes[nvars][repr] = family
#             else:
#                 print(f'Degenerate class found: {repr}')
#     with open('/Users/brainkz/Documents/GitHub/qca_exploration/npn_classes_filtered.json', 'w') as f:
#         dump(filtered_classes, f)             
                
                