import numpy as np
from operator import and_
from operator import or_
from itertools import product, combinations, chain
from functools import reduce

inv = lambda a, bitlen: (((1<<bitlen)-1)^a)

def pset(iterable, minlen = 1, maxlen = None):
    """Yields all powersets of the iterable of size between minlen and maxlen.  """
    s = list(iterable)
    if maxlen is None:
        maxlen = len(s)
    return chain.from_iterable(combinations(s, r) for r in range(minlen, maxlen + 1))

def bin_codes(nvars, with_neg = True):
    ''' Generate binary codes for {nvars} variables.
    By default, also produce the negated codes.
    Examples: 
    >> codes = bin_codes(2, with_neg = True ) -> [0b0101, 0b1010, 0b0011, 0b1100]
    >> codes = bin_codes(2, with_neg = False) -> [0b01010101, 0b00110011, 0b00001111]
    '''
    bitlen = 1 << nvars
    codes = []
    for i in range(nvars):
        base = (1 << (1 << i)) - 1
        deg = np.arange(0, bitlen, 1 << (i + 1), dtype=int)
        num = np.sum(base * (1 << deg))
        # print(f"{num:0{bitlen}b} -> {num}")
        codes.append(num)
        if with_neg:
            codes.append(inv(num, bitlen))
    return codes

def or_up_to_depth(codes, maxdepth, nvars):
    '''Returns all possible truth tables obtainable from given codes using an OR tree with maximum depth of {maxdepth}'''
    
    for i in range(1, maxdepth+1):
        ncodes = len(codes)
        print(f'\tIteration {i}')
        print('\tGenerating code array')
        code_array = np.bitwise_or.outer(codes, codes)[np.tri(ncodes,dtype=bool)]
        print(f'\tCode array size is {len(code_array)}')
        print('\tFinding uniques')
        new_codes = np.unique(code_array)
        print(len(new_codes))
        if len(new_codes) == len(codes) or len(new_codes) == 1<<(1<<nvars):
            print(f'Aborted at depth {i}')
            return new_codes
        else:
            codes = new_codes
    return codes
nbits = 4
bitlen = 1 << nbits
mask = (1 << bitlen) - 1
codes = bin_codes(nbits, with_neg = False)

print(sorted(codes))
depth = 5
outputs = set()
for i in range(1, depth + 1):
    print(f'Depth = {i}')
    inputs = or_up_to_depth(codes, maxdepth = 5, nvars = nbits)
    if len(inputs) == mask+1:
        print(f'All possible functions are achievable with at most {i} cycles')
        break
    idx = np.tri(len(inputs),len(inputs),-1,dtype=bool)
    idxm1 = np.tri(len(inputs),len(inputs),-1,dtype=bool)
    ORed  = np.bitwise_or.outer(inputs, inputs)
    ANDed = np.bitwise_and.outer(inputs, inputs)
    # IMPed = np.bitwise_and.outer(inputs, np.bitwise_xor(inputs, mask))
    NORed = np.bitwise_xor(ORed[idx], mask)
    # XORed = np.bitwise_xor.outer(inputs, inputs)
    # IMPed, XORed,
    outputs.update(n for q in chain(ORed[idx], ANDed[idxm1], NORed,) for n in q.ravel())

    remaining = mask+1 - len(outputs)
    print(f'Remaining: {remaining}')
    if not remaining:
        print(f'All possible functions are achievable with at most {i} cycles')
        break
    # for num in set(range(1 << bitlen)).difference(sorted(outputs2)):
    #     print(f"{num:0{bitlen}b} -> {num}")
    codes = list(set().union(codes, outputs))
