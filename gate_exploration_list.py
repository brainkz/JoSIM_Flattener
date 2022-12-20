import numpy as np
from operator import and_
from operator import or_
from itertools import product, combinations, chain
from functools import reduce
import numba
from logic import batch_check, p_repr

inv = lambda a, bitlen: (((1<<bitlen)-1)^a)

def pset(iterable, minlen = 1, maxlen = None):
    """Yields all powersets of the iterable of size between minlen and maxlen.  """
    s = list(iterable)
    if maxlen is None:
        maxlen = len(s)
    return chain.from_iterable(combinations(s, r) for r in range(minlen, maxlen + 1))

def bin_codes(nvars, with_neg = True, add_const = True):
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
    if add_const:
        codes.extend((0, (1<<bitlen)-1))
    return np.sort(codes)


def or_up_to_depth(codes, depths, maxdepth, seq_depth, nvars):
    inputs = depths
    for i in range(1, maxdepth+1):
        ncodes = len(codes)
        idx = np.tri(ncodes, ncodes, -1, dtype=bool)
        print(f'\tIteration {i}:\n\t\tGenerating code array')
        code_array = np.bitwise_or.outer(codes, codes)[idx]
        inputs[code_array] = seq_depth

        print(f'\t\tCode array size is {len(code_array)}')
        if np.all(inputs < 255):
            print(f'Aborted at depth {i}')
            return np.where(inputs)[0]
        else:
            codes = np.where(inputs)[0]
    return codes

'''
@numba.njit
def outer_or(arr):
    n = len(arr)
    res = np.zeros((n,n), arr.dtype) + 255
    for ii in range(n):
        for jj in range(ii+1):
            res[ii,jj] = arr[ii] | arr[jj]
    return res

@numba.njit
def outer_add_plus(arr):
    n = len(arr)
    res = np.zeros((n,n), arr.dtype) + 255
    for ii in range(n):
        for jj in range(ii+1):
            res[ii,jj] = arr[ii] + arr[jj] + 1
    return res
'''
from numba import int_

@numba.njit
def comb_or(codes, comb_costs):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        for j in range(i+1, ncodes):
            cj = codes[j]
            new_code = ci | cj
            new_cost = comb_costs[ci] + comb_costs[cj] + 1
            if new_cost < comb_costs[new_code]:
                comb_costs[new_code] = new_cost
    return comb_costs

def combinational(codes, ntt,):
    comb_depths = np.zeros(ntt, dtype=np.uint8) - 1
    comb_depths[codes] = 0
    comb_costs = np.zeros(ntt, dtype=np.uint8) - 1
    comb_costs[codes] = 0
    ncodes = len(codes)
    for i in range(1, max_comb_depth+1):
        mask = np.tri(ncodes, ncodes, -1, dtype=bool)
        print(f'\tIteration {i}:\n\t\t#inputs={len(codes)}\n\t\tGenerating code array')
        
        comb_costs = comb_or(codes, comb_costs)
        discovered_codes = np.where(comb_costs < 255)
        new_codes = np.setdiff1d(discovered_codes, codes)
        if len(new_codes) == 0:
            print(f'\t\tNo new functions discovered at depth {i}')
            break
        else:
            print(f'\t\tDiscovered {len(new_codes)} new functions at depth {i}')
            
        
        comb_depths[new_codes] = i
        codes = np.where(comb_depths < 255)[0]
        ncodes = len(codes)
        if ncodes == ntt: # all possible functions discovered
            # No need to update depths or codes
            print(f'\t\tAll functions discovered at depth {i}')
            break
        
    return codes, comb_depths, comb_costs

@numba.njit
def seq_or(codes, seq_depths, seq_depth):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        for j in range(i+1, ncodes):
            cj = codes[j]
            new_code = ci | cj
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
    return seq_depths

@numba.njit
def seq_and(codes, seq_depths, seq_depth):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        for j in range(i+1, ncodes):
            cj = codes[j]
            new_code = ci & cj
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
    return seq_depths

@numba.njit
def seq_nor(codes, seq_depths, seq_depth, ones):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        for j in range(i+1, ncodes):
            cj = codes[j]
            new_code = (ci | cj) ^ ones
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
    return seq_depths


# @numba.njit
def seq_self(codes, seq_depths, seq_depth):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        if seq_depths[ci] == 255:
            seq_depths[ci] = seq_depth
    return seq_depths

@numba.njit
def seq_all(codes, seq_depths, seq_depth, ones):
    ncodes = len(codes)
    for i in range(ncodes):
        ci = codes[i]
        for j in range(i+1, ncodes):
            cj = codes[j]
            new_code = ci | cj
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
            new_code = new_code ^ ones
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
            new_code = ci & cj
            if seq_depths[new_code] > seq_depth:
                seq_depths[new_code] = seq_depth
    return seq_depths

if __name__ == '__main__':
    from time import perf_counter
    max_comb_depth = 5
    max_seq_depth = 5
    
    nvars = 4
    bitlen = 1 << nvars
    ntt = (1 << bitlen)
    ones = ntt - 1
    
    codes = bin_codes(nvars, with_neg = True, add_const = False)
    seq_depths = np.zeros(ntt, dtype=np.uint8) - 1
    seq_depths[codes] = 0
    print(codes)
    
    
    for seq_depth in range(max_seq_depth):
        print(f'ANALYZING {seq_depth}')
        discovered_codes, comb_depths, comb_costs = combinational(codes, ntt,)
        seq_depths = seq_self(discovered_codes, seq_depths, seq_depth)
        
        assert(np.setdiff1d(discovered_codes, np.where(seq_depths < 255)[0]).size == 0)
        if len(discovered_codes) == ntt:
            print(f'All possible {nvars}-input functions are achievable with at most {seq_depth} cycles')
            break
        else:
            seq_depths = seq_all(discovered_codes, seq_depths, seq_depth, ones)
                
            codes = np.where(seq_depths < 255)[0]
            if len(codes) == ntt:
                print(f'All possible {nvars}-input functions are achievable with at most {seq_depth} cycles')
                break    
    
    tt_check = {}
    batch_check(list(range(1<<(1<<nvars))), tt_check, nvars)
    p_classes = {k:v[3] for k,v in tt_check.items()}
    npn_classes = {k:v[1] for k,v in tt_check.items()}
    onecycle_p = {p_classes[q] for q in np.where(seq_depths==0)[0]}
    onecycle_npn = {npn_classes[q] for q in np.where(seq_depths==0)[0]}
    
    
    
    
    
    # start = perf_counter()
    # elapsed = perf_counter() - start
    # print(elapsed)
if False:
    
    max_seq_depth = 5

    outputs = set()
    for i in range(1, max_seq_depth + 1):
        print(f'Depth = {i}')
        inputs = or_up_to_depth(codes, maxdepth = 5, nvars = nvars)
        if len(inputs) == ones+1:
            print(f'All possible functions are achievable with at most {i} cycles')
            break
        idx = np.tri(len(inputs),len(inputs),-1,dtype=bool)
        idxm1 = np.tri(len(inputs),len(inputs),-1,dtype=bool)
        ORed  = np.bitwise_or.outer(inputs, inputs)
        ANDed = np.bitwise_and.outer(inputs, inputs)
        # IMPed = np.bitwise_and.outer(inputs, np.bitwise_xor(inputs, mask))
        NORed = np.bitwise_xor(ORed[idx], ones)
        # XORed = np.bitwise_xor.outer(inputs, inputs)
        # IMPed, XORed,
        outputs.update(n for q in chain(ORed[idx], ANDed[idxm1], NORed,) for n in q.ravel())

        remaining = ones+1 - len(outputs)
        print(f'Remaining: {remaining}')
        if not remaining:
            print(f'All possible functions are achievable with at most {i} cycles')
            break
        # for num in set(range(1 << bitlen)).difference(sorted(outputs2)):
        #     print(f"{num:0{bitlen}b} -> {num}")
        codes = list(set().union(codes, outputs))
