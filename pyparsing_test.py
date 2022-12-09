# Simple expression parse base on SimpleCalc.py from the
# pyparsing module.
# Supports out-of-order definitions.

from pyparsing import ParseException, Word, alphas, alphanums
from collections import deque
from fourFn import BNF, exprStack, evaluate_stack


def parse_line(input_string, variables, debug_flag = False):
    # Reset to an empty exprStack
    del exprStack[:]

    arithExpr = BNF()
    ident = Word(alphas, alphanums).setName("identifier")
    assignment = ident("varname") + "=" + arithExpr
    pattern = assignment | arithExpr

    if input_string:
        # try parsing the input string
        try:
            L = pattern.parseString(input_string, parseAll=True)
        except ParseException as err:
            L = ["Parse Failure", input_string, (str(err), err.line, err.column)]

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

if __name__ == "__main__":
    debug_flag = False
    lines = deque([
    'Phi0=2.067833848E-15',
    'B0=1',
    'Ic0=0.0001',
    'IcRs=100e-6*6.859904418',
    'B0Rs=IcRs/Ic0*B0',
    'Rsheet=2',
    'Lsheet=1.13e-12',
    'LP=0.2e-12',
    'IC=2.5',
    'LB=2e-12',
    'IB1=BiasCoef*Ic0*IC',
    'BiasCoef=0.70',])

    vars = {}
    while lines:
        input_string = lines.pop().strip().upper()
        status, varname, value = parse_line(input_string, vars)
        if status:
            vars[varname] = value
        else:
            # Parse the string later
            lines.appendleft(input_string)
