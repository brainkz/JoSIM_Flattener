__Netlist flattener__

Netlist flattening utility converting the file with nested subcircuits into an equivalent netlist without subcircuits.

This Python tool enables the use of named subcircuit parameters in [JoSIM](https://github.com/JoeyDelp/JoSIM), Superconductor Circuit Simulator.

The tool replaces a netlist with a parameterized subcircuit

```hspice
* Input netlist: amp_test.cir
X1 TEST 1 0 AMP=1.0 FREQ=1e9
I1      0 1 1.0

.subckt TEST IN GND AMP=1.0 FREQ=1e9
    R1 IN GND 10.0
    I1 GND IN sin(0 AMP FREQ)
.ends

.tran 0.25p 5n 0 1.0p
.print v(1)
.end
```

into a flattened netlist with parameters replaced by the specified values

```hspice
* Flattened netlist: amp_test_temp.cir
I1 0 1  DC 1.0
R1|1 1 0  10.0
I1|1 0 1  SIN(0 1.0 1000000000.0)
.TRAN 0.25P 5N 0 1.0P
.PRINT V(1)
.END
```

By default, JoSIM is called next to simulate the flattened netlist.

__Usage__:

```
python sim_flatten.py path-to-netlist.cir [-t (--tempfile) path-to-temporary-file.cir] [-n (--no_sim)] [-d (--delete)]
```
where 
- `path-to-netlist.cir` is the path to netlist file to be flattened
- `-t` denotes location of temporary file. The default is `{path-to-netlist}_temp.cir`;
- `-n` asks the tool to skip simulation and only to generate the temporary file
- `-d` asks the tool to delete the temporary file at the end of the script.

__Dependencies__:
- python3 
- pyparsing
- josim-cli

