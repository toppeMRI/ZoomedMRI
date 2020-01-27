# Design TOPPE sequence for IV imaging

Design small-tip fast recovery (STFR) sequence.

For a description of TOPPE, see https://toppemri.github.io/

# REQUIREMENTS

## Module files that must exist in the (StackOfSpirals) subfolder:

- `tipdown.mod`: 3D tailored inner-volume excitation pulse. See ../RFdesign/tipdown/
- `tipup.mod`: Spectrally-selective tip-up pulse. See ../RFdesign/tipup/

## TOPPE Matlab toolbox

Available at https://github.com/toppeMRI/toppe

```
$ git clone https://github.com/toppeMRI/toppe.git
```

Add the `toppe` folder to your Matlab path.


# USAGE

## Create the TOPPE sequence files

```
>> cd StackOfSpirals;
>> % edit main.m as needed (FOV, matrix size, etc)
>> main;
```

Creates the following files:

- `readout.mod`: Contains cartesian/spiral readout gradients
- `spoiler.mod`: Contains a trapezoidal spoiler gradient
- `modules.txt`: A brief text file listing the various .mod files to be used in the sequence 
- `scanloop.txt`: Each row in this file contains instructions for executing a module once.

## View the .mod files and check PNS

```
>> toppe.plotmod('all');
```

## Preview the sequence in movie (loop) mode

```
>> nModulesPerTR = 4;
>> toppe.playseq(nModulesPerTR);
```

## Execute
