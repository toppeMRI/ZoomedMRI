# Design TOPPE sequence for IV imaging

Design small-tip fast recovery (STFR) sequence.

For a description of TOPPE, see https://toppemri.github.io/

# Requirements

## Module files that must reside in this folder:

- `tipdown.mod`: 3D tailored inner-volume excitation pulse. See ../RFdesign/tipdown/
- `tipup.mod`: Spectrally-selective tip-up pulse. See ../RFdesign/tipup/

## TOPPE Matlab toolbox

Available at https://github.com/toppeMRI/toppe

```
$ git clone https://github.com/toppeMRI/toppe.git
```

Add the `toppe` folder to your Matlab path.


# Usage

```
>> main;
```

Creates the following files:

- `readout.mod`
- `spoiler.mod`
- `modules.txt`

Under construction

