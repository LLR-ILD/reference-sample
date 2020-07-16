# Reference Sample

TODO: Proper Readme! So far just copy-pasted.

## Steering

Steering functionality was outsourced into its own project, [`pySteer`](https://github.com/kunathj/pySteer).

## Processors

Some key phrases can tell what the processor is written for:

- `Check*`: Indicates that this processor/folder is not needed in a final
  analysis. Instead, it is a processor that shall help understand a data file or
  (3rd party) processor better.
- `DraftProcessor`: Can be used to try out ideas without the necessity to write
  a full new processor.
- `ExampleProcessor`: Provided as a copy-paste skeleton for new processors.
- `*Util`: The class is not a processor itself, but provides functionality for
  use in (multiple) other processors. The `_util.h/.cc` files should be in a
  folder  called `shared_functionality`.

## Installation

Explain here:

- what are the package dependencies (iLCSoft, others ?)
- how to compile your package. Should normally be something like:

```shell
source /path/to/ilcsoft/init_ilcsoft.sh
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
```

## How to run the analysis

Explain here:

- where to find data needed for your analysis or how to produce them
- how to run you analysis:
  - Marlin processors to run ?
  - ROOT macros to run ?
  - Shell scripts ?
  - Run the analysis on grid if you provide scripts for that

Example:

```shell
export MARLIN_DLL=./lib/libZH_ILD_ReferenceSample.so
Marlin ./scripts/ExampleProcessor.xml
```

If you want to provide a lot of details on your analysis, use the doc/Readme.md and point to it from this Readme.md file:

More documentation available here in [doc/Readme.md](doc/Readme.md) !

## Issues and contact

Explain here how can people reach you:

- via the Github issue interface. For the skeleton package: https://github.com/ILDAnaSoft/ZH_ILD_ReferenceSample/issues
- **not mandatory**:
  - email address
  - working institute
