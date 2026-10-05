#!/bin/bash
# Build the QSS fork with the GCC matching the Python toolchain
module load GCC/13.3.0 >/dev/null 2>&1
cd /work/mehlc/0_CODES_GITS/cantera-qss
exec /work/mehlc/0_CODES_GITS/MechanismReduction/venv/bin/scons build -j16 "$@"
