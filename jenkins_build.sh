#!/usr/bin/env bash

# Source common shell script variables
source set_vars.sh

# Source the Intel compilers
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

# Make bash script more like high-level languages.
# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
# Have to do this after sourcing ifortvars.sh becuase the shell script has unbound variables
set -Eeuxo pipefail

# Clone dependencies
cd ..
for deprepodir in "${!deprepo[@]}"; do
    if [ ! -d ${deprepodir} ]; then
        all_proxy=${proxyout} git clone ${deprepo[$deprepodir]}
    else
        cd ${deprepodir} 
        if [ ${deprepodir} == "eigen" ]; then
            all_proxy=${proxyout} git checkout master
        else
            all_proxy=${proxyout} git checkout dev
        fi
        all_proxy=${proxyout} git pull --ff-only
        cd ..
    fi
done
# Perform repo tests
cd ${workdir}/src/cpp/tests/${repo}/
if [ -f ${repo}.o ]; then
    make clean
fi
make
