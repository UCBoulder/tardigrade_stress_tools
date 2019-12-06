Tools for computing stress-strain behaviors. The stress functions should be 
able to return the stress and the jacobian w.r.t. the strain metric of 
interest.

Note: In order to use the Intel compiler one must run the following command 
in a bash prompt:
source /apps/intel2016/bin/ifortvars.sh -arch intel64 -platform linux

This is the same command that the abaqus command issues. It may be that 
this command will change on different platforms.

---

---

Dependencies: 

These tools have several dependencies that must be available in the same parent
directory as this repo. 

* eigen: https://gitlab.com/libeigen/eigen
* constitutive\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/constitutive_tools
* error\_tools: https://xcp-stash.lanl.gov/projects/MM/repos/error_tools
