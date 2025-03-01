#!/bin/csh
#
# run a simple test receiver function example
#
#  note that the iterative code is set up to
#    produce a result that when convolved with
#    the vertical, reproduces the radial. This
#    is not the case in the FDomain rftn. The
#    averaging function of the FDomain result
#    should have unit area, but it does not.
#    if you normalize the FDomain rftn, lac.eqr,
#    by the area in the rftn, you will get the
#    same amplitudes as the time domain.
#    See the macro compare.macro for details
#
../iterdecon <<end
1992.156.04.04.03.EDM_sp.r 
1992.156.04.04.03.EDM_sp.z 
400      * nbumps
5.0      * phase delay for result
0.001    * min error improvement to accept
2.5      * Gaussian width factor
1        * 1 allows negative bumps
0        * 0 form minimal output (1) will output lots of files
end
#
# rename the result (decon.out) to a more useful value
#
mv decon.out 1992.156.04.04.03.EDM.i.eqr
