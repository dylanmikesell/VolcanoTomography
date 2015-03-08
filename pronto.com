# Unix execution script for program PRONTO
#
# Primary function of this script is to associate Fortran source code
# numeric designators of input and output files with file pathnames.
# Pathnames can be relative or absolute (i.e., relative to root directory).
# PRONTO execution command is also included (preceded by a UNIX run time
# command).

rm -f fort.*
ln -s pronto.run      fort.10
ln -s data.obs        fort.11
ln -s vel.initial     fort.12
ln -s vel.ref         fort.13
ln -s topo.surf       fort.14
ln -s vel.final       fort.15
ln -s ray.final       fort.16
ln -s data.prd        fort.17
time pronto.exe
rm -f fort.*

# Brief description of input/output files follows.  Each file has a
# numeric designation (fort.#) (used in Fortran source code READ or 
# WRITE statements) and a suggested pathname (relative to current
# directory above).
#
# Files input to PRONTO:
#    1) fort.10: pronto.run:  algorithm execution control parameters
#    2) fort.11: data.obs:    observed data (observed traveltimes, source
#                                and receiver position coordinates, data
#                                weights)
#    3) fort.12: vel.initial: velocity model used to initiate tomographic
#                                inversion (optional)
#    4) fort.13: vel.ref:     reference velocity model, associated with
#                                regularization constraints imposed during
#                                tomographic inversion (optional)
#    5) fort.14: topo.surf:   surface topography model (optional)
# 
# Files output from PRONTO:
#    1) fort.15: vel.final:   tomographically determined velocity model
#    2) fort.16: ray.final:   raypath density distribution for final
#                                velocity model
#    3) fort.17: data.prd:    predicted traveltime data (or traveltime
#                                residuals) for final velocity model
#
