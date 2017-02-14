# makefile of kondo effect

##################################################################################

clean:
	rm -r *.out *.o *.mod

#################################################
# Fortran
#
F90= /opt/pgi/linux86-64/15.9/mpi/mpich/bin/mpif90

F90FLAGS=  -tp istanbul-64 -w -g -O0 -mcmodel=medium
################################################

# INCLUDES
# #MPIINCS = -I/pkg/mpi/pgi/11.6/mvapich-1.2/include
MPIINCS = -I/opt/pgi/15.9/linux86-64/mpi/mpich/include
# #MPIINCS = 
LAPINCS =
FFTINCS = -I/home/CHEN/bch/tools/64fftw/include
CCINCS  = $(MPIINCS) $(FFTINCS)
MISCINCS =
INCS = $(MISCINCS) $(MPIINCS) $(LAPINCS) $(FFTINCS)
#
#
#

#CURRENT part

shot_bch1: shot_bch1.f90 GaussLeg.f90
	$(F90) -o shot_bch1 shot_bch1.f90 GaussLeg.f90 glob_var.f90 


shot_bch: shot_bch.f90 GaussLeg.f90
	$(F90) -o shot_bch shot_bch.f90 GaussLeg.f90 glob_var.f90 

shot_bch2k: shot_bch2k.f90 GaussLeg.f90
	$(F90) -o shot_bch2k shot_bch2k.f90 GaussLeg.f90 glob_var.f90 


shot: shot.f90 GaussLeg.f90
	$(F90) -o shot shot.f90 GaussLeg.f90 glob_var.f90


shot1: shot1.f90 GaussLeg.f90
	$(F90) -o shot1 shot1.f90 GaussLeg.f90 glob_var.f90
.f90.o:
	$(F90) $(F90FLAGS) $(INCS) -c $<
