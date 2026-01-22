# Makefile for Parallel M Solver

FC = mpif90
FFLAGS = -O3 -march=native
LDFLAGS = 

# Main targets
all: parallel_block_solver lu_solver

parallel_block_solver: parallel_block_solver.f90
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

lu_solver: lu_solver.f90
	gfortran $(FFLAGS) -o $@ $<

block_schur_solver: block_schur_solver.f90
	gfortran $(FFLAGS) -o $@ $<

# Utility targets
analyze_structure: analyze_structure.f90
	gfortran -O2 -o $@ $<

debug_solver: debug_solver.f90
	gfortran -O2 -o $@ $<

clean:
	rm -f parallel_block_solver lu_solver block_schur_solver *.o *.mod
	rm -f analyze_structure debug_solver check_upper check_coupling
	rm -f solution_*.txt

run: parallel_block_solver
	mpirun -np 4 ./parallel_block_solver

run_serial: lu_solver
	./lu_solver

.PHONY: all clean run run_serial
