

gfortran -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -c file_a.F90
# gfortran -fdefault-integer-8 -c file_a.F90
c++ file_b.cpp file_a.o -lgfortran

