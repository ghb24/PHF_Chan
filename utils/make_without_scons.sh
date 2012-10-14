#!/bin/bash
# -lboost_program_options might need a package to be installed (part of
#   boost C++ libraries).
export LIBS_BASE="-lboost_program_options "

# note on math libraries:
#  o program assumes by default that a BLAS with 64bit integers is linked.
#    for MKL, that means using the mkl_intel_ilp64 variant (NOT the lp64 one!),
#    for ACML that means using the version from the _int64 subdirectory
#    (e.g., gfortran64_fma4_int64/)
#  o If a 64bit integer blas is not available, you can re-define "FORTINT"
#    from 'long' to 'int' at the top of CxFortranInt.h, or add -D_I4_ to
#    the compiler command line.
#  o using the -Wl,-rpath,, you can embed an "rpath" into the executable. This
#    directory is used by the dynamic linker to look up missing libraries.
#    This makes setting an LD_LIBRARY_PATH for BLAS unnecessary.
#  o If you have a Intel Sandy-Bridge or newer processor, you should use MKL
#    and -lmkl_avx. This can make executions up to twice as fast as when not
#    using avx (note: this does not actually work in 10.x versions of MKL).
#    Similarly, if you have a AMD processor with FMA4 capabilities, you should
#    use the _fma4 versions of ACML with similar benefits. If you have a
#    processor not supporting those features but try to use them anyway,
#    you will get 'illegal instruction' errors when trying to execute the
#    program.
# The following lines are examples of BLAS/LAPACK lines you might take as
# suggestion on how to construct one for your own system:

# export LIBS_LA="-L/opt/acml5.1.0/gfortran64_int64/lib/ -lacml -Wl,-rpath,/opt/acml5.1.0/gfortran64_int64/lib/"
# export LIBS_LA="-L/opt/intel/mkl/10.2.4.032/lib/em64t/ -lmkl_core -lmkl_sequential -lmkl_intel_ilp64 -Wl,-rpath,/opt/intel/mkl/10.2.4.032/lib/"
export LIBS_LA="-L/opt/intel/composerxe/mkl/lib/intel64 -lmkl_gf_ilp64 -lmkl_core -lmkl_avx -lmkl_sequential -Wl,-rpath,/opt/intel/composerxe/mkl/lib/intel64"


# compile bfint with optimizations
c++ -o bfint -O3 -DNDEBUG -DAIC_LOCAL_INCLUDES *.cpp $LIBS_BASE $LIBS_LA
