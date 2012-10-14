import os

# make a list of source files: take everything in the current directory
# which looks like C++ or Fortran
FileList = [s for s in os.listdir(".") if os.path.splitext(s)[1] in [".cpp",'.F90','.F']]

# take BLAS/LAPACK library flags from environment variable if present.
# otherwise, substitute a default
if "LAPACK_LINK" in os.environ:
   BlasLapack = os.environ["LAPACK_LINK_I8"]
else:
   if " avx " in open("/proc/cpuinfo").read():
      avx = "-lmkl_avx"
   else:
      avx = ""
   BlasLapack = "-L/opt/intel/composerxe/mkl/lib/intel64 -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core %s -Wl,-rpath,/opt/intel/composerxe/mkl/lib/intel64" % avx

Libs1 = ["boost_program_options"]

def MakeTarget(NameOut, ObjDir, Env ):
   Env.VariantDir(ObjDir, '.', duplicate=0)
   return Env.Program(NameOut, [ObjDir + s for s in FileList])

Incl = "-DAIC_LOCAL_INCLUDES "

FortranFlags = '-fdefault-integer-8 -fdefault-real-8 -fdefault-double-8'
FortranFlagsDbg = FortranFlags + " -g -D_DEBUG "
FortranFlagsOpt = FortranFlags + " -O3 -DNDEBUG "

LinkFlags = BlasLapack

dbg = MakeTarget("phfd", "build/debug/", Environment( CCFLAGS = "-g -O0 -Wall -Wno-unknown-pragmas -D_DEBUG -DFCI_NO_MAIN " + Incl, F77FLAGS = FortranFlagsDbg, F90FLAGS = FortranFlagsDbg, LINKFLAGS = LinkFlags, LIBS = Libs1))
opt = MakeTarget("phf", "build/release/", Environment( CCFLAGS = "-O3 -Wall -Wno-unknown-pragmas -DNDEBUG -DFCI_NO_MAIN " + Incl, F77FLAGS = FortranFlagsOpt, F90FLAGS = FortranFlagsOpt, LINKFLAGS = LinkFlags, LIBS = Libs1))

Default(dbg)
Default(opt)

