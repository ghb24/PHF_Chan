all:
	scons -Q --jobs=4
opt:
	scons -Q --jobs=4 phf
dbg:
	scons -Q --jobs=4 phfd
clean:
	scons --clean
