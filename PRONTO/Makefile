# Define the makefile commands

all: inverse forward

# compile the tomography inversion software
inverse: pronto.f
	gfortran pronto.f -O3 -o pronto.exe
	mv -f pronto.exe /hammer/bin/

# compile the forward modeling software for checkerboard tests
# fmod is a "broken" version of pronto that does forward modeling
forward: fmod.f
	gfortran fmod.f -O3 -o fmod.exe
	mv -f fmod.exe /hammer/bin/

clean: 
	rm -f fmod.exe pronto.exe

# NOTE:
# I use gfortran on my Mac. I included the exe file on the FTP site but it may
# not work on your machine. The machine I compiled on is my machine, a Mac
# laptop running snow OSX Yosemite. Gfortran is freely available and has both Mac
# and Linux versions. The old "f77" or "g77" fortran compilers will work as
# well and interchangeably with "gfrotran".

# After compiling these codes make sure to add /hammer/bin/ to your PATH in
# .bashrc or .bash_profile. Or where ever you prefer to compile the binaries.