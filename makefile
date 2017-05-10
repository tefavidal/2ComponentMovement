all:SimGoldbeter
SimGoldbeter:	
	gfortran -o 2CompCellMov main-Source.f anfang.f rs-Source2D.f ZeroPillar.f \
	out.f ODE-Merson.f ic-Source.f StartingTime.f flow.f \
	-L/usr/bmp/pgplot-5.2/ -lpgplot \
	-L/usr/bmp/slatec-4.1/lib -lslatec \
	-L/usr/bmp/lapack-3.4.0 \
	-lX11 \
	-lpng
