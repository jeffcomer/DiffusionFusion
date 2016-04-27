ROOT=diffusionFusion
# DON'T USE -ffast-math, it seems to keep things from converging

default: gnu

gnu:
	g++ -O3 -Wall -march=native -mfpmath=sse -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

linux:
	g++ -O3 -Wall -mfpmath=sse -DHAVE_INLINE $(ROOT).C -o $(ROOT)_LINUX -lm -lgsl -lgslcblas -fopenmp

linuxdebug:
	g++ -ggdb -DFUSION_DEBUG  -Wall -DHAVE_INLINE $(ROOT).C -o $(ROOT)_LINUX -lm -lgsl -lgslcblas -fopenmp

stable:
	g++ -O1 -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

debug:
	g++ -ggdb -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

fdebug:
	g++ -DFUSION_DEBUG -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

clean:
	rm $(ROOT)

mac:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

macquick:
	g++ -I/opt/local/include -L/opt/local/lib -O0 -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

mac_debug:
macdebug:
	g++ -I/opt/local/include -L/opt/local/lib -DFUSION_DEBUG -ggdb -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

debug_no_gsl:
	g++ -DNO_GSL -ggdb -Wall $(ROOT).C -o $(ROOT) -lm -fopenmp

dx:
	g++ -O2 -Wall dxToArray2d.C -o dxToArray2d
integ:
	g++ -O2 -Wall cubicIntegrate.C -o cubicIntegrate
insert:
	g++ -O2 -Wall insertDataDx2d1.C -o insertDataDx2d
map:
	g++ -O2 -Wall mapField.C -o mapField

ccgCost:
	g++ -O2 -Wall ccgCost.C -o ccgCost
gnuplot:
	g++ -O2 -Wall dxToGnuplot2d.C -o dxToGnuplot2d
brown:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE brownTown2d.C -o brownTown2d -lm -lgsl -lgslcblas
passage:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE brownTown2dPassage.C -o brownTown2dPassage -lm -lgsl -lgslcblas

fox:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE testFox.C -o testFox -lm -lgsl -lgslcblas

fox1:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE testFox1.C -o testFox1 -lm -lgsl -lgslcblas

smol:
	g++ -O2 -Wall doSmoluchowskiReflect.C -o doSmoluchowskiReflect -lm -fopenmp
smoldebug:
	g++ -ggdb doSmoluchowskiReflect.C -o doSmoluchowskiReflect -lm -fopenmp

crank:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doCrank.C -o doCrank -lm -lgsl -lgslcblas -fopenmp
crankdir:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doCrankDirichlet.C -o doCrankDirichlet -lm -lgsl -lgslcblas -fopenmp
gamma:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF gammaTest.C -o gammaTest -lm -lgsl -lgslcblas
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF gTest.C -o gTest -lm -lgsl -lgslcblas
timefrac:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doTimeFracCrank.C -o doTimeFracCrank -lm -lgsl -lgslcblas -fopenmp
timefrac_debug:
	g++ -ggdb -Wall doTimeFracCrank.C -o doTimeFracCrank -lm -lgsl -lgslcblas -fopenmp
fracdir:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doFracCrankDirichlet.C -o doFracCrankDirichlet -lm -lgsl -lgslcblas -fopenmp
fracdir_debug:
	g++ -ggdb -Wall doFracCrankDirichlet.C -o doFracCrankDirichlet -lm -lgsl -lgslcblas -fopenmp

fracfick:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doFracCrankDirichletFick.C -o doFracCrankDirichletFick -lm -lgsl -lgslcblas -fopenmp

crankfick:
	g++ -I/opt/local/include -L/opt/local/lib -O3 -Wall -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF doCrankDirichletFick.C -o doCrankDirichletFick -lm -lgsl -lgslcblas -fopenmp

