ROOT=diffusionFusion
# DON'T USE -ffast-math, it seems to keep things from converging

default: gnu

gnu:
	g++ -O3 -Wall -march=native -mfpmath=sse $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas -fopenmp

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
