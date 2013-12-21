///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"
#include "DiffusionFusion.H"
#include "CommandLineReader.H"
#include <omp.h>

void printUsage(const char* argv0) {
  DiffusionFusion::printUsage();
  printf("\nUsage: %s config_file [-nt numThreads] [-o outputPrefix]\n", argv0);
  printf("\n");
}

int main(int argc, char* argv[]) {
  printf("*****DiffusionFusion*****\n");
  printf("Written by Jeff Comer <jeffcomer at gmail dot com>\n");
  printf("Implementation of a Bayesian inference scheme for calculating position-dependent diffusivities and free-energy landscapes from MD trajectories.\n");
  printf("The approach is described in Comer, Chipot and Gonzalez-Nilo (2013) J. Chem. Theory Comput., 9(2) pp 876-882.\n\n");
  
  // Initialize the clock.
  double clockInit = omp_get_wtime();

  // Read the parameters and options. 
  CommandLineReader cmd(argc, argv);

  if (cmd.getParamNum() != 1) {
    printUsage(argv[0]);
    fprintf(stderr,"ERROR A configuration file was not specified on the command line.\n");
    exit(-1);
  }

  int numThreads = -1;
  String outPreCmdLine;
  const int optN = cmd.getOptionNum();
  for (int o = 0; o < optN; o++) {
    String opt = cmd.getOption(o);
    String val = cmd.getOptionValue(o);
    if (opt == String("nt")) {
      numThreads = atoi(val.cs());
    } else if(opt == "o") {
      outPreCmdLine = val;
    } else {
      fprintf(stderr,"ERROR Unrecognized option `%s'.\n", opt.cs());
    }
  }

  
  // Get the number of threads if it's specified.
  if (numThreads > 0) {
    omp_set_num_threads(numThreads);
    printf("OMP: Number of OpenMP threads was set to %d.\n", numThreads);
  }
#pragma omp parallel
  {
    printf("Threads: %d\n", omp_get_num_threads());
  }

  // Make the main object.
  DiffusionFusion fusion(cmd.getParam(0), outPreCmdLine);
  printf("\nDone with initialization of DiffusionFusion object.\n");
  printf("Running...\n");
  fusion.run();
  
  // Report the run time.
  double runTime = omp_get_wtime() - clockInit;
  printf("\nRun time: %.10g s\n", runTime);

  return 0;
}
