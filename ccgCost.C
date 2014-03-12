///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jeffcomer at gmail>

#include "useful.H"

// The heart of the method: Eq 2 of Comer, Chipot, Gonzalez.
inline double ccgCost(double dt, double dx, double frc, double dif, double gradDif, double beta) {
  // Store the data that can be used to reconstruct gt.
  double gtNumer0 = dx - beta*dif*frc*dt - gradDif*dt;
  double gtVar0 = 2.0*dif*dt;
  // Avoid race conditions.
  //gtNumer = gtNumer0;
  //gtVar = gtVar0;

  //printf("kT %g dt %g\n", 1.0/beta, dt);

  return 0.5*log(2.0*M_PI*gtVar0) + 0.5*(gtNumer0*gtNumer0)/gtVar0;
}

int main(int argc, char* argv[]) {
  if (argc != 7) {
    printf("Usage: %s dt dx frc dif gradDif kT\n", argv[0]);
    exit(0);
  }
  
  double dt = strtod(argv[1],NULL);
  double dx = strtod(argv[2],NULL);
  double frc = strtod(argv[3],NULL);
  double dif = strtod(argv[4],NULL);
  double gradDif = strtod(argv[5],NULL);
  double beta = 1.0/strtod(argv[6],NULL);
 
  double cost = ccgCost(dt, dx, frc, dif, gradDif, beta);

  printf("%.10g\n", cost);
  
  return 0;
}
