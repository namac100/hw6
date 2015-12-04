#include <iostream>
#include <fstream>
using namespace std;

void rungekutta(double* r, double& dt, double& t, double& a, double& b, double& c);
void lorenz(double* r, double* K, double t, double& a, double& b, double& c);

int main(){
  double a = 10.0;
  double b = 28.0;
  double c = 8/3.0;
  double tstart = 0.0;
  double tend = 1.0;
  double dt = 0.1;
  double r[3] = {1.0, 1.0, 1.0}; 
  
  for(double t = tstart; t < tend; t += dt){  
    cout << t << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;  
    rungekutta(r, dt, t, a, b, c);
  }
  
  return 0;
}

void rungekutta(double* r, double& dt, double& t, double& a, double& b, double& c){
  double K1[3], K2[3], K3[3], K4[3], rtemp[3];
  
  lorenz(r, K1, t, a, b, c);
  rtemp[0] = r[0] + 0.5*K1[0];
  rtemp[1] = r[1] + 0.5*K1[1];
  rtemp[2] = r[2] + 0.5*K1[2];
      
  lorenz(rtemp, K2, (t + 0.5*dt), a, b, c);
  rtemp[0] = r[0] + 0.5*K2[0];
  rtemp[1] = r[1] + 0.5*K2[1];
  rtemp[2] = r[2] + 0.5*K2[2];
  
  lorenz(rtemp, K3, (t + 0.5*dt), a, b, c);
  rtemp[0] = r[0] + K3[0];
  rtemp[1] = r[1] + K3[1];
  rtemp[2] = r[2] + K3[2];
    
  lorenz(rtemp, K4, (t + dt), a, b, c);
 
  r[0] += dt*(K1[0]/6.0 + K2[0]/3.0 + K3[0]/3.0 + K4[0]/6.0);
  r[1] += dt*(K1[1]/6.0 + K2[1]/3.0 + K3[1]/3.0 + K4[1]/6.0);
  r[2] += dt*(K1[2]/6.0 + K2[2]/3.0 + K3[2]/3.0 + K4[2]/6.0);
  
}
  

void lorenz(double* r, double* K, double t, double& a, double& b, double& c){  
  K[0] = a*(r[1]-r[0]);
  K[1] = r[0]*(b-r[2])-r[1];
  K[2] = r[0]*r[1] - c*r[2]; 
}