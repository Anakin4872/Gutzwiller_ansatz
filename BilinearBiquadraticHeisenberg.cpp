#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <iomanip>
#include <sstream>
#include <string.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <lapack.h>
#include <matrixtypes.h>


using namespace std;
using namespace ula;

#define Pi 3.14159265358979323846



void ITEp(ComplexMatrix  &Phi, ComplexMatrix &dPhi, int n, int s, double D, double th);//Iterations with periodic boundary conditions
void RK4( double h, ComplexMatrix &Phi,int n,int s, double D, double th);//
double DP(ComplexMatrix &Phi, int n, int s, int m);
double M(ComplexMatrix &Phi, int n, int s);
void Energy(double &E, ComplexMatrix &Phi, int n, int s, double D, double th);

void F_IT( ComplexMatrix &f, ComplexMatrix &dPhi,int n,int s, double D, double th );

int g(int k,int s){
  if(k<0){
    k=s-1;
  }
  return k;
}


int main(){

  int n=3, s=1; //n: value of magnetic projection  m->{-1,0,1}s:number of sites  jmax: max number of iterations
  double Ei=0, Ef=0, h=0.001; //h: step size

for(int th_=-60;th_<=-59;th_+=1){
//#pragma omp parallel for num_threads(3)
	for(int D_=-200;D_<=200;D_+=5){

	double th=th_*Pi/100.;
	double D=D_/100.;

ofstream fout,fEnergyTime,fDensity;

ComplexMatrix f(s,n);//Gutzwiller's coefficients
RealVector norm(s);


//cx_mat f(s,n); f.zeros(); //Gutzwiller's coefficients
//vec norm(s); //norm vector

/*string ssth = to_string(th/Pi);
string ssD = to_string(D);
string sss = to_string(s);
ssth.erase(ssth.find_last_not_of('0') + 1, std::string::npos);
ssD.erase(ssD.find_last_not_of('0') + 1, std::string::npos);
	fEnergyTime.open("EnergyTime"+sss+"sites/EnergyTime th="+ssth+"Pi"+" D="+ssD+".dat", ios::out);
	fout.open("GutzwillerCoeff"+sss+"sites/GutzwillerCoeff th="+ssth+"Pi"+" D="+ssD+".dat", ios::out);*/


/*string ssth = to_string(th/Pi);
string ssD = to_string(D);
string sss = to_string(s);
ssth.erase(ssth.find_last_not_of('0') + 1, std::string::npos);
ssD.erase(ssD.find_last_not_of('0') + 1, std::string::npos);
	fEnergyTime.open("EnergyTime th="+ssth+"Pi"+" D="+ssD+".dat", ios::out);*/
	//fout.open("GutzwillerCoeff th="+ssth+"Pi"+" D="+ssD+".dat", ios::out);

	for(int i=0; i<s; i++){ //Initial conditions: balanced mixture P(m=-1)=P(m=0)=P(m=+1)=1/3
 		f(i,0) = sqrt(0.3);
 		f(i,1) = sqrt(0.3);
 		f(i,2) = sqrt(0.3);

norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) + conj(f(i,2))*f(i,2)));

	if(norm(i) != 1 ){ //Normalization of initial conditions
	f(i,0)=f(i,0)/norm(i);
	f(i,1)=f(i,1)/norm(i);
	f(i,2)=f(i,2)/norm(i);
	}
}

bool B = true;

for(int j=1;B;j++){


	Energy(Ei,f,n,s,D,th);

	RK4(h,f,n,s,D,th); //RK4 iteration

	for(int i=0;i<s; i++){
	norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) +conj(f(i,2))*f(i,2)));

		if(norm(i) != 1 ){
		f(i,0)=f(i,0)/norm(i);
		f(i,1)=f(i,1)/norm(i);
		f(i,2)=f(i,2)/norm(i);
		}

	}

	Energy(Ef,f,n,s,D,th);

	if(abs(Ef-Ei)<=pow(10,-5)){B=false;}

	fEnergyTime<<j<<" "<<Ef<<" "<<endl;

cout<<Ef<<" "<<Ef-Ei<<endl;
}

	fEnergyTime.close();

	//#pragma omp critical (cout)
	//{
	for(int i=0; i<s; i++){
	//fout<<real(f(i,0))<<" "<<real(f(i,1))<<" "<<real(f(i,2))<<endl;
	fDensity<<i<<" "<<real(conj(f(i,0))*f(i,0))<<" "<<real(conj(f(i,1))*f(i,1))<<" "<<real(conj(f(i,2))*f(i,2))<<" "<<real(conj(f(i,0))*f(i,0))+real(conj(f(i,1))*f(i,1))+real(conj(f(i,2))*f(i,2))<<endl;
	}
	//cout<<th/Pi<<" "<<D<<endl;
	//}

	//fout.close();
	fDensity.close();

		}
	}
	return 0;
}

//////////////////////
//////Functions///////
//////////////////////

void Energy(double &E, ComplexMatrix &Phi, int n, int s, double D, double th){
	E=0;
	for(int i=0; i<s; i++){

	E += real( cos(th)*( ( conj(Phi(i,1))*Phi(i,0) + conj(Phi(i,2))*Phi(i,1) )*( conj(Phi((i+1)%s,0))*Phi((i+1)%s,1) + conj(Phi((i+1)%s,1))*Phi((i+1)%s,2) )
	+ ( conj(Phi(i,0))*Phi(i,1) + conj(Phi(i,1))*Phi(i,2) )*( conj(Phi((i+1)%s,1))*Phi((i+1)%s,0) + conj(Phi((i+1)%s,2))*Phi((i+1)%s,1) )
	+ ( - conj(Phi(i,0))*Phi(i,0) + conj(Phi(i,2))*Phi(i,2) )*( -conj(Phi((i+1)%s,0))*Phi((i+1)%s,0) + conj(Phi((i+1)%s,2))*Phi((i+1)%s,2) ) )//first term on the eq.
	+ sin(th)*( conj(Phi(i,2))*Phi(i,0)*conj(Phi((i+1)%s,0))*Phi((i+1)%s,2)
	+ conj(Phi(i,0))*Phi(i,2)*conj(Phi((i+1)%s,2))*Phi((i+1)%s,0)
	+ ( conj(Phi(i,0))*Phi(i,0) + conj(Phi(i,2))*Phi(i,2) )*( conj(Phi((i+1)%s,0))*Phi((i+1)%s,0) + conj(Phi((i+1)%s,2))*Phi((i+1)%s,2) )
	+ ( conj(Phi(i,0))*Phi(i,0) + conj(Phi(i,1))*Phi(i,1) )*( conj(Phi((i+1)%s,1))*Phi((i+1)%s,1) + conj(Phi((i+1)%s,2))*Phi((i+1)%s,2) )
	+ ( conj(Phi(i,1))*Phi(i,1) + conj(Phi(i,2))*Phi(i,2) )*( conj(Phi((i+1)%s,0))*Phi((i+1)%s,0) + conj(Phi((i+1)%s,1))*Phi((i+1)%s,1) )
	- conj(Phi(i,1))*Phi(i,0)*conj(Phi((i+1)%s,1))*Phi((i+1)%s,2)
	- conj(Phi(i,1))*Phi(i,2)*conj(Phi((i+1)%s,1))*Phi((i+1)%s,0)
	- conj(Phi(i,2))*Phi(i,1)*conj(Phi((i+1)%s,0))*Phi((i+1)%s,1)
	- conj(Phi(i,0))*Phi(i,1)*conj(Phi((i+1)%s,2))*Phi((i+1)%s,1) )
	+ D*( conj(Phi(i,0))*Phi(i,0) + conj(Phi(i,2))*Phi(i,2) ) );

	}
}

void ITEp(ComplexMatrix &f, ComplexMatrix &df, int n, int s, double D, double th){ //Iterations with periodic boundary conditions

  for(int i=0; i<s; i++){//Está interpolado el cambio de signo debido al tiempo imaginario.

	df(i,0) = cos(th)*( f(i,0)*( conj(f((i+1)%s,0))*f((i+1)%s,0) + f(g(i-1,s),0)*conj(f(g(i-1,s),0)) - conj(f((i+1)%s,2))*f((i+1)%s,2) - f(g(i-1,s),2)*conj(f(g(i-1,s),2)))
                    + f(i,1)*( conj(f((i+1)%s,1))*f((i+1)%s,0) + conj(f((i+1)%s,2))*f((i+1)%s,1) )
                    + f(i,0)*(f(g(i-1,s),0)*conj(f(g(i-1,s),1)) + f(g(i-1,s),1)*conj(f(g(i-1,s),2))) )
         	+ sin(th)*( f(i,0)*(f((i+1)%s,0)*conj(f((i+1)%s,0)) + f((i+1)%s,2)*conj(f((i+1)%s,2)) + f(g(i-1,s),0)*conj(f(g(i-1,s),0)) + f(g(i-1,s),2)*conj(f(g(i-1,s),2)))
         	          - f(i,1)*f(g(i-1,s),1)*conj(f(g(i-1,s),2))
        	          - f(i,1)*f((i+1)%s,1)*conj(f((i+1)%s,2))
        	          + f(i,2)*f(g(i-1,s),0)*conj(f(g(i-1,s),2))
        	          + f(i,0)*(conj(f(g(i-1,s),1))*f(g(i-1,s),1) + f(g(i-1,s),2)*conj(f(g(i-1,s),2)))
        	          + f(i,0)*(f((i+1)%s,1)*conj(f((i+1)%s,1)) + f((i+1)%s,2)*conj(f((i+1)%s,2)))
        	          + f(i,2)*f((i+1)%s,0)*conj(f((i+1)%s,2)));
        	          + D*f(i,0);

	df(i,1) = cos(th)*( f(i,1)*( f(g(i-1,s),1)*conj(f(g(i-1,s),2)) + f(g(i-1,s),0)*conj(f(g(i-1,s),1)) + conj(f(g(i-1,s),0))*f(g(i-1,s),1) + f(g(i-1,s),2)*conj(f(g(i-1,s),1)))
                    + f(i,0)*( conj(f((i+1)%s,0))*f((i+1)%s,1) + conj(f((i+1)%s,1))*f((i+1)%s,2))
                    + f(i,2)*( conj(f((i+1)%s,1))*f((i+1)%s,0) + conj(f((i+1)%s,2))*f((i+1)%s,1)) )
        	+ sin(th)*( f(i,0)*f((i+1)%s,2)*conj(f((i+1)%s,1))
        	          - f(i,2)*f(g(i-1,s),0)*conj(f(g(i-1,s),1))
        	          - f(i,2)*f((i+1)%s,0)*conj(f((i+1)%s,1))
       	            + f(i,0)*f(g(i-1,s),2)*conj(f(g(i-1,s),1))
        	          + f(i,1)*(f((i+1)%s,0)*conj(f((i+1)%s,0)) + f((i+1)%s,1)*conj(f((i+1)%s,1)) + f(g(i-1,s),1)*conj(f(g(i-1,s),1)) + f(g(i-1,s),2)*conj(f(g(i-1,s),2)))
        	          + f(i,1)*(f((i+1)%s,1)*conj(f((i+1)%s,1)) + f((i+1)%s,2)*conj(f((i+1)%s,2)) + f(g(i-1,s),1)*conj(f(g(i-1,s),1)) + f(g(i-1,s),0)*conj(f(g(i-1,s),0))) );

	df(i,2) = cos(th)*( f(i,1)*( f(i,2)*(conj(f((i+1)%s,2))*f((i+1)%s,2) - conj(f((i+1)%s,0))*f((i+1)%s,0) + f(g(i-1,s),2)*conj(f(g(i-1,s),2)) - f(g(i-1,s),0)*conj(f(g(i-1,s),0)))
                    + conj(f((i+1)%s,0))*f((i+1)%s,1) + conj(f((i+1)%s,1))*f((i+1)%s,2)
                    + f(i,2)*( f(g(i-1,s),1)*conj(f(g(i-1,s),0)) + f(g(i-1,s),2)*conj(f(g(i-1,s),1)))) )
         	+ sin(th)*( f(i,2)*(f((i+1)%s,0)*conj(f((i+1)%s,0)) + f(g(i-1,s),0)*conj(f(g(i-1,s),0)) + f((i+1)%s,2)*conj(f((i+1)%s,2)) + f(g(i-1,s),2)*conj(f(g(i-1,s),2)))
		                - f(i,1)*f((i+1)%s,1)*conj(f((i+1)%s,0))
		                - f(i,1)*f(g(i-1,s),1)*conj(f(g(i-1,s),0))
		                + f(i,0)*f((i+1)%s,2)*conj(f((i+1)%s,0))
		                + f(i,2)*(f((i+1)%s,0)*conj(f((i+1)%s,0)) + f((i+1)%s,1)*conj(f((i+1)%s,1)))
		                + f(i,2)*(f(g(i-1,s),0)*conj(f(g(i-1,s),0)) + f(g(i-1,s),1)*conj(f(g(i-1,s),1)))
		                + f(i,0)*f(g(i-1,s),2)*conj(f(g(i-1,s),0)));
                    + D*f(i,2);
  }

}

void RK4( double h, ComplexMatrix &Phi,int n,int s, double D, double th){
  ComplexMatrix k1(s,n); ComplexMatrix k2(s,n); ComplexMatrix k3(s,n); ComplexMatrix k4(s,n); ComplexMatrix Phi_(s,n);
  Phi_=Phi;

    ITEp(Phi_,k1,n,s,D,th);

    Phi_=Phi+0.5*h*k1;
    ITEp(Phi_,k2,n,s,D,th);

    Phi_=Phi+0.5*h*k2;
    ITEp(Phi_,k3,n,s,D,th);

    Phi_=Phi+h*k3;
    ITEp(Phi_,k4,n,s,D,th);

 Phi +=  h*(k1+2.*k2+2.*k3+k4)/6.; //Average of slopes

 }
