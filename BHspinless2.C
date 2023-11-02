#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <armadillo>
#include <omp.h>

using namespace std;
using namespace arma;

void ITE(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V );//Iterations with no periodic boundary conditions
void ITEp(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V );//Iterations with periodic boundary conditions
void ITEpV(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V );//Iterations with periodic boundary conditions 									and a density density interaction term
void RK4( double h, cx_mat &Phi,int n,int s, double mu, double t, double V );
void ITEpV2(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V );
int main(){

double V = 0.4;

int MUmin=-20 , MUmax=120+2;  //Min and Max mu values

	ofstream fout2;
	fout2.open("CompleteDataV=0.4.dat");

for(int MU=MUmin ; MU<=MUmax ; MU=MU+1){ //Iteration over all the 'mu' values


int Tmin=0 , Tmax=400+5 ;  //Min and Max t values

//	ofstream fout;
//	string dir  ("Data_mu/"), name ("GA_mu_"), ext(".dat");

//	stringstream ssMU;
//	ssMU << MU << endl;
//	fout.open(dir+name+ssMU.str()+ext, ios::out); //This writtes into a file .dat all data corresponding to a fixed mu value

#pragma omp parallel for num_threads(4)
for(int T=Tmin ; T<=Tmax ; T=T+5){ //Iteration over all the hopping 't' values


int n=3, s=10, jmax=100000; //n:maximun number of particles per site	s:number of sites	jmax: max number of iterations
double h=0.001; //h: step size

	double t=T/1000.;
	double mu=MU/100.;

cx_mat f(s,n);f.zeros(); //Gutzwiller's coefficients

vec norm(s); //norm vector

for(int i=0; i<s; i++){ //Initial conditions
  f(i,0) = 0;
  f(i,1) = sqrt(2)/2;
  f(i,2) = sqrt(2)/2;

	norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) +conj(f(i,2))*f(i,2)));
//cout<<"norma1="<<norm(i)<<endl;

if(norm(i) != 1 ){ //Normalization of initial conditions
	f(i,0)=f(i,0)/norm(i);
	f(i,1)=f(i,1)/norm(i);
	f(i,2)=f(i,2)/norm(i);
	}

//norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) +conj(f(i,2))*f(i,2)));
//cout<<"norma2="<<norm(i)<<endl;

}

//cout<<f<<endl;

for(int j=1;j<=jmax;j++){

RK4(h,f,n,s,mu,t,V); //RK4 iteration

for(int i=0; i<s; i++){
	norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) +conj(f(i,2))*f(i,2)));

//cout<<"norma1="<<norm(i)<<endl;

if(norm(i) != 1 ){
	f(i,0)=f(i,0)/norm(i);
	f(i,1)=f(i,1)/norm(i);
	f(i,2)=f(i,2)/norm(i);
	}
	
//norm(i)=sqrt(real(conj(f(i,0))*f(i,0) + conj(f(i,1))*f(i,1) +conj(f(i,2))*f(i,2)));
//cout<<"norma2="<<norm(i)<<endl;
//cout<<norm<<endl;

}
}

vec Npart(s); //Number of particles per site
vec N2part(s);//Squared number of particles per site

for(int i=0;i<s;i++){
Npart(i)=real(conj(f(i,1))*f(i,1))+2*real(conj(f(i,2))*f(i,2));}
//cout<<Npart<<endl;

double Nns=0; //Expectation value of number of particles
for(int i=0; i<s; i++)
Nns += Npart(i)/s;
//cout<<Nns<<endl; //Total average of number of particles

for(int i=0;i<s;i++){
N2part(i)=real(conj(f(i,1))*f(i,1))+4*real(conj(f(i,2))*f(i,2));}

double SDns=0, N2ns=0; //SDns: Standar deviation of number of particles
for(int i=0; i<s; i++){
N2ns += N2part(i)/s;}
SDns=N2ns-pow(Nns,2);
//cout<<SDns<<endl;

//fout<<MU<<"	"<<T<<"		"<<SDns<<"	"<<Nns<<endl;
fout2<<mu<<"	"<<t<<"		"<<SDns<<"	"<<Nns<<endl;
cout<<mu<<"	"<<t<<"		"<<SDns<<"	"<<Nns<<endl;

}
//fout.close();
	#pragma omp critical (cout)
	{
	cout<<endl;
	fout2<<endl;
	}
}

fout2.close();


	return 0;
}

////////////////////
//Functions//////////
////////////////////

void ITE(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V ){//Iterations with no periodic boundary conditions

	dPhi(0,0) = t*Phi(0,1)*( conj(Phi(1,1))*Phi(1,0) + sqrt(2)*conj(Phi(1,2))*Phi(1,1) );
	dPhi(0,1) = mu*Phi(0,1) + t*Phi(0,0)*( conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) ) + 					  t*Phi(0,2)*sqrt(2)*( conj(Phi(1,1))*Phi(1,0) + sqrt(2)*conj(Phi(1,2))*Phi(1,1) );
	dPhi(0,2) = 2*mu*Phi(0,2) - Phi(0,2) + t*Phi(0,1)*sqrt(2)*( conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) );

  for(int i=1; i<s-1; i++){

	dPhi(i,0) = t*Phi(i,1)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) + conj(Phi(i+1,1))*Phi(i+1,0) + 			    sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,1) = mu*Phi(i,1) + t*Phi(i,0)*( conj(Phi(i-1,0))*Phi(i-1,1) + sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + 			    conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) ) +                                           			    t*Phi(i,2)*sqrt(2)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) +		            	                conj(Phi(i+1,1))*Phi(i+1,0) + sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,2) = 2*mu*Phi(i,2) - Phi(i,2) + t*Phi(i,1)*sqrt(2)*( conj(Phi(i-1,0))*Phi(i-1,1) +                     			    sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) );
}

	dPhi(s-1,0) = t*Phi(s-1,1)*( conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) );
	dPhi(s-1,1) = mu*Phi(s-1,1) + t*Phi(s-1,0)*( conj(Phi(s-2,0))*Phi(s-2,1) +sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) ) + 			      t*Phi(s-1,2)*sqrt(2)*(conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) );
	dPhi(s-1,2) = 2*mu*Phi(s-1,2) - Phi(s-1,2) + t*Phi(s-1,1)*sqrt(2)*( conj(Phi(s-2,0))*Phi(s-2,1) + 			      sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) );

}

void ITEp(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V ){ //Iterations with periodic boundary conditions

	dPhi(0,0) = t*Phi(0,1)*( conj(Phi(s-1,1))*Phi(s-1,0) + sqrt(2)*conj(Phi(s-1,2))*Phi(s-1,1) + conj(Phi(1,1))*Phi(1,0) + 			    sqrt(2)*conj(Phi(1,2))*Phi(1,1) );

	dPhi(0,1) = mu*Phi(0,1) + t*Phi(0,0)*( conj(Phi(s-1,0))*Phi(s-1,1) + sqrt(2)*conj(Phi(s-1,1))*Phi(s-1,2) + 			    conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) ) +                                           			    t*Phi(0,2)*sqrt(2)*( conj(Phi(s-1,1))*Phi(s-1,0) + sqrt(2)*conj(Phi(s-1,2))*Phi(s-1,1) +		            	                conj(Phi(1,1))*Phi(1,0) + sqrt(2)*conj(Phi(1,2))*Phi(1,1) );

	dPhi(0,2) = 2*mu*Phi(0,2) - Phi(0,2) + t*Phi(0,1)*sqrt(2)*( conj(Phi(s-1,0))*Phi(s-1,1) +                     			    sqrt(2)*conj(Phi(s-1,1))*Phi(s-1,2) + conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) );

  for(int i=1; i<s-1; i++){

	dPhi(i,0) = t*Phi(i,1)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) + conj(Phi(i+1,1))*Phi(i+1,0) + 			    sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,1) = mu*Phi(i,1) + t*Phi(i,0)*( conj(Phi(i-1,0))*Phi(i-1,1) + sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + 			    conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) ) +                                           			    t*Phi(i,2)*sqrt(2)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) +		            	                conj(Phi(i+1,1))*Phi(i+1,0) + sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,2) = 2*mu*Phi(i,2) - Phi(i,2) + t*Phi(i,1)*sqrt(2)*( conj(Phi(i-1,0))*Phi(i-1,1) +                     			    sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) );
}

	dPhi(s-1,0) = t*Phi(s-1,1)*( conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) + conj(Phi(0,1))*Phi(0,0) + 			      sqrt(2)*conj(Phi(0,2))*Phi(0,1) );

	dPhi(s-1,1) = mu*Phi(s-1,1) + t*Phi(s-1,0)*( conj(Phi(s-2,0))*Phi(s-2,1) + sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) + 			      conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) ) +                                           			      t*Phi(s-1,2)*sqrt(2)*( conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) +	      	                    conj(Phi(0,1))*Phi(0,0) + sqrt(2)*conj(Phi(0,2))*Phi(0,1) );

	dPhi(s-1,2) = 2*mu*Phi(s-1,2) - Phi(s-1,2) + t*Phi(s-1,1)*sqrt(2)*( conj(Phi(s-2,0))*Phi(s-2,1) +               		      sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) + conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) );

}

void ITEpV(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V ){ //Iterations with periodic boundary conditions 												and a density density interaction ter

  for(int i=1; i<s-1; i++){

	dPhi(i,0) = t*Phi(i,1)*(conj(Phi(i+1,1))*Phi(i+1,0) +  sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,1) = mu*Phi(i,1) + t*Phi(i,0)*( conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) ) + t*Phi(i,2)*sqrt(2)*(conj(Phi(i+1,1))*Phi(i+1,0) + sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) ) - V*Phi(i,1)*(conj(Phi(i+1,1))*Phi(i+1,1) +2.*conj(Phi(i+1,2))*Phi(i+1,2)) ;

	dPhi(i,2) = 2*mu*Phi(i,2) - Phi(i,2) + t*Phi(i,1)*sqrt(2)*(conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) )- 2.*V*Phi(i,2)*(conj(Phi(i+1,1))*Phi(i+1,1) +2.*conj(Phi(i+1,2))*Phi(i+1,2));
}

	dPhi(s-1,0) = t*Phi(s-1,1)*(conj(Phi(0,1))*Phi(0,0)+sqrt(2)*conj(Phi(0,2))*Phi(0,1) );

	dPhi(s-1,1) = mu*Phi(s-1,1) + t*Phi(s-1,0)*(conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) )+ t*Phi(s-1,2)*sqrt(2)*(conj(Phi(0,1))*Phi(0,0) + sqrt(2)*conj(Phi(0,2))*Phi(0,1) ) - V*Phi(s-1,1)*(conj(Phi(0,1))*Phi(0,1)+2.*conj(Phi(0,2))*Phi(0,2)) ;

	dPhi(s-1,2) = 2*mu*Phi(s-1,2) - Phi(s-1,2) + t*Phi(s-1,1)*sqrt(2)*(conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) )- 2.*V*Phi(s-1,2)*(conj(Phi(0,1))*Phi(0,1) +2.*conj(Phi(0,2))*Phi(0,2));

}

void ITEpV2(cx_mat &Phi, cx_mat &dPhi, int n, int s, double mu, double t, double V ){ //Iterations with periodic boundary conditions 												and a density density interaction term

	dPhi(0,0) = (t/2)*Phi(0,1)*( conj(Phi(s-1,1))*Phi(s-1,0) + sqrt(2)*conj(Phi(s-1,2))*Phi(s-1,1) + conj(Phi(1,1))*Phi(1,0) + 			    sqrt(2)*conj(Phi(1,2))*Phi(1,1) );

	dPhi(0,1) = mu*Phi(0,1) + (t/2)*Phi(0,0)*( conj(Phi(s-1,0))*Phi(s-1,1) + sqrt(2)*conj(Phi(s-1,1))*Phi(s-1,2) + 			    conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) ) +                                           			    (t/2)*Phi(0,2)*sqrt(2)*( conj(Phi(s-1,1))*Phi(s-1,0) + sqrt(2)*conj(Phi(s-1,2))*Phi(s-1,1) +		            	                conj(Phi(1,1))*Phi(1,0) + sqrt(2)*conj(Phi(1,2))*Phi(1,1) ) - V*Phi(0,0)*(conj(Phi(s-1,1))*Phi(s-1,1)+2.*conj(Phi(s-1,2))*Phi(s-1,2) + conj(Phi(1,1))*Phi(1,1) +2.*conj(Phi(1,2))*Phi(1,2));

	dPhi(0,2) = 2*mu*Phi(0,2) - Phi(0,2) + (t/2)*Phi(0,1)*sqrt(2)*( conj(Phi(s-1,0))*Phi(s-1,1) +                     			    sqrt(2)*conj(Phi(s-1,1))*Phi(s-1,2) + conj(Phi(1,0))*Phi(1,1) + sqrt(2)*conj(Phi(1,1))*Phi(1,2) )- 2*V*Phi(0,2)*(conj(Phi(s-1,1))*Phi(s-1,1)+2.*conj(Phi(s-1,2))*Phi(s-1,2) + conj(Phi(1,1))*Phi(1,1) +2.*conj(Phi(1,2))*Phi(1,2));

  for(int i=1; i<s-1; i++){

	dPhi(i,0) = (t/2)*Phi(i,1)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) + conj(Phi(i+1,1))*Phi(i+1,0) + 			    sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) );

	dPhi(i,1) = mu*Phi(i,1) + (t/2)*Phi(i,0)*( conj(Phi(i-1,0))*Phi(i-1,1) + sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + 			    conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) ) +                                           			    (t/2)*Phi(i,2)*sqrt(2)*( conj(Phi(i-1,1))*Phi(i-1,0) + sqrt(2)*conj(Phi(i-1,2))*Phi(i-1,1) +		            	                conj(Phi(i+1,1))*Phi(i+1,0) + sqrt(2)*conj(Phi(i+1,2))*Phi(i+1,1) ) - V*Phi(i,1)*(conj(Phi(i-1,1))*Phi(i-1,1)+2.*conj(Phi(i-1,2))*Phi(i-1,2) + conj(Phi(i+1,1))*Phi(i+1,1) +2.*conj(Phi(i+1,2))*Phi(i+1,2)) ;

	dPhi(i,2) = 2*mu*Phi(i,2) - Phi(i,2) + (t/2)*Phi(i,1)*sqrt(2)*( conj(Phi(i-1,0))*Phi(i-1,1) +                     			    sqrt(2)*conj(Phi(i-1,1))*Phi(i-1,2) + conj(Phi(i+1,0))*Phi(i+1,1) + sqrt(2)*conj(Phi(i+1,1))*Phi(i+1,2) )- 2.*V*Phi(i,2)*(conj(Phi(i-1,1))*Phi(i-1,1)+2.*conj(Phi(i-1,2))*Phi(i-1,2) + conj(Phi(i+1,1))*Phi(i+1,1) +2.*conj(Phi(i+1,2))*Phi(i+1,2));
}

	dPhi(s-1,0) = (t/2)*Phi(s-1,1)*( conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) + conj(Phi(0,1))*Phi(0,0) + 			      sqrt(2)*conj(Phi(0,2))*Phi(0,1) );

	dPhi(s-1,1) = mu*Phi(s-1,1) + (t/2)*Phi(s-1,0)*( conj(Phi(s-2,0))*Phi(s-2,1) + sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) + 			      conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) ) +                                           			      (t/2)*Phi(s-1,2)*sqrt(2)*( conj(Phi(s-2,1))*Phi(s-2,0) + sqrt(2)*conj(Phi(s-2,2))*Phi(s-2,1) +	      	                    conj(Phi(0,1))*Phi(0,0) + sqrt(2)*conj(Phi(0,2))*Phi(0,1) ) - V*Phi(s-1,1)*(conj(Phi(s-2,1))*Phi(s-2,1)+2.*conj(Phi(s-2,2))*Phi(s-2,2) + conj(Phi(0,1))*Phi(0,1) +2.*conj(Phi(0,2))*Phi(0,2)) ;

	dPhi(s-1,2) = 2*mu*Phi(s-1,2) - Phi(s-1,2) + (t/2)*Phi(s-1,1)*sqrt(2)*( conj(Phi(s-2,0))*Phi(s-2,1) +               		      sqrt(2)*conj(Phi(s-2,1))*Phi(s-2,2) + conj(Phi(0,0))*Phi(0,1) + sqrt(2)*conj(Phi(0,1))*Phi(0,2) )- 2.*V*Phi(s-1,2)*(conj(Phi(s-2,1))*Phi(s-2,1)+2.*conj(Phi(s-2,2))*Phi(s-2,2) + conj(Phi(0,1))*Phi(0,1) +2.*conj(Phi(0,2))*Phi(0,2));

}

void RK4( double h, cx_mat &Phi,int n,int s, double mu, double t, double V ){
  cx_mat k1(s,n); cx_mat k2(s,n); cx_mat k3(s,n); cx_mat k4(s,n); cx_mat Phi_(s,n);
  Phi_=Phi;

    ITEpV2(Phi_,k1,n,s,mu,t,V);
    
    Phi_=Phi+0.5*h*k1;
    ITEpV2(Phi_,k2,n,s,mu,t,V);

    Phi_=Phi+0.5*h*k2;
    ITEpV2(Phi_,k3,n,s,mu,t,V);

    Phi_=Phi+h*k3;
    ITEpV2(Phi_,k4,n,s,mu,t,V);

 Phi+=h*(k1+2.*k2+2.*k3+k4)/6.; //Average of slopes

 }
