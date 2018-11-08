#ifndef GUARD_deriv_magnetic_h
#define GUARD_deriv_magnetic_h

#include <iostream>

#include <cmath>
#include <vector>
#include "vecmanip.h"
#include "box_muller.h"

struct Deriv {
	public:

	Deriv(int NN, double BB, double omegaa,  
		 Ranq2 ranNRR):
		N(NN),B(BB), omega(omegaa),ranNR(ranNRR)
		  {sqrt2 = std::sqrt(2.);}


	// evolves the state of the system to t+dt
	void  operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double& t, double dt);

	private:
		int N;		// number of particles
		double B;
		double omega;
		Ranq2 ranNR;	// the ranom number generator

		double sqrt2;

		double w(const std::vector<double>& ri) {
					return B*sin(ri[1]*omega);}
		double wp(const std::vector<double>& ri) {
					return B*omega*cos(ri[1]*omega);}

};



void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double& t, double dt)
{
	double sqrt_dt = std::sqrt(dt);
	double etaX, etaY;
	double wi; // magnetic field at position ri
	double wpi;
	double A,B,Ap,Bp;
	for(int i=0;i<N;++i) {

		wi = w(r[i]);
		wpi = wp(r[i]);

		
		B = 1./(1+wi*wi);
		A = wi*B;
		Bp = -2*A*B*wpi;
		Ap = wpi*(1-wi*wi)*B*B;

		etaX = ndist(ranNR)*sqrt_dt*sqrt2;
		etaY = ndist(ranNR)*sqrt_dt*sqrt2;

		dr[i][0] = (Ap + B*etaX)*dt + A*etaY; 
		dr[i][1] = (Bp + B*etaY)*dt - A*etaX;
		dr[i][2] = ndist(ranNR)*sqrt_dt*sqrt2;	

		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

	}
	
	t += dt;
}





#endif


