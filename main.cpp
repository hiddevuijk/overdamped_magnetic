
#include "ConfigFile.h"

#include "nr3.h"
#include "ran.h"

#include "vecmanip.h"
#include "init_random.h"
#include "deriv.h"
#include "integrate.h"
#include "density.h"
#include "orientation.h"
#include "flux.h"

#include <iostream>
#include <vector>
#include <vector>
#include <string>


using namespace std;

typedef vector<double> M1;
typedef vector<vector<double> > M2;
typedef vector<vector<vector<double> > > M3;


int main()
{

	int N;		// number of particles
	int navg;	// number of averages
	int seed;	// seed for the random number gen.
	int nbin;	// number of bins
	double a;   // where calculate flux
	double B;	// magnetic field strength
	double w0;	// number of periods of the magnetic field
	double dt;	// time step for the integration
	double tsamp; // time between calculating averages
	double teq;	// equilibration time
	double L;	// box size
	int nprint; // print status every nprint

	
	const string inputName = "input.txt";

	ConfigFile config(inputName);
	N = config.read<int>("N");
	navg = config.read<int>("navg");
	seed = config.read<int>("seed");
	nbin = config.read<int>("nbin");
	a = config.read<double>("a");
	B = config.read<double>("B");
	w0 = config.read<double>("w0");
	dt = config.read<double>("dt");
	tsamp = config.read<double>("tsamp");
	teq = config.read<double>("teq");
	L = config.read<double>("L");
	nprint = config.read<int>("nprint");
	double omega = w0*2*acos(-1.)/L;
	double bs = L/(nbin);


	// density in the xy plane
	M2 rho(nbin,M1(nbin,0.));

	// flux in the xy plane
	// in the x direction
	M2 fx(nbin,M1(nbin,0.));
	// in the y direction
	M2 fy(nbin,M1(nbin,0.));


	// bins
	M1 bins(nbin,0.);
	for(int i=0;i<nbin;++i)
		bins[i] = (i+.5)*bs;


	// the random number generator
	Ranq2 ranNR(seed);

	// vectors describing the state of the system
	M2 r(N,M1(3));
	M2 dr(N,M1(3));


	// the deriv objecs takes care of 
	// a single step of dt
	Deriv deriv(N,B,omega,ranNR);

	// integrates uses the deriv object
	// to increment the system teq in time,
	// in steps of dt

	double t = 0;
	integrate(r,dr,deriv,t,t+teq,dt);


	vector<double> x(navg,0);
	vector<double> y(navg,0);
	vector<double> z(navg,0);
	for(int n=0;n<navg;++n) {
		if(n%nprint==0)	cout << n << '\t' << navg << endl;
		x[n] = r[0][0];
		y[n] = r[0][1];
		z[n] = r[0][2];

		integrate(r,dr,deriv,t,t+tsamp-dt,dt);
		// increment time by dt s.t. 
		// the last distplacement is dt
		deriv(r,dr,t,dt);
		fluxXY(fx,fy,r,dr,dr,bs,L,nbin,a);	


		// calculate density
		density(r,rho,bs,L,nbin);

		
	}

	// normalize the flux and density
	for(int i=0;i<nbin;++i) {
		for(int j=0;j<nbin;++j){
			fx[i][j] /= dt*navg*N*bs*bs;
			fy[i][j] /= dt*navg*N*bs*bs;

			rho[i][j] /= navg*N*bs*bs;
		}
	}



	write_vec("x.dat",x);
	write_vec("y.dat",y);
	write_vec("z.dat",z);
	write_vec("bins.dat",bins);
	write_mat("rho.dat",rho);
	write_mat("fx.dat",fx);
	write_mat("fy.dat",fy);
	
	return 0;

}








