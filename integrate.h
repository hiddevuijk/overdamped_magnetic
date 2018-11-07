#ifndef GUARD_integrate_h
#define GUARD_integrate_h

#include <vector>

// increments the state of the system over time T,
// in steps of dt

template<class derivobj>
void integrate(
	std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	derivobj& deriv,double& ti, double tf, double dt)
{
	
	while( (ti+dt) < tf) { 
		deriv(r,dr,ti,dt);
	}

	// integrate the remaining time
	if(ti<tf) 
		deriv(r,dr,ti,tf-ti); 
}





#endif
