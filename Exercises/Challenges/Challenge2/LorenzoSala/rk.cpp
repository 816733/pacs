#include "rk.hpp"
#include <cmath>
#include <algorithm> // for max
namespace ODE
{
  // f(t,y)
  
// step for RK 45
double rk45_step
  (
   std::function<double (double const &, double const &)> const & f,
   double const & y0, 
   double const & t0, 
   double const & h, 
   double & error
   )
  {
    // Butcher array (row 1 is trivial)
    constexpr double a21 = 1./4.;
    constexpr double c2  = a21;
    constexpr double a31 = 3./32;
    constexpr double a32 = 9./32;
    constexpr double c3  = a31 + a32;
    constexpr double a41 = 1932./2197.;
    constexpr double a42 =-7200./2197.;
    constexpr double a43 = 7296./2197.;
    constexpr double c4  = a41 + a42 + a43;
    constexpr double a51 = 439./216.;
    constexpr double a52 =- 8.;
    constexpr double a53 =  3680./513.;
    constexpr double a54 =-845./4104.;
    constexpr double c5  = a51 + a52 + a53 + a54;
    constexpr double a61 =-8./27.;
    constexpr double a62 = 2.;
    constexpr double a63 =-3544./2565.;
    constexpr double a64 = 1859./4104.;
    constexpr double a65 =-11./40.;
    constexpr double c6  = a61 + a62 + a63 + a64 + a65;
    
    constexpr double b41 = 25./216.;
    //    constexpr double b42 = 0.;
    constexpr double b43 = 1408./2565.;
    constexpr double b44 = 2197./4104.;
    constexpr double b45 =-1./5.;
    constexpr double b51 = 16./135.;
    //    constexpr double b52 = 0.;
    constexpr double b53 = 6656./12825.;
    constexpr double b54 = 28561./56430.;
    constexpr double b55 =-9./50.;
    constexpr double b56 = 2./55.;
    
    double F1 = h * f(t0, y0);
    double F2 = h * f(t0 + c2 * h, y0 + a21 * F1);
    double F3 = h * f(t0 + c3 * h, y0 + a31 * F1 + a32 * F2);
    double F4 = h * f(t0 + c4 * h, y0 + a41 * F1 + a42 * F2 + a43 * F3);
    double F5 = h * f(t0 + c5 * h, y0 + a51 * F1 + a52 * F2 + a53 * F3 + a54 * F4 );
    double F6 = h * f(t0 + c6 * h, y0 + a61 * F1 + a62 * F2 + a63 * F3 + a64 * F4 + a65 * F5);

    double y4 =   y0 + b41 * F1 + b43 * F3 + b44 * F4 + b45 * F5;
    double y5 =   y0 + b51 * F1 + b53 * F3 + b54 * F4 + b55 * F5 + b56 * F6;
    error = std::abs(y5 - y4);
    return y5;
  }

// step for RK 23
double rk23_step
  (
   std::function<double (double const &, double const &)> const & f,
   double const & y0, 
   double const & t0, 
   double const & h, 
   double & error
   )
  {
    // Butcher array (row 1 is trivial)
    constexpr double a21 = 1./2.;
    constexpr double c2  = a21;
    constexpr double a31 = 0;
    constexpr double a32 = 3./4.;
    constexpr double c3  = a32;
    constexpr double a41 = 2./9.;
    constexpr double a42 = 1./3.;
    constexpr double a43 = 4./9.;
    constexpr double c4  = a41 + a42 + a43;
    
    constexpr double b41 = 2./9.;
    constexpr double b42 = 1./3.;
    constexpr double b43 = 4./9.;
    // constexpr double b44 = 0;

    constexpr double b51 = 7./24.;
    constexpr double b52 = 1./4.;
    constexpr double b53 = 1./3.;
    constexpr double b54 = 1./8.;
    

    double F1 = h * f(t0, y0);
    double F2 = h * f(t0 + c2 * h, y0 + a21 * F1);
    double F3 = h * f(t0 + c3 * h, y0 + a31 * F1 + a32 * F2);
    double F4 = h * f(t0 + c4 * h, y0 + a41 * F1 + a42 * F2 + a43 * F3);
    
    double y2 =   y0 + b41 * F1 + b42 * F2 + b43 * F3;
    double y3 =   y0 + b51 * F1 + b52 * F2 + b53 * F3 + b54 * F4;
    error = std::abs(y3 - y2);
    return y3;
  }

  std::vector<std::pair<double,double>> 
    rk(std::function<double (double const &, double const &)> const & dy,
	 Parametri const & data,
	 int & status,
	 Flag const & flag)
  {
    double t0 = data.t0;
    double T = data.T;
    double y0 = data.y0;
    double h_initial = data.h_init;
    double h_max = (T-t0)/4.0;
    double final_error = data.errorDesired;
    status=0;
    const std::size_t maxReduction = data.maxsteps;
    // parameters for decreasing/increasing time step
    double const c1 = 0.1;
    // I need to have a sufficient decrease of the local error
    // to allow time step coarsening
    double c2;
    if (flag.RKmethod == RK45)
    	c2=1./64.;
    else 
	c2=1./16.;

    double length=T-t0;
    //! Make sure that h allows to reach T
    std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
    double h=length/initialNSteps;
    // To avoid underflow we need in any case to limit the time step to a positive number
    // Here I allow h to become 128 time smaller than that giving the maximal number of steps
    double h_min = length/(128*data.maxsteps);
    // Some counters
    std::size_t count(0);
    std::size_t stepsCounter(0);
    // Initial data
    double time(t0);
    double y(y0);
    double errorPerTimeStep=final_error/initialNSteps;
    if (initialNSteps>=data.maxsteps) throw std::runtime_error("RK45: initial time step h too small!");
    std::vector<std::pair<double,double>> solution;
    solution.emplace_back(std::make_pair(t0,y0));
    double localError;
    double newy;
    double s;
    while (time<T && stepsCounter < data.maxsteps)
      {
	//Do a step
	//adjust h if needed for the last step
	if (time + h > T) h = T-time;
	// Choose the right method
	if (flag.RKmethod == RK45)
		newy = rk45_step(dy,y,time,h,localError);
	else 
		newy = rk23_step(dy,y,time,h,localError);

	while (localError > c1*errorPerTimeStep && count<maxReduction)
	  {
	    if( h> h_min)
	      {
		if(flag.timestep == CHTS){ // half time step
			h /=2;
			errorPerTimeStep /=2;
		}	
		else if(flag.timestep == OTS){ // optimal step
			s = std::sqrt(std::sqrt((0.84*(errorPerTimeStep*h/(2*localError)))));
			h *= s;
			errorPerTimeStep *= s;
		}
		
		++count;
		// Choose the right method
		if (flag.RKmethod == RK45)
			newy = rk45_step(dy,y,time,h,localError);
		else 
			newy = rk23_step(dy,y,time,h,localError);
	      }
	    else status=3;
	    break;
	  }
	if (count>=maxReduction)status=1;
	count=0;
	//! advance
	y = newy;
	time +=h;
	++stepsCounter;
	solution.emplace_back(std::make_pair(time,y));
	//! check if we reached end
	if(localError<c2*errorPerTimeStep && h<h_max)
	  {
	    // Double step
	    h *= 2;
	    errorPerTimeStep *= 2;
	  }
      }
    //handle exceptions
    if(stepsCounter>= data.maxsteps && time < T)
      {
	status=2;
	// Choose the right method
	if (flag.RKmethod == RK45)
		throw std::runtime_error("RK45: Max number of time steps exceeded");	
	else 
		throw std::runtime_error("RK23: Max number of time steps exceeded");
      }
    return solution;
  }

}// end namespace
