
namespace CGAL {

/*!
\ingroup PkgProfilingTools

The class `Timer` is a timer class for measuring user process time. 
A timer `t` of type `Timer` is an object with a state. It is 
either *running* or it is *stopped*. The state is controlled 
with `Timer::start()` and `Timer::stop()`. The timer counts the 
time elapsed since its creation or last reset. It counts only the time 
where it is in the running state. The time information is given in seconds. 
The timer counts also the number of intervals it was running, i.e. it 
counts the number of calls of the `Timer::start()` member function since the 
last reset. If the reset occures while the timer is running it counts as the 
first interval. 

### Implementation ###

The timer class is based on the C function `std::clock()` on 
PC systems and the C function `getrusage()` on standard 
POSIX systems. The counter for the `std::clock()` based 
solution might wrap around (overflow) after only about 36 
minutes. This won't happen on POSIX systems. The system calls to these 
timers might fail, in which case a warning message will be issued 
through the \cgal error handler and the functions return with the 
error codes indicated above. The `Timer::precision()` method computes the 
precision dynamically at runtime at its first invocation. 

*/

class Timer {
public:

/// \name Creation 
/// @{

/*! 
state is *stopped*. 
*/ 
Timer(); 

/// @} 

/// \name Operations 
/// @{

/*! 
\pre state is *stopped*. 
*/ 
void start(); 

/*! 
\pre state is *running*. 
*/ 
void stop (); 

/*! 
reset timer to zero. The state is unaffected. 
*/ 
void reset(); 

/*! 
`true` if the current state is running. 
*/ 
bool is_running(); 

/*! 
user process time in seconds, or 0 if the 
underlying system call failed. 
*/ 
double time(); 

/*! 
number of start/stop-intervals since 
the last reset. 
*/ 
int intervals(); 

/*! 
smallest possible time step in seconds, 
or -1 if the system call failed. 
*/ 
double precision(); 

/*! 
maximal representable time in seconds. 
*/ 
double max(); 

/// @}

}; /* end Timer */
} /* end namespace CGAL */
