
namespace CGAL {

/*!
\ingroup PkgProfilingTools

The class `Real_timer` is a timer class for measuring real time. 
A timer `t` of type `Real_timer` is an object with a state. It is 
either <I>running/</I> or it is <I>stopped</I>. The state is controlled 
with `t.start()` and `t.stop()`. The timer counts the 
time elapsed since its creation or last reset. It counts only the time 
where it is in the running state. The time information is given in seconds. 
The timer counts also the number of intervals it was running, i.e. it 
counts the number of calls of the `start()` member function since the 
last reset. If the reset occures while the timer is running it counts as the 
first interval. 

Implementation 
-------------- 

The timer class is based in the C function `gettimeofday()` on 
POSIX systems and the C function `_ftime()` on MS Visual C++. 
The system calls to these 
timers might fail, in which case a warning message will be issued 
through the \cgal error handler and the functions return with the 
error codes indicated above. The `precision` method computes the 
precision dynamically at runtime at its first invocation. 

*/

class Real_timer {
public:

/// \name Creation 
/// @{

/*! 
state is <I>stopped.</I> 
*/ 
Real_timer(); 

/// @} 

/// \name Operations 
/// @{

/*! 
\pre state is <I>stopped</I>. 
*/ 
void start(); 

/*! 
\pre state is <I>running</I>. 
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
real time in seconds, or 0 if the 
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

}; /* end Real_timer */
} /* end namespace CGAL */
