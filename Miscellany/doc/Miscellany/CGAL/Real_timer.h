
namespace CGAL {

/*!
\ingroup MiscellanyRef

The class `Real_timer` is a timer class for measuring real time.
A timer `t` of type `Real_timer` is an object with a state. It is
either *running* or it is *stopped*. The state is controlled
with `Real_timer::start()` and `Real_timer::stop()`. The timer counts the
time elapsed since its creation or last reset. It counts only the time
where it is in the running state. The time information is given in seconds.
The timer counts also the number of intervals it was running, i.e.\ it
counts the number of calls of the `Real_timer::start()` member function since the
last reset. If the reset occures while the timer is running it counts as the
first interval.

\cgalHeading{Implementation}

The timer class is based on the C function `gettimeofday()` on
POSIX systems and the C function `_ftime()` on MS Visual C++.
The system calls to these
timers might fail, in which case a warning message will be issued
through the \cgal error handler and the functions return with the
error codes indicated above. The `Real_timer::precision()` method computes the
precision dynamically at runtime at its first invocation.

*/

class Real_timer {
public:

/// \name Creation
/// @{

/*!
state is *stopped*.
*/
Real_timer();

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
bool is_running() const;

/*!
real time in seconds, or 0 if the
underlying system call failed.
*/
double time() const;

/*!
number of start/stop-intervals since
the last reset.
*/
int intervals() const;

/*!
smallest possible time step in seconds,
or -1 if the system call failed.
*/
double precision() const;

/*!
maximal representable time in seconds.
*/
double max const();

/// @}

}; /* end Real_timer */
} /* end namespace CGAL */
