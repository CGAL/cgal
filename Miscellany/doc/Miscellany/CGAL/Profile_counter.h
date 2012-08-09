
namespace CGAL {

/*!
\ingroup PkgProfilingTools

The class `Profile_counter` provides a way to count the number of times a given 
line of code is executed during the execution of a program, and print this 
number at the end of the execution of the program. Such counters can be added 
at critical places in the code, and at the end of the execution of a program, 
the count is printed on `std::cerr`, together with an identification 
string passed to the constructor. The macro `CGAL_PROFILER` can be 
used to conveniently place these counters anywhere. They are disabled by 
default and activated by the global macro `CGAL_PROFILE`. 

Operations 
-------------- 

`\#define CGAL_PROFILER(MSG)` 
If `CGAL_PROFILE` is not defined, then `CGAL_PROFILER` is defined 
to an empty statement. Otherwise, it is defined to `static 
CGAL::Profile_counter tmp(MSG); ++tmp;`. 

\cgalexample{Profile_counter.cpp} 

will print at exit: 

\code{.cpp} 

[CGAL::Profile_counter] 10 iterations of the for-loop 

\endcode 

*/

class Profile_counter {
public:

/// \name Creation 
/// @{

/*! 
The internal counter is 
initialized to 0, and the string `s` is stored for further printing 
by the destructor. 
*/ 
Profile_counter(std::string s); 

/*! 
The value of the counter is printed 
to `std::cerr` together with the string. 
*/ 
~Profile_counter(); 

/// @} 

/// \name Operations 
/// @{

/*! 
Increments the internal counter. 
*/ 
void operator++(); 

/// @}

}; /* end Profile_counter */
} /* end namespace CGAL */
