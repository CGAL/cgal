namespace CGAL {

/*!
\ingroup PkgGenerators
The class `Random` is a random numbers generator. It generates 
uniformly distributed random `bool`, `int` and `double`. 
It can be used as the random number generating function object in the 
\stl algorithm `std::random_shuffle`. 

Instances of `Random` can be seen as input streams. Different 
streams are <I>independent</I> of each other, i.e.\ the sequence of 
numbers from one stream does <I>not</I> depend upon how many numbers 
were extracted from the other streams. At each time, an instance has 
a <I>state</I> that uniquely determines the subsequent numbers being 
produced. 

It can be very useful, e.g. for debugging, to reproduce a sequence of 
random numbers. This can be done by either initialising with a fixed 
seed, or by using the state functions as described below. 

\note A `Random` object is not deterministic when used by several threads at 
the same time, even if a fixed seed has been chosen.

\cgalHeading{Implementation}

We use the boost random library function `boost::rand48` to generate the random 
numbers. 

\sa `CGAL::get_default_random` 

*/

class Random {
public:

/// \name Types 
/// @{

/*!
The State type. 
*/ 
typedef unspecified_type State; 

/// @} 

/// \name Creation 
/// @{

/*!

%Default constructor. The 
seed is chosen based on the system time. 
*/ 
Random( ); 

/*!

Constructor initializing its internal state using `seed`. Equal 
values for `seed` result in equal sequences of random 
numbers. 
*/ 
Random( unsigned int seed); 

/// @} 

/// \name Operations 
/// @{

/*!

returns a random `bool`. 
*/ 
bool get_bool( ); 

/*!

returns a random `int` value from the interval 
\f$[0,2^b)\f$. This is supposed to 
be efficient. 
*/ 
template <int b> int get_bits(); 

/*!

returns a random `int` from the interval 
`[lower,upper)`. 
*/ 
int get_int( int lower, int upper); 

/*!

returns a random `double` from the interval 
`[lower,upper)`. 
*/ 
double get_double( double lower = 0.0, 
double upper = 1.0); 

/// @} 

/// \name Distributions 
/// The following member functions are a 1-to-1 correspondence to some distributions from the boost random library.
/// @{

/*!

returns a random `IntType` from the interval 
`[lower,upper)`. `IntType` can be an integral type 
as `int`, `std::ptrdiff_t`, `std::size_t`,etc. 
\warning In contrast to `get_int` this function may return `upper`.
*/ 
template <typename IntType> IntType uniform_smallint( IntType lower=0, IntType upper=9); 

/*!

returns a random `IntType` from the interval 
`[lower,upper)`. `IntType` can be an integral type 
as `int`, `std::ptrdiff_t`, `std::size_t`,etc. 
\warning In contrast to `get_int` this function may return `upper`.
*/ 
template <typename IntType> IntType uniform_int( IntType lower=0, IntType upper=9); 

/*!

returns a random `RealType` from the interval 
`[lower,upper)`. `RealType` can be `float`, `double`, etc. 
*/ 
template <typename RealType> Realtype uniform_real( RealType lower = 0.0, 
RealType upper = 1.0); 

/*!

returns a random `RealType` from the interval 
`[0,1)`. `RealType` can be `float`, `double`, etc. 
*/ 
template <typename RealType> RealType uniform_01(); 

/*!

returns `random``uniform_int<IntType>( 0, upper-1)`. 
*/ 
template <typename IntType> IntType operator() ( IntType upper); 

/// @} 

/// \name Seed and State Functions 
/// @{

/*!

returns the seed used for initialization. 
*/ 
unsigned int get_seed() const; 

/*!

saves the current internal state in `state`. 
*/ 
void save_state( State& state) const; 

/*!

restores the internal state from `state`. 
*/ 
void restore_state( State const& state); 

/// @} 

/// \name Equality Test 
/// @{

/*!

returns `true`, iff the random object and `random2` have equal 
internal states. 
*/ 
bool operator == ( Random const& random2) const; 

/// @}

}; /* end Random */

/*!
  \ingroup PkgGenerators
  The global function `get_default_random()` returns the default random 
  numbers generator used for the generator functions and classes.
*/
Random &get_default_random();

/*!
  \ingroup PkgGenerators
  \deprecated The variable `default_random` is the default random
  numbers generator used for the generator functions and
  classes. Deprecated. Use `get_default_random()` instead.
*/
extern CGAL::Random default_random;

} /* end namespace CGAL */
