
namespace CGAL {

/*!
\ingroup nt_gmp

An object of the class `mpzf`  is a multiple-precision floating-point number which can represent 
numbers of the form \f$ m*2^e\f$, where \f$ m\f$ is an arbitrary precision integer 
based on the <span class="textsc">Gmp</span> library, and \f$ e\f$ 
is of type `int`. This type can be considered exact, even if the 
exponent is not a multiple-precision number. This number type offers 
functionality very similar to `MP_Float` and `Gmpzf` but is faster. 

\cgalModels `EuclideanRing` 
\cgalModels `RealEmbeddable` 

\cgalHeading{Implementation}

\todo need to add something

*/

class mpzf {
public:

/// \name Creation 
/// @{

/*! 
creates a `mpzf` initialized with `0`. 
*/ 
mpzf(); 

/*! 
creates a `mpzf` initialized with `i`. 
*/ 
mpzf(int i); 

/*! 
creates a `mpzf` initialized with `l`. 
*/ 
mpzf(long int l); 

/*! 
creates a `mpzf` initialized with `i`.
\todo to implement
*/ 
mpzf(const Gmpz& i); 

/*! 
creates a `mpzf` initialized with `d`. 
*/ 
mpzf(double d); 

/// @}

}; /* end mpzf */

/*! 
writes a double approximation of `f` to the ostream `out`. 
\relates mpzf 
*/ 
std::ostream& operator<<(std::ostream& out, const mpzf& f); 

/*! 
reads a `double` from `in`, then converts it to a `mpzf`. 
\relates mpzf 
*/ 
std::istream& operator>>(std::istream& in, mpzf& f); 

} /* end namespace CGAL */

