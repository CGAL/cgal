
/*!
\ingroup PkgKernelDLinAlgConcepts
\cgalconcept

An instance of data type `Vector` is a vector of variables of 
number type `NT`. Together with the type `Matrix` it realizes 
the basic operations of linear algebra. 

Operations 
-------------- 

The same operations `begin()`, `end()` exist for 
`const_iterator`. 

*/

class Vector {
public:

/// \name Types 
/// @{

/*! 
the ring type of the components. 
*/ 
typedef Hidden_type NT; 

/*! 
the iterator type for accessing components. 
*/ 
typedef Hidden_type iterator; 

/*! 
the const iterator type for accessing 
components. 
*/ 
typedef Hidden_type const_iterator; 

/// @} 

/// \name Creation 
/// @{

/*! 
creates an instance `v` of type 
`Vector`. 
*/ 
Vector(); 

/*! 
creates an instance `v` of type 
`Vector`. `v` is initialized to a vector of dimension \f$ d\f$. 

*/ 
Vector(int d); 

/*! 
creates an instance `v` of 
type `Vector`. `v` is initialized to a vector of dimension 
\f$ d\f$ with entries `x`. 
*/ 
Vector(int d, NT x); 

/*! 
creates an 
instance `v` of type `Vector`; `v` is initialized to the 
vector with entries `set [first,last)`. 
\requires `Forward_iterator` has value type `NT`. 
*/ 
template <class Forward_iterator> 
Vector(Forward_iterator first, Forward_iterator last); 

/// @} 

/// \name Operations 
/// @{

/*! 
returns the dimension of `v`. 
*/ 
int dimension() ; 

/*! 
returns true iff `v` is the zero 
vector. 
*/ 
bool is_zero() ; 

/*! 
returns the \f$ i\f$-th component of `v`. 

\pre \f$ 0\le i \le v.dimension()-1\f$. 
*/ 
NT& operator[](int i) ; 

/*! 
iterator to the first component. 
*/ 
iterator begin() ; 

/*! 
iterator beyond the last component. 
*/ 
iterator end() ; 

/*! 
Addition. 

\pre `v.dimension() == v1.dimension()`. 
*/ 
Vector operator+(const Vector& v1) ; 

/*! 
Subtraction.

\pre `v.dimension() = v1.dimension()`. 
*/ 
Vector operator-(const Vector& v1) ; 

/*! 
Inner Product.

\pre `v.dimension() = v1.dimension()`. 
*/ 
NT operator*(const Vector& v1) ; 

/*! 
Negation. 
*/ 
Vector operator-() ; 

/*! 
Addition plus assignment.

\pre `v.dimension() == v1.dimension()`. 
*/ 
Vector& operator+=(const Vector& v1); 

/*! 
Subtraction plus assignment.

\pre `v.dimension() == v1.dimension()`. 
*/ 
Vector& operator-=(const Vector& v1); 

/*! 
Scalar multiplication plus 
assignment. 
*/ 
Vector& operator*=(const NT& s); 

/*! 
Scalar division plus assignment. 

*/ 
Vector& operator/=(const NT& s); 

/*! 
Component-wise multiplication with number \f$ r\f$. 
*/ 
Vector operator*(const NT& r, const Vector& v); 

/*! 
Component-wise multiplication with number \f$ r\f$. 
*/ 
Vector operator*(const Vector& v, const NT& r); 

/// @}

}; /* end Vector */

