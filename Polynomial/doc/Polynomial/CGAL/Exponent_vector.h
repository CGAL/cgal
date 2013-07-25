
namespace CGAL {

/*!
\ingroup PkgPolynomialClasses

For a given (multivariate) monomial the vector of its exponents is called the 
exponent vector. The class `Exponent_vector` is meant to represent 
such a vector. 

A vector is considered as valid, in case it represents a valid monomial. 
that is, it should not contain negative exponents. 
We decided to use int as the value type, 
since negative exponents may appear in intermediate results. 
The set of exponent vectors with elementwise 
addition forms an <I>Abelian Group</I>. 

Beside the constructors `Exponent_vector` has almost the same interface as an 
`std::vector<int>`. Moreover the comparison is changed such that 
the lexicographic order starts the comparison at the last entry. 
This reflects the fact that the last entry corresponds to the 
outermost variable of a multivariate polynomial. 

\cgalModels `RandomAccessContainer`
\cgalModels `DefaultConstructible`
\cgalModels `Assignable`
\cgalModels `CopyConstructible`
\cgalModels `EqualityComparable`
\cgalModels `LessThanComparable`

\sa `Polynomial_d`
\sa `PolynomialTraits_d`

*/
class Exponent_vector {
public:

/// \name Creation 
/// @{

/*!
%Default constructor.

*/ 
Exponent_vector(); 

/*!
The copy constructor. 

*/ 
Exponent_vector(const Exponent_vector & ev_); 

/*!
Creates a vector containing the given element. 

*/ 
Exponent_vector(int e1); 

/*!
Creates a vector containing the given elements. 

*/ 
Exponent_vector(int e1, int e2); 

/*!
Creates a vector containing the given elements. 

*/ 
Exponent_vector(int e1, int e2, int e3); 

/*!
Creates a vector containing the given elements. 

*/ 
Exponent_vector(int e1, int e2, int e3, int e4); 

/*!

Creates a vector with a copy of the given range. 
\pre `InputIterator` must allow the value type `int`. 

*/ 
template < class InputIterator > 
Exponent_vector(InputIterator begin, InputIterator end); 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Exponent_vector operator+(const Exponent_vector &ev1); 

/*!

*/ 
Exponent_vector operator-(const Exponent_vector &ev1); 

/*!
\pre ev1.size() == ev2.size(). 

*/ 
Exponent_vector 
operator+(const Exponent_vector &ev1, const Exponent_vector &ev2); 

/*!

\pre ev1.size() == ev2.size() 

*/ 
Exponent_vector 
operator-(const Exponent_vector &ev1, const Exponent_vector &ev2); 

/*!

\pre `fo`.size() == ev2.size() 

*/ 
Exponent_vector 
operator+=(const Exponent_vector &ev2); 

/*!

\pre `fo`.size() == ev2.size() 

*/ 
Exponent_vector 
operator-=(const Exponent_vector &ev2); 

/*!

*/ 
bool 
operator==(const Exponent_vector &ev1, const Exponent_vector &ev2); 

/*!

*/ 
bool 
operator!=(const Exponent_vector &ev1, const Exponent_vector &ev2); 

/*!

Lexicographic compare, starting with the <I>last</I> variable. 
*/ 
bool 
operator<(const Exponent_vector &ev1, const Exponent_vector &ev2); 

/// @}

}; /* end Exponent_vector */

/*!
Returns true if all entries of exponent vector `ev` are not negative. 
\relates Exponent_vector 
*/ 
bool is_valid(const Exponent_vector& ev); 


} /* end namespace CGAL */
