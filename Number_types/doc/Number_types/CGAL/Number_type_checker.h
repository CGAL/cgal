
namespace CGAL {

/*!
\ingroup nt_cgal

`Number_type_checker` is a number type whose instances store two numbers 
of types `NT1` and `NT2`. It forwards all arithmetic operations to 
them, and calls the binary predicate `Comparator` to check the equality of 
the instances after each modification, as well as for each comparison. 

This is a debugging tool which is useful when dealing with number types. 


\tparam NT1, NT2 must be a model of the same algebraic structure concept
and they must be `FromDoubleConstructible`. 
\tparam Comparator must be a model of a binary predicate taking `NT1` 
as first argument, and `NT2` as second. The `Comparator` parameter 
has a default value which is a functor calling `operator==` between 
the two arguments. 

\cgalModels `IntegralDomainWithoutDivision` (same as `NT1`) 
\cgalModels `RealEmbeddable` 

\cgalHeading{Operations}

Some operations have a particular behavior documented here. 

*/
template< typename NT1, typename NT2, typename Comparator >
class Number_type_checker {
public:

/// \name Creation 
/// @{

/*!
introduces an uninitialized variable `c`. 
*/ 
Number_type_checker(); 

/*!
introduces the integral value i. 
*/ 
Number_type_checker(int i); 

/*!
introduces the floating point value d. 
*/ 
Number_type_checker(double d); 

/*!
introduces a variable storing the pair `n1, n2`. 
*/ 
Number_type_checker(const NT1 &n1, const NT2 &n2); 

/// @} 

/// \name Operations 
/// @{

/*!
returns a const reference to the object of type `NT1`. 
*/ 
const NT1 & n1() const; 

/*!
returns a const reference to the object of type `NT2`. 
*/ 
const NT2 & n2() const; 

/*!
returns a reference to the object of type `NT1`. 
*/ 
NT1 & n1(); 

/*!
returns a reference to the object of type `NT2`. 
*/ 
NT2 & n2(); 

/*!
calls the `Comparator` binary predicate on the two stored objects 
and returns its result. 
*/ 
bool is_valid() const; 

/// @}

}; /* end Number_type_checker */

/*!
writes `c.n1()` to the ostream `out`. 
\relates Number_type_checker 
*/ 
std::ostream& operator<<(std::ostream& out, const Number_type_checker& c); 

/*!
reads an `NT1` from `in`, then converts it to an `NT2`, 
so a conversion from `NT1` to `NT2` is required here. 
\relates Number_type_checker 
*/ 
std::istream& operator>>(std::istream& in, Number_type_checker& c); 

} /* end namespace CGAL */
