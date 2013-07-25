
namespace CGAL {

/*!
\ingroup PkgPolynomialClasses

This class template provides a convenient way to obtain the type representing a multivariate polynomial with `d` variables, where `T` is the innermost coefficient type. In case `T` happens to be a `CGAL::Polynomial` the generator will add `d` variables to `T`. 

`T` must be a model of `IntegralDomainWithoutDivision`. 

`d` must be of type int. 

\sa `CGAL::Polynomial<Coeff>` 

*/
template< typename T, typename d >
class Polynomial_type_generator {
public:

/// \name Types 
/// @{

/*!
The generated type. 
*/ 
typedef unspecified_type Type; 

/// @}

}; /* end Polynomial_type_generator */
} /* end namespace CGAL */
