
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the degree 
of a `PolynomialTraits_d::Polynomial_d` with respect to a certain variable. 

The degree of \f$ p\f$ with respect to a certain variable \f$ x_i\f$, 
is the highest power \f$ e\f$ of \f$ x_i\f$ such that the coefficient of \f$ x_i^{e}\f$ in 
\f$ p\f$ is not zero. 

For instance the degree of \f$ p = x_0^2x_1^3+x_1^4\f$ with respect to \f$ x_1\f$ is \f$ 4\f$. 

The degree of the zero polynomial is set to \f$ 0\f$. 
From the mathematical point of view this should 
be \f$ -infinity\f$, but this would imply an inconvenient return type. 

\cgalRefines `AdaptableUnaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::TotalDegree`
\sa `PolynomialTraits_d::DegreeVector`

*/

class PolynomialTraits_d::Degree {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef int result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes the degree of \f$ p\f$ with respect to the outermost variable \f$ x_{d-1}\f$. 
*/ 
result_type operator()(argument_type p); 

/*!
Computes the degree of \f$ p\f$ with respect to variable \f$ x_i\f$. 
\pre \f$ 0 \leq i < d\f$. 

*/ 
result_type operator()(argument_type p, int i); 

/// @}

}; /* end PolynomialTraits_d::Degree */

