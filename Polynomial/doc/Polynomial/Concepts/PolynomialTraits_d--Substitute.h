
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `Functor` substitutes all variables of a given multivariate 
`PolynomialTraits_d::Polynomial_d` by the values given in the 
iterator range, where begin refers the value for the innermost variable. 

\cgalRefines Assignable
\cgalRefines CopyConstructible
\cgalRefines DefaultConstructible

\cgalHeading{Types}

Note that the `result_type` is the coercion type of the value type of the 
given iterator range and `PolynomialTraits_d::Innermost_coefficient_type`. 
In particular `std::iterator_traits<Input_iterator>::%value_type` must be 
`ExplicitInteroperable` with `PolynomialTraits_d::Innermost_coefficient_type`. 
Hence, it can not be provided as a public type in advance. 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::SubstituteHomogeneous`
\sa  `CGAL::Coercion_traits`

*/

class PolynomialTraits_d::Substitute {
public:

/// \name Operations 
/// @{

/*!

Substitutes each variable of \f$ p\f$ by the values given in the iterator range, 
where begin refers to the innermost variable \f$ x_0\f$. 
\pre (`end-begin` == `PolynomialTraits_d::d`) 

*/ 
template<class Input_iterator> 
result_type operator()(PolynomialTraits_d::Polynomial_d p, 
Input_iterator begin, Input_iterator end); 

/// @}

}; /* end PolynomialTraits_d::Substitute */

