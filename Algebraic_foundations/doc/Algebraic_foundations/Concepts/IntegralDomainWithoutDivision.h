
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalconcept

This is the most basic concept for algebraic structures considered within CGAL. 

A model `IntegralDomainWithoutDivision` represents an integral domain, 
i.e.\ commutative ring with 0, 1, +, * and unity free of zero divisors. 

<B>Note:</B> A model is not required to offer the always well defined integral division. 

It refines `Assignable`, `CopyConstructible`, `DefaultConstructible` 
and `FromIntConstructible`. 

It refines `EqualityComparable`, where equality is defined w.r.t. 
the ring element being represented. 

The operators unary and binary plus +, unary and binary minus -, 
multiplication * and their compound forms +=, -=, *= are required and 
implement the respective ring operations. 

Moreover, `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >` is a model of 
`AlgebraicStructureTraits` providing: 

- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Algebraic_type` derived from `CGAL::Integral_domain_without_division_tag` 
- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Is_zero`  which is a model of `AlgebraicStructureTraits::IsZero`
- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Is_one`  which is a model of `AlgebraicStructureTraits::IsOne`
- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Square`  which is a model of `AlgebraicStructureTraits::Square`
- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Simplify`  which is a model of `AlgebraicStructureTraits::Simplify`
- `CGAL::Algebraic_structure_traits< IntegralDomainWithoutDivision >::Unit_part`  which is a model of `AlgebraicStructureTraits::UnitPart`

\refines `Assignable` 
\refines `CopyConstructible` 
\refines `DefaultConstructible` 
\refines `EqualityComparable` 
\refines `FromIntConstructible` 

\sa ::IntegralDomainWithoutDivision 
\sa ::IntegralDomain 
\sa ::UniqueFactorizationDomain 
\sa ::EuclideanRing 
\sa ::Field 
\sa ::FieldWithSqrt 
\sa ::FieldWithKthRoot 
\sa ::FieldWithRootOf 
\sa ::AlgebraicStructureTraits 

*/

class IntegralDomainWithoutDivision {
public:

/// \name Operations 
/// @{

/*! 
unary plus 
*/ 
IntegralDomainWithoutDivision 
operator+(const IntegralDomainWithoutDivision &a); 

/*! 
unary minus 
*/ 
IntegralDomainWithoutDivision 
operator-(const IntegralDomainWithoutDivision &a); 

/*! 

*/ 
IntegralDomainWithoutDivision 
operator+(const IntegralDomainWithoutDivision &a, 
const IntegralDomainWithoutDivision &b); 

/*! 

*/ 
IntegralDomainWithoutDivision 
operator-(const IntegralDomainWithoutDivision &a, 
const IntegralDomainWithoutDivision &b); 

/*! 

*/ 
IntegralDomainWithoutDivision 
operator*(const IntegralDomainWithoutDivision &a, 
const IntegralDomainWithoutDivision &b); 


/*! 

*/ 
IntegralDomainWithoutDivision 
operator+=(const IntegralDomainWithoutDivision &b); 

/*! 

*/ 
IntegralDomainWithoutDivision 
operator-=(const IntegralDomainWithoutDivision &b); 

/*! 

*/ 
IntegralDomainWithoutDivision 
operator*=(const IntegralDomainWithoutDivision &b); 

/*! 
The `result_type` is convertible to `bool`. 
*/ 
result_type 
operator==(const IntegralDomainWithoutDivision &a, const IntegralDomainWithoutDivision &b); 

/*! 
The `result_type` is convertible to `bool`. 
*/ 
result_type 
operator!=(const IntegralDomainWithoutDivision &a, const IntegralDomainWithoutDivision &b); 

/// @}

}; /* end IntegralDomainWithoutDivision */

