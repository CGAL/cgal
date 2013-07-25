
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `AlgebraicStructureTraits` reflects the algebraic structure 
of an associated type `Type`. 

Depending on the concepts that `Type` fulfills, 
it contains various functors and descriptive tags. 
Moreover it gives access to the several possible 
algebraic operations within that structure. 

\sa `IntegralDomainWithoutDivision`
\sa `IntegralDomain`
\sa `UniqueFactorizationDomain`
\sa `EuclideanRing`
\sa `Field`
\sa `FieldWithSqrt`
\sa `FieldWithKthRoot`
\sa `FieldWithRootOf`
\sa `CGAL::Integral_domain_without_division_tag`
\sa `CGAL::Integral_domain_tag`
\sa `CGAL::Unique_factorization_domain_tag`
\sa `CGAL::Euclidean_ring_tag`
\sa `CGAL::Field_tag`
\sa `CGAL::Field_with_sqrt_tag`
\sa `CGAL::Field_with_kth_root_tag`
\sa `CGAL::Field_with_root_of_tag`

\cgalHasModel `CGAL::Algebraic_structure_traits<T>`

*/

class AlgebraicStructureTraits {
public:

/// \name Types 
/// A model of `AlgebraicStructureTraits` is supposed to provide:
/// @{

/*!
The associated type. 
*/ 
typedef unspecified_type Type; 

/*!
Tag indicating the algebraic structure of the associated type. 

<TABLE CELLSPACING=5 > 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
Tag is: 
<TD ALIGN=LEFT NOWRAP> 
`Type` is model of: 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Null_tag` 
<TD ALIGN=LEFT NOWRAP> 
no algebraic concept 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Integral_domain_without_division_tag` 
<TD ALIGN=LEFT NOWRAP> 
`IntegralDomainWithoutDivision` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Integral_domain_tag` 
<TD ALIGN=LEFT NOWRAP> 
`IntegralDomain` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Unique_factorization_domain_tag` 
<TD ALIGN=LEFT NOWRAP> 
`UniqueFactorizationDomain` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Euclidean_ring_tag` 
<TD ALIGN=LEFT NOWRAP> 
`EuclideanRing` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Field_tag` 
<TD ALIGN=LEFT NOWRAP> 
`Field` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Field_with_sqrt_tag` 
<TD ALIGN=LEFT NOWRAP> 
`FieldWithSqrt` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Field_with_kth_root_tag` 
<TD ALIGN=LEFT NOWRAP> 
`FieldWithKthRoot` 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
`CGAL::Field_with_root_of_tag` 
<TD ALIGN=LEFT NOWRAP> 
`FieldWithRootOf` 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
</TABLE> 

*/ 
typedef unspecified_type Algebraic_category; 

/*!
Tag indicating whether `Type` is exact. 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 

An algebraic structure is considered exact, if all operations 
required by its concept are computed such that a comparison 
of two algebraic expressions is always correct. 
The exactness covers only those operations that are required by 
the algebraic structure concept. 

e.g. an exact `Field` may have a `Sqrt` functor that 
is not exact. 

*/ 
typedef unspecified_type Is_exact; 

/*!
Tag indicating whether `Type` is numerical sensitive. 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 

An algebraic structure is considered as numerically sensitive, 
if the performance of the type is sensitive to the condition 
number of an algorithm. 

*/ 
typedef unspecified_type Is_numerical_sensitive; 

/*!
This type specifies the return type of the predicates provided 
by this traits. The type must be convertible to `bool` and 
typically the type indeed maps to `bool`. However, there are also 
cases such as interval arithmetic, in which it is `Uncertain<bool>` 
or some similar type. 

*/ 
typedef unspecified_type Boolean; 

/// @} 

/// \name Functors 
/// In case a functor is not provided, it is set to `CGAL::Null_functor`.
/// @{

/*!

A model of `AlgebraicStructureTraits_::IsZero`. 

Required by the concept `IntegralDomainWithoutDivision`. 
In case `Type` is also model of `RealEmbeddable` this is a 
model of `RealEmbeddableTraits_::IsZero`. 

*/ 
typedef unspecified_type Is_zero; 

/*!

A model of `AlgebraicStructureTraits_::IsOne`. 

Required by the concept `IntegralDomainWithoutDivision`. 

*/ 
typedef unspecified_type Is_one; 

/*!

A model of `AlgebraicStructureTraits_::Square`. 

Required by the concept `IntegralDomainWithoutDivision`. 

*/ 
typedef unspecified_type Square; 

/*!

A model of `AlgebraicStructureTraits_::Simplify`. 

Required by the concept `IntegralDomainWithoutDivision`. 

*/ 
typedef unspecified_type Simplify; 

/*!

A model of `AlgebraicStructureTraits_::UnitPart`. 

Required by the concept `IntegralDomainWithoutDivision`. 

*/ 
typedef unspecified_type Unit_part; 

/*!

A model of `AlgebraicStructureTraits_::IntegralDivision`. 

Required by the concept `IntegralDomain`. 

*/ 
typedef unspecified_type Integral_division; 

/*!

A model of `AlgebraicStructureTraits_::Divides`. 

Required by the concept `IntegralDomain`. 

*/ 
typedef unspecified_type Divides; 

/*!

A model of `AlgebraicStructureTraits_::IsSquare`. 

Required by the concept `IntegralDomainWithoutDivision`. 

*/ 
typedef unspecified_type Is_square; 

/*!

A model of `AlgebraicStructureTraits_::Gcd`. 

Required by the concept `UniqueFactorizationDomain`. 

*/ 
typedef unspecified_type Gcd; 

/*!

A model of `AlgebraicStructureTraits_::Mod`. 

Required by the concept `EuclideanRing`. 

*/ 
typedef unspecified_type Mod; 

/*!

A model of `AlgebraicStructureTraits_::Div`. 

Required by the concept `EuclideanRing`. 

*/ 
typedef unspecified_type Div; 

/*!

A model of `AlgebraicStructureTraits_::DivMod`. 

Required by the concept `EuclideanRing`. 

*/ 
typedef unspecified_type Div_mod; 

/*!

A model of `AlgebraicStructureTraits_::Inverse`. 

Required by the concept `Field`. 

*/ 
typedef unspecified_type Inverse; 

/*!

A model of `AlgebraicStructureTraits_::Sqrt`. 

Required by the concept `FieldWithSqrt`. 

*/ 
typedef unspecified_type Sqrt; 

/*!

A model of `AlgebraicStructureTraits_::KthRoot`. 

Required by the concept `FieldWithKthRoot`. 

*/ 
typedef unspecified_type Kth_root; 

/*!

A model of `AlgebraicStructureTraits_::RootOf`. 

Required by the concept `FieldWithRootOf`. 

*/ 
typedef unspecified_type Root_of; 

/// @}

}; /* end AlgebraicStructureTraits */

