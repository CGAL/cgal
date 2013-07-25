
/*!
\ingroup PkgNumberTypesConcepts
\cgalConcept

Concept to represent algebraic numbers of degree up to 2 over a `RealEmbeddable` `IntegralDomain` `RT`. 

A model of this concept is associated to this `RT` via `CGAL::Root_of_traits<RT>`, which 
provides `Root_of_2` as a public type. Moreover, `CGAL::Root_of_traits<RT>` provides 
the public type `Root_of_1`, which is the quotient field of `RT`. 
We refer to `Root_of_1` as FT (for field type). 

The model of `RootOf_2 ` is a `RealEmbeddable` `IntegralDomain`, 
which is `ImplicitInteroperable` with `RT`, `FT`. 
In particular, it provides the comparison operators `==, !=, <, >, <=, >=` as well as the `sign` 
and `compare` functions needed to compare elements of types `RootOf_2, RT` and `FT`. 
It also provides all arithmetic operators `+,-,*,/` among elements of type `RootOf_2` as well as mixed forms with `RT` and `FT`. 

However, it is important to note that arithmetic operations among elements of `RootOf_2` 
are only allowed in the special case when they have been constructed from equations having the 
same discriminant, that is, if they are defined in the same algebraic extension of degree 2. 

Besides construction from `int, RT` and `FT` the following functions provide 
special construction for extensions of degree 2: 

- `CGAL::make_root_of_2()` 

- `CGAL::make_sqrt()` 

\cgalRefines `DefaultConstructible` 
\cgalRefines `CopyConstructible` 
\cgalRefines `FromIntConstructible` 
\cgalRefines `ImplicitInteroperable` with RT 
\cgalRefines `ImplicitInteroperable` with FT 

\cgalHasModel `double` (not exact) 
\cgalHasModel `CGAL::Sqrt_extension` 

\sa `CGAL::make_root_of_2<RT>` 
\sa `CGAL::make_sqrt<RT>` 
\sa `CGAL::compute_roots_of_2<RT,OutputIterator>` 
\sa `CGAL::Root_of_traits<RT>` 
\sa `AlgebraicKernelForCircles::PolynomialForCircles_2_2` 
\sa `AlgebraicKernelForCircles` 

*/

class RootOf_2 {
public:

/// \name Operations 
/// Same for operator `-,*,/,!=,<=,>,>=` as well as mixed forms with
/// `RT` and `FT`.
/// @{

/*!
\pre `*this` and `a` are defined in the same extension. 
*/ 
RootOf_2 & operator+=(const RootOf_2& a); 

/*!
\pre `a` and `b` are defined in the same extension. 
*/ 
RootOf_2 operator+(const RootOf_2&a,const RootOf_2& b); 

/*!

*/ 
bool operator==(const RootOf_2&a,const RootOf_2& b); 

/*!

*/ 
bool operator< (const RootOf_2&a,const RootOf_2& b); 

/// @}

}; /* end RootOf_2 */

