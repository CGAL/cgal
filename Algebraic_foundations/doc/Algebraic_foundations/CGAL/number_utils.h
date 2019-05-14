namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `abs()` returns the absolute value of a number. 

The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::Abs` 

*/
template <class NT> NT abs(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `compare()` compares the first argument with respect to 
the second, i.e.\ it returns `CGAL::LARGER` if \f$ x\f$ is larger then \f$ y\f$. 

In case the argument types `NT1` and `NT2` differ, 
`compare` is performed with the semantic of the type determined via 
`Coercion_traits`. 
The function is defined if this type 
is a model of the `RealEmbeddable` concept. 

The `result_type` is convertible to `CGAL::Comparison_result`. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::Compare` 
*/
template <class NT1, class NT2> 
result_type compare(const NT &x, const NT &y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `div()` computes the integral quotient of division 
with remainder. 

In case the argument types `NT1` and `NT2` differ, 
the `result_type` is determined via `Coercion_traits`. 

Thus, the `result_type` is well defined if `NT1` and `NT2` 
are a model of `ExplicitInteroperable`. 

The actual `div` is performed with the semantic of that type. 

The function is defined if `result_type` 
is a model of the `EuclideanRing` concept. 

\sa `EuclideanRing` 
\sa `AlgebraicStructureTraits_::Div` 
\sa `CGAL::mod()`
\sa `CGAL::div_mod()`

*/
template< class NT1, class NT2> 
result_type 
div(const NT1& x, const NT2& y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

computes the quotient \f$ q\f$ and remainder \f$ r\f$, such that \f$ x = q*y + r\f$ 
and \f$ r\f$ minimal with respect to the Euclidean Norm of the 
`result_type`.

The function `div_mod()` computes the integral quotient and remainder of 
division with remainder. 

In case the argument types `NT1` and `NT2` differ, 
the `result_type` is determined via `Coercion_traits`. 

Thus, the `result_type` is well defined if `NT1` and `NT2` 
are a model of `ExplicitInteroperable`. 

The actual `div_mod` is performed with the semantic of that type. 

The function is defined if `result_type` 
is a model of the `EuclideanRing` concept. 

\sa `EuclideanRing` 
\sa `AlgebraicStructureTraits_::DivMod`
\sa `CGAL::mod()`
\sa `CGAL::div()`

*/
template <class NT1, class NT2> 
void 
div_mod(const NT1& x, const NT2& y, result_type& q, result_type& r);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `gcd()` computes the greatest common divisor of two values. 

In case the argument types `NT1` and `NT2` differ, 
the `result_type` is determined via `Coercion_traits`. 

Thus, the `result_type` is well defined if `NT1` and `NT2` 
are a model of `ExplicitInteroperable`. 

The actual `gcd` is performed with the semantic of that type. 

The function is defined if `result_type` 
is a model of the `UniqueFactorizationDomain` concept. 

\sa `UniqueFactorizationDomain` 
\sa `AlgebraicStructureTraits_::Gcd` 

*/
template <class NT1, class NT2> result_type 
gcd(const NT1& x, const NT2& y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `integral_division()` (a.k.a.\ exact division or division without remainder) 
maps ring elements \f$ (x,y)\f$ to ring element \f$ z\f$ such that \f$ x = yz\f$ if such a \f$ z\f$ 
exists (i.e.\ if \f$ x\f$ is divisible by \f$ y\f$). Otherwise the effect of invoking 
this operation is undefined. Since the ring represented is an integral domain, 
\f$ z\f$ is uniquely defined if it exists. 

In case the argument types `NT1` and `NT2` differ, 
the `result_type` is determined via `Coercion_traits`. 

Thus, the `result_type` is well defined if `NT1` and `NT2` 
are a model of `ExplicitInteroperable`. 

The actual `integral_division` is performed with the semantic of that type. 

The function is defined if `result_type` 
is a model of the `IntegralDomain` concept. 

\sa `IntegralDomain` 
\sa `AlgebraicStructureTraits_::IntegralDivision` 

*/
template <class NT1, class NT2> result_type 
integral_division(const NT1& x, const NT2& y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `inverse()` returns the inverse element with respect to multiplication. 

The function is defined if the argument type 
is a model of the `Field` concept. 

\pre \f$ x \neq0\f$.

\sa `Field` 
\sa `AlgebraicStructureTraits_::Inverse` 

*/
template <class NT> NT inverse(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `is_negative()` determines if a value is negative or not. 
The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

The `result_type` is convertible to `bool`. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::IsNegative` 

*/
result_type is_negative(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `is_one()` determines if a value is equal to 1 or not. 

The function is defined if the argument type 
is a model of the `IntegralDomainWithoutDivision` concept. 

The `result_type` is convertible to `bool`. 

\sa `IntegralDomainWithoutDivision` 
\sa `AlgebraicStructureTraits_::IsOne` 

*/
template <class NT> result_type is_one(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `is_positive()` determines if a value is positive or not. 
The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

The `result_type` is convertible to `bool`. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::IsPositive` 

*/
result_type is_positive(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

An ring element \f$ x\f$ is said to be a square iff there exists a ring element 
\f$ y\f$ such 
that \f$ x= y*y\f$. In case the ring is a `UniqueFactorizationDomain`, 
\f$ y\f$ is uniquely defined up to multiplication by units. 

The function `is_square()` is available if 
`Algebraic_structure_traits::Is_square` is not the `CGAL::Null_functor`. 

The `result_type` is convertible to `bool`. 

\sa `UniqueFactorizationDomain` 
\sa `AlgebraicStructureTraits_::IsSquare` 

*/
template <class NT> result_type is_square(const NT& x);

/*!
\ingroup PkgAlgebraicFoundations

An ring element \f$ x\f$ is said to be a square iff there exists a ring element 
\f$ y\f$ such 
that \f$ x= y*y\f$. In case the ring is a `UniqueFactorizationDomain`, 
\f$ y\f$ is uniquely defined up to multiplication by units. 

The function `is_square()` is available if 
`Algebraic_structure_traits::Is_square` is not the `CGAL::Null_functor`. 

The `result_type` is convertible to `bool`.

\sa `UniqueFactorizationDomain` 
\sa `AlgebraicStructureTraits_::IsSquare` 

*/
template <class NT> result_type is_square(const NT& x, NT& y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `is_zero()` determines if a value is equal to 0 or not. 

The function is defined if the argument type 
is a model of the `RealEmbeddable` or of 
the `IntegralDomainWithoutDivision` concept. 

The `result_type` is convertible to `bool`. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::IsZero` 
\sa `IntegralDomainWithoutDivision` 
\sa `AlgebraicStructureTraits_::IsZero` 
*/
template <class NT> result_type is_zero(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `kth_root()` returns the k-th root of a value. 

The function is defined if the second argument type 
is a model of the `FieldWithKthRoot` concept. 

\sa `FieldWithKthRoot` 
\sa `AlgebraicStructureTraits_::KthRoot` 

*/
template <class NT> NT kth_root(int k, const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `mod()` computes the remainder of division with remainder. 

In case the argument types `NT1` and `NT2` differ, 
the `result_type` is determined via `Coercion_traits`. 

Thus, the `result_type` is well defined if `NT1` and `NT2` 
are a model of `ExplicitInteroperable`. 

The actual `mod` is performed with the semantic of that type. 

The function is defined if `result_type` 
is a model of the `EuclideanRing` concept. 

\sa `EuclideanRing` 
\sa `AlgebraicStructureTraits_::DivMod` 
\sa `CGAL::div_mod()`
\sa `CGAL::div()`

*/
template< class NT1, class NT2> 
result_type 
mod(const NT1& x, const NT2& y);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

returns the k-th real root of the univariate polynomial, which is
defined by the iterator range, where begin refers to the constant
term.

The function `root_of()` computes a real root of a square-free univariate 
polynomial. 

The function is defined if the value type, `NT`, 
of the iterator range is a model of the `FieldWithRootOf` concept. 

\pre The polynomial is square-free.

\sa `FieldWithRootOf` 
\sa `AlgebraicStructureTraits_::RootOf` 

*/
template <class InputIterator> NT 
root_of(int k, InputIterator begin, InputIterator end);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `sign()` returns the sign of its argument. 

The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

The `result_type` is convertible to `CGAL::Sign`. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::Sgn` 

*/
template <class NT> result_type sign(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `simplify()` may simplify a given object. 

The function is defined if the argument type 
is a model of the `IntegralDomainWithoutDivision` concept. 

\sa `IntegralDomainWithoutDivision` 
\sa `AlgebraicStructureTraits_::Simplify` 

*/
template <class NT> void simplify(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `sqrt()` returns the square root of a value. 

The function is defined if the argument type 
is a model of the `FieldWithSqrt` concept. 

\sa `FieldWithSqrt` 
\sa `AlgebraicStructureTraits_::Sqrt` 

*/
template <class NT> NT sqrt(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `square()` returns the square of a number. 

The function is defined if the argument type 
is a model of the `IntegralDomainWithoutDivision` concept. 

\sa `IntegralDomainWithoutDivision` 
\sa `AlgebraicStructureTraits_::Square` 

*/
template <class NT> NT square(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `to_double()` returns a double approximation of a number.
Note that in general, the value returned is not guaranteed to be the same
when called several times on the same number. For example, if `NT` is a lazy
number type (such as an instance of `CGAL::Lazy_exact_nt`), the double approximation
returned might be affected by an exact computation internally triggered
(that might have improved the double approximation).

The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

Remark: In order to control the quality of approximation one has to resort to methods that are specific to NT. There are no general guarantees whatsoever. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::ToDouble` 

*/
template <class NT> double to_double(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The template function `to_interval()` computes for a given real embeddable 
number \f$ x\f$ a double interval containing \f$ x\f$. 
This interval is represented by a `std::pair<double,double>`. 
The function is defined if the argument type 
is a model of the `RealEmbeddable` concept. 

\sa `RealEmbeddable` 
\sa `RealEmbeddableTraits_::ToInterval` 

*/
template <class NT> 
std::pair<double,double> to_interval(const NT& x);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgAlgebraicFoundations

The function `unit_part()` computes the unit part of a given ring 
element. 

The function is defined if the argument type 
is a model of the `IntegralDomainWithoutDivision` concept. 

\sa `IntegralDomainWithoutDivision` 
\sa `AlgebraicStructureTraits_::UnitPart` 

*/
template <class NT> NT unit_part(const NT& x);

} /* namespace CGAL */

