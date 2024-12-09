
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableFunctor` computes both integral quotient and remainder
of division with remainder. The quotient \f$ q\f$ and remainder \f$ r\f$ are computed
such that \f$ x = q*y + r\f$ and \f$ |r| < |y|\f$ with respect to the proper integer norm of
the represented ring.
\cgalFootnote{For integers this norm is the absolute value.
For univariate polynomials this norm is the degree.}
In particular, \f$ r\f$ is chosen to be \f$ 0\f$ if possible.
Moreover, we require \f$ q\f$ to be minimized with respect to the proper integer norm.

Note that the last condition is needed to ensure a unique computation of the
pair \f$ (q,r)\f$. However, an other option is to require minimality for \f$ |r|\f$,
with the advantage that
a <I>mod(x,y)</I> operation would return the unique representative of the
residue class of \f$ x\f$ with respect to \f$ y\f$, e.g. \f$ mod(2,3)\f$ should return \f$ -1\f$.
But this conflicts with nearly all current implementation
of integer types. From there, we decided to stay conform with common
implementations and require \f$ q\f$ to be computed as \f$ x/y\f$ rounded towards zero.

The following table illustrates the behavior for integers:

<TABLE CELLSPACING=5 >
<TR>
<TD ALIGN=CENTER NOWRAP>
<TABLE CELLSPACING=5 >
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
<TR>
<TD ALIGN=CENTER NOWRAP>
x
<TD ALIGN=CENTER NOWRAP>
y
<TD ALIGN=CENTER NOWRAP>
q
<TD ALIGN=CENTER NOWRAP>
r
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
<TR>
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
1
<TD ALIGN=CENTER NOWRAP>
0
<TR>
<TD ALIGN=CENTER NOWRAP>
2
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
2
<TR>
<TD ALIGN=CENTER NOWRAP>
1
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
1
<TR>
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
0
<TR>
<TD ALIGN=CENTER NOWRAP>
-1
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
-1
<TR>
<TD ALIGN=CENTER NOWRAP>
-2
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
-2
<TR>
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
-1
<TD ALIGN=CENTER NOWRAP>
0
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
</TABLE>

<TD ALIGN=CENTER NOWRAP>
-
<TD ALIGN=CENTER NOWRAP>
<TABLE CELLSPACING=5 >
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
<TR>
<TD ALIGN=CENTER NOWRAP>
x
<TD ALIGN=CENTER NOWRAP>
y
<TD ALIGN=CENTER NOWRAP>
q
<TD ALIGN=CENTER NOWRAP>
r
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
<TR>
<TD ALIGN=CENTER NOWRAP>
3
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
-1
<TD ALIGN=CENTER NOWRAP>
0
<TR>
<TD ALIGN=CENTER NOWRAP>
2
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
2
<TR>
<TD ALIGN=CENTER NOWRAP>
1
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
1
<TR>
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
0
<TR>
<TD ALIGN=CENTER NOWRAP>
-1
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
-1
<TR>
<TD ALIGN=CENTER NOWRAP>
-2
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
0
<TD ALIGN=CENTER NOWRAP>
-2
<TR>
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
-3
<TD ALIGN=CENTER NOWRAP>
1
<TD ALIGN=CENTER NOWRAP>
0
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=4><HR>
</TABLE>

</TABLE>

\cgalRefines{AdaptableFunctor}

\sa `AlgebraicStructureTraits`
\sa `AlgebraicStructureTraits_::Mod`
\sa `AlgebraicStructureTraits_::Div`

*/

class DivMod {
public:

/// \name Types
/// @{

/*!
Is void.
*/
typedef unspecified_type result_type;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type first_argument_type;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type second_argument_type;

/*!
Is `AlgebraicStructureTraits::Type&`.
*/
typedef unspecified_type third_argument_type;

/*!
Is `AlgebraicStructureTraits::Type&`.
*/
typedef unspecified_type fourth_argument_type;

/// @}

/// \name Operations
/// @{

/*!

computes the quotient \f$ q\f$ and remainder \f$ r\f$, such that \f$ x = q*y + r\f$
and \f$ r\f$ minimal with respect to the Euclidean Norm on
`Type`.
*/
result_type operator()( first_argument_type x,
second_argument_type y,
third_argument_type q,
fourth_argument_type r);

/*!
This operator is defined if `NT1` and `NT2` are `ExplicitInteroperable`
with coercion type `AlgebraicStructureTraits::Type`.
*/
template <class NT1, class NT2> result_type
operator()(NT1 x, NT2 y, third_argument_type q, fourth_argument_type r);

/// @}

}; /* end DivMod */

} /* end of namespace AlgebraicStructureTraits_ */
