namespace CGAL {

/*!
\ingroup rational_rotation_approximation_grp

computes integers `sin_num`, `cos_num` and `denom`, such
that `sin_num`/`denom` approximates the sine of direction
\f$ (\f$`dirx`,`diry`\f$ )\f$. The difference between the sine and
the approximating rational is bounded by `eps_num`/`eps_den`.
\pre `eps_num` \f$ \neq0\f$.

\cgalHeading{Implementation}

The approximation is based on Farey sequences as described in
the rational rotation method presented by Canny and Ressler at the
8th SoCG 1992. We use a slower version which needs no division operation
in the approximation.

\sa `CGAL::Aff_transformation_2<Kernel>`
*/
template <RingNumberType>
void
rational_rotation_approximation( const RingNumberType & dirx,
                                 const RingNumberType & diry,
                                 RingNumberType & sin_num,
                                 RingNumberType & cos_num,
                                 RingNumberType & denom,
                                 const RingNumberType & eps_num,
                                 const RingNumberType & eps_den );

} /* namespace CGAL */

