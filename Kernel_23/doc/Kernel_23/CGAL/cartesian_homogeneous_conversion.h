namespace CGAL {

/*!
\addtogroup kernel_conversion Cartesian/Homogenous Conversion
\ingroup PkgKernel23

Functions to convert between Cartesian and Homogeneous Kernels.

\sa `CGAL::Cartesian<FieldNumberType>` 
\sa `CGAL::Cartesian_converter<K1, K2, NTConverter>` 
\sa `CGAL::Homogeneous<RingNumberType>` 
\sa `CGAL::Homogeneous_converter<K1, K2, RTConverter, FTConverter>` 
\sa `CGAL::Simple_cartesian<FieldNumberType>` 
\sa `CGAL::Simple_homogeneous<RingNumberType>` 
*/
/// @{

/*!

converts 2d point `cp` with Cartesian representation 
into a 2d point with homogeneous representation with the same
number type.


*/
Point_2< Homogeneous<RT> >
cartesian_to_homogeneous(const Point_2< Cartesian<RT> >& cp);

/*!

converts 3d point `cp` with Cartesian representation 
into a 3d point with homogeneous representation with the same
number type.
*/
Point_3< Homogeneous<RT> >
cartesian_to_homogeneous(const Point_3< Cartesian<RT> >& cp);


/*!

converts 2d point `hp` with homogeneous representation 
into a 2d point with Cartesian representation with the same
number type.

*/
Point_2< Cartesian<FT> >
homogeneous_to_cartesian(const Point_2< Homogeneous<FT> >& hp);

/*!

converts 3d point `hp` with homogeneous representation 
into a 3d point with Cartesian representation with the same
number type.
*/
Point_3< Cartesian<FT> >
homogeneous_to_cartesian(const Point_3< Homogeneous<FT> >& hp);


/*!

converts the 2d point `hp` with homogeneous representation 
with number type `RT` into a 2d point with Cartesian 
representation with number type `Quotient<RT>`.
*/
Point_2< Cartesian<Quotient<RT> > >
homogeneous_to_quotient_cartesian(const Point_2<Homogeneous<RT> >& hp);

/*!

converts the 3d point `hp` with homogeneous representation 
with number type `RT` into a 3d point with Cartesian 
representation with number type `Quotient<RT>`.
*/
Point_3< Cartesian<Quotient<RT> > >
homogeneous_to_quotient_cartesian(const Point_3<Homogeneous<RT> >& hp);


/*!

converts 2d point `cp` with Cartesian representation 
with number type `Quotient<RT>` into a 2d point 
with homogeneous representation with number type `RT`.
*/
Point_2< Homogeneous<RT> >
quotient_cartesian_to_homogeneous(
const Point_2< Cartesian< Quotient<RT> > >& cp);

/*!

converts 3d point `cp` with Cartesian representation 
with number type `Quotient<RT>` into a 3d point 
with homogeneous representation with number type `RT`.

*/
Point_3< Homogeneous<RT> >
quotient_cartesian_to_homogeneous(
const Point_3< Cartesian< Quotient<RT> > >& cp);

} /* namespace CGAL */

