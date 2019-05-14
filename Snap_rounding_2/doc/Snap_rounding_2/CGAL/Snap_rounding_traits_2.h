
namespace CGAL {

/*!
\ingroup PkgSnapRounding2

The class `Snap_rounding_traits_2<Kernel>` is a model of the 
`SnapRoundingTraits_2` concept, and is the only traits class supplied 
with the package. 
This class should be instantiated with an exact geometric kernel that conforms 
to the \cgal kernel-concept, such as the 
`Exact_predicates_exact_constructions_kernel`, or `Cartesian<Gmpq>`. 

This geometric kernel must provide an (arbitrary-precision) rational number type 
(`FT`), `Point_2`, `Segment_2` and `Iso_rectangle_2`. 
It should be possible to cast numbers of the number type `FT` to 
double-precision representation. That is, the function 
`CGAL::to_double(FT)` must be supported. 

The `CGAL::to_double()` function is used to implement the operation that 
rounds the coordinates of a point to a center of a pixel. This operation is one 
of the traits-concept requirement. The coordinates are converted to double, 
rounded down to the nearest grid point, and finally adjusted to lie on a center 
of a pixel. Notice that if 
`CGAL::to_double()` returns the closet value, then when it rounds up a given 
coordinate, the resulting ISR, may be imprecise, and the distance between some 
vertex and a non-incident edge can be slightly less than the guaranteed 
half-the width-of-a-pixel. 

\cgalModels `SnapRoundingTraits_2`

*/
template< typename Kernel >
class Snap_rounding_traits_2 {
public:

/// @}

}; /* end Snap_rounding_traits_2 */
} /* end namespace CGAL */
