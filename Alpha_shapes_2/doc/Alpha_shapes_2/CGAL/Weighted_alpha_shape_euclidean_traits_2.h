
namespace CGAL {

/*!
\ingroup PkgAlphaShape2

\deprecated The class is deprecated since \cgal 4.10, as the weighted point and the function
objects for weighted points are part of the concept `Kernel`. The class is kept for backward
compatibility.

The class `Weighted_alpha_shape_euclidean_traits_2` was the default model for the concept
`AlphaShapeTraits_2` for the regular version of Alpha Shapes.

\tparam K must be a model of `Kernel`.

\cgalModels `AlphaShapeTraits_2`

*/
template< typename K >
class Weighted_alpha_shape_euclidean_traits_2
  : public K
{
public:

}; /* end Weighted_alpha_shape_euclidean_traits_2 */
} /* end namespace CGAL */
