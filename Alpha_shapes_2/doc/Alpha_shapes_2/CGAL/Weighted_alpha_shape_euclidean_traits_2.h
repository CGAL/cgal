
namespace CGAL {

/*!
\ingroup PkgAlphaShape2

The class `Weighted_alpha_shape_euclidean_traits_2` is the default model for the concept 
`AlphaShapeTraits_2` for the regular version of Alpha Shapes. 

\tparam K must be a `Kernel`. 

\cgalModels `AlphaShapeTraits_2`

*/
template< typename K >
class Weighted_alpha_shape_euclidean_traits_2 :
    public Regular_triangulation_euclidean_traits_2<K, typename K::FT> {
public:

}; /* end Weighted_alpha_shape_euclidean_traits_2 */
} /* end namespace CGAL */
