
namespace CGAL {

/*!
\ingroup PkgTriangulation3TraitsClasses

\deprecated The class is deprecated since \cgal 4.10, as the weighted point and the function
objects for weighted points are part of the concept `Kernel`. The class is kept for backward
compatibility, but ignores the template parameter `Weight`.

\tparam K must be a model of the `Kernel` concept.

\tparam Weight This template parameter is ignored, as `Kernel::Weighted_point_3`
               uses the type `Kernel::FT`.

\cgalModels{RegularTriangulationTraits_3}

*/
template< typename K, typename Weight >
class Regular_triangulation_euclidean_traits_3 : public K {
public:

}; /* end Regular_triangulation_euclidean_traits_3 */

} /* end namespace CGAL */
