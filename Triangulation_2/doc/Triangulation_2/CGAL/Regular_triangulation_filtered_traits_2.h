
namespace CGAL {

/*!
\ingroup PkgTriangulation2TraitsClasses

The class `Regular_triangulation_filtered_traits_2` is designed as a traits class for the 
class `Regular_triangulation_2<RegularTriangulationTraits_2,TriangulationDataStructure_2>`. 
Its difference with `Regular_triangulation_euclidean_traits_2` is that it 
provides filtered predicates which are meant to be fast and exact. 

The first argument `FK` must be a model of the `Kernel` concept, and 
it is also restricted to be an instance of the `Filtered_kernel` template. 

\cgalModels `RegularTriangulationTraits_2`

\sa `CGAL::Regular_triangulation_euclidean_traits_2`

*/
template< typename FK >
class Regular_triangulation_filtered_traits_2 : public Regular_triangulation_euclidean_traits_2<FK> {
public:

/// @}

}; /* end Regular_triangulation_filtered_traits_2 */
} /* end namespace CGAL */
