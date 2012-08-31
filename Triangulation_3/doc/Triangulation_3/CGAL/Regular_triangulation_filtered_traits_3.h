
namespace CGAL {

/*!
\ingroup PkgTriangulation3TraitsClasses

\deprecated This class is deprecated since \cgal 3.6. 
The class `CGAL::Regular_triangulation_euclidean_traits_3` should be used instead. 
Filtered predicates are automatically used if the boolean `Has_filtered_predicates` 
in the kernel provided as template parameter of that class is set to `true`. 

The class `Regular_triangulation_filtered_traits_3` is designed as a traits class for the 
class `Regular_triangulation_3<RegularTriangulationTraits_3,TriangulationDataStructure_3>`. 
Its difference with `Regular_triangulation_euclidean_traits_3` is that it 
provides filtered predicates which are meant to be fast and exact. 

The first argument `FK` must be a model of the `Kernel` concept, and 
it is also restricted to be a instance of the `Filtered_kernel` template. 

\models ::RegularTriangulationTraits_3 

\sa `CGAL::Regular_triangulation_euclidean_traits_3`. 

*/
template< typename FK >
class Regular_triangulation_filtered_traits_3 : public Regular_triangulation_euclidean_traits_3<FK> {
public:

/// @}

}; /* end Regular_triangulation_filtered_traits_3 */
} /* end namespace CGAL */
