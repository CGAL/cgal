
namespace CGAL {

/*!
\ingroup PkgTriangulation2TraitsClasses

`Regular_triangulation_euclidean_traits_2` is a model for the concept `RegularTriangulationTraits_2` 
This traits class is templated by a kernel class `K` 
and a weight type `Weight`. 
This class inherits from `K` 
and uses a `Weighted_point` type 
derived from the type `K::Point_2`. 

Note that this template class is specialized for 
`Exact_predicates_inexact_constructions_kernel`, so that it is as if 
`Regular_triangulation_filtered_traits_2` was used, i.e.\ you get 
filtered predicates automatically. 

\cgalModels `RegularTriangulationTraits_2`

\sa `RegularTriangulationTraits_2` 
\sa `CGAL::Regular_triangulation_filtered_traits_2` 
\sa `CGAL::Regular_triangulation_2` 

*/
template< typename K, typename Weight >
class Regular_triangulation_euclidean_traits_2 : public K {
public:

/// @}

}; /* end Regular_triangulation_euclidean_traits_2 */
} /* end namespace CGAL */
