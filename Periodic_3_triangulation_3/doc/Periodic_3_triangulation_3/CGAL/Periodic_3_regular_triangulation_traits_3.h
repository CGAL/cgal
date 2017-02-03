
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3TraitsClasses

The class `Periodic_3_regular_triangulation_traits_3` is designed as a default traits class for the
class `Periodic_3_regular_triangulation_3<Periodic_3RegularTriangulationTraits_3,TriangulationDataStructure_3>`.

\tparam K must be a model of the `RegularTriangulationTraits_3` concept.
\tparam Periodic_3Offset_3 must be a model of the concept `Periodic_3Offset_3` and defaults to `Periodic_3_offset_3`. 

\cgalModels Periodic_3RegularTriangulationTraits_3

This template class is specialized for 
`CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Filtered_kernel>`, 
so that it automatically provides 
filtered predicates. This holds implicitly for 
`CGAL::Regular_triangulation_euclidean_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel>`, 
as `CGAL::Exact_predicates_inexact_constructions_kernel` is an
instantiation of `CGAL::Filtered_kernel`. 

*/
template< typename K, typename Periodic_3Offset_3 >
class Periodic_3_regular_triangulation_traits_3 : 
  public Regular_triangulation_euclidean_traits_3< K > {
public:
}; /* end Periodic_3_regular_triangulation_traits_3 */
} /* end namespace CGAL */
