
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3TraitsClasses

The class `Periodic_3_Delaunay_triangulation_traits_3` is designed as a default traits class for the
class `Periodic_3_Delaunay_triangulation_3<Periodic_3DelaunayTriangulationTraits_3,TriangulationDataStructure_3>`.

\tparam Traits must be a model of the `DelaunayTriangulationTraits_3` concept.
\tparam Periodic_3Offset_3 must be a model of the concept `Periodic_3Offset_3` and defaults to `Periodic_3_offset_3`. 

Note that this template class is specialized for 
`CGAL::Filtered_kernel`, so that it automatically provides 
filtered predicates. This holds implicitly for 
`CGAL::Exact_predicates_inexact_constructions_kernel`, as it is an 
instantiation of `CGAL::Filtered_kernel`. 

\cgalModels Periodic_3DelaunayTriangulationTraits_3
*/
template< typename Traits, typename Periodic_3Offset_3 >
class Periodic_3_Delaunay_triangulation_traits_3 : public Traits {
public:
}; /* end Periodic_3_Delaunay_triangulation_traits_3 */
} /* end namespace CGAL */
