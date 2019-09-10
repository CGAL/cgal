
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3TraitsClasses

The class `Periodic_3_Delaunay_triangulation_traits_3` is designed as a default traits class for the
class `Periodic_3_Delaunay_triangulation_3<Periodic_3DelaunayTriangulationTraits_3,TriangulationDataStructure_3>`.

\tparam Traits must be a model of the `DelaunayTriangulationTraits_3` concept.
\tparam Offset must be a model of the concept `Periodic_3Offset_3` and defaults to `Periodic_3_offset_3`.

If `Traits` is a `CGAL::Filtered_kernel` (detected when `Traits::Has_filtered_predicates` exists
and is `true`), this class automatically provides filtered predicates. Similarly, statically filtered predicates
will be used if the flag `Traits::Has_static_filters` exists and is `true`.
By default, this holds for `CGAL::Exact_predicates_inexact_constructions_kernel` and
`CGAL::Exact_predicates_exact_constructions_kernel`.

\cgalModels Periodic_3DelaunayTriangulationTraits_3
*/
template< typename Traits, typename Offset >
class Periodic_3_Delaunay_triangulation_traits_3
  : public Periodic_3_triangulation_traits_base_3 < Traits, Offset> {
public:
}; /* end Periodic_3_Delaunay_triangulation_traits_3 */
} /* end namespace CGAL */
