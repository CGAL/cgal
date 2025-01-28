// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2TraitsClasses

The class `Periodic_2_triangulation_traits_2` is designed as a default
traits class for the class `CGAL::Periodic_2_triangulation_2<Periodic_2TriangulationTraits_2,TriangulationDataStructure_2>`.

\tparam Traits must be a model of the `TriangulationTraits_2` concept.
\tparam Periodic_2Offset_2 must be a model of the concept
`Periodic_2Offset_2` and defaults to `CGAL::Periodic_2_offset_2`.

If `Traits` is a `CGAL::Filtered_kernel` (detected when `Traits::Has_filtered_predicates` exists
and is `true`), this class automatically provides filtered predicates. Similarly, statically filtered predicates
will be used if the flag `Traits::Has_static_filters` exists and is `true`.
By default, this holds for `CGAL::Exact_predicates_inexact_constructions_kernel` and
`CGAL::Exact_predicates_exact_constructions_kernel`.

\cgalModels{Periodic_2TriangulationTraits_2}

*/
template< typename Traits, typename Periodic_2Offset_2 >
class Periodic_2_triangulation_traits_2 : public Traits
{
}; /* end Periodic_2_triangulation_traits_2 */
} /* end namespace CGAL */
