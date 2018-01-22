// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2TraitsClasses

The class `Periodic_2_triangulation_traits_2` is designed as a default
traits class for the class
`Periodic_2_triangulation_2<Periodic_2TriangulationTraits_2,TriangulationDataStructure_2>`.

The argument \tparam Traits must be a model of the
`TriangulationTraits_2` concept. The argument
\tparam Periodic_2Offset_2 must be a model of the concept
`Periodic_2Offset_2` and defaults to `Periodic_2_offset_2`.

Note that this template class is specialized for
`CGAL::Filtered_kernel`, so that it automatically provides filtered
predicates. This holds implicitly for
`CGAL::Exact_predicates_inexact_constructions_kernel`, as it is an
instantiation of `CGAL::Filtered_kernel`.

\cgalModels `Periodic_2TriangulationTraits_2` and
\cgalModels `Periodic_2DelaunayTriangulationTraits_2` if the template parameter `Traits` is a model of the
`DelaunayTriangulationTraits_2` concept.

*/
template< typename Traits, typename Periodic_2Offset_2 >
class Periodic_2_triangulation_traits_2 : public Traits
{
}; /* end Periodic_2_triangulation_traits_2 */
} /* end namespace CGAL */
