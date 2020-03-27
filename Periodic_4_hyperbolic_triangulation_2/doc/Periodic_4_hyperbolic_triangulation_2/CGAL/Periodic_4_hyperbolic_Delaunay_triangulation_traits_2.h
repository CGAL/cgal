// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2TraitsClasses

The class `Periodic_4_hyperbolic_Delaunay_triangulation_traits_2` is
the default traits class for the class `Periodic_4_hyperbolic_Delaunay_triangulation_2`.

\cgalModels `Periodic_4HyperbolicDelaunayTriangulationTraits_2`

\sa `Periodic_4_hyperbolic_Delaunay_triangulation_2`
*/
template< class Kernel = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt >
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 :
        public Hyperbolic_Delaunay_triangulation_traits_2<Kernel>
{
public:
        /// \name Types
        /// @{

                typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>                                  Base;

                typedef typename Base::FT                                                                                                     FT;
                typedef Hyperbolic_octagon_translation<FT>                                                                         Hyperbolic_translation;

                typedef typename Base::Hyperbolic_point_2                                                                     Hyperbolic_point_2;
                typedef Hyperbolic_point_2                                                                                                         Hyperbolic_Voronoi_point_2;
                typedef typename Base::Hyperbolic_segment_2                                                                   Hyperbolic_segment_2;

        /// @}
}; /* end Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 */
} /* end namespace CGAL */
