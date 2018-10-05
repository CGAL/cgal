// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2TraitsClasses

The class `Periodic_4_hyperbolic_Delaunay_triangulation_traits_2` is 
a model of the concept `Periodic_4HyperbolicDelaunayTriangulationTraits_2`.

It is the default traits class for the class `Periodic_4_hyperbolic_Delaunay_triangulation_2`.

\cgalModels `Periodic_4HyperbolicDelaunayTriangulationTraits_2`

\sa `Periodic_4_hyperbolic_Delaunay_triangulation_2`
*/
template< class Kernel = CGAL::Cartesian<CORE::Expr>>
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 
{
public:
	/// \name Types
	/// @{
		typedef typename Kernel::FT     				FT;
		typedef Hyperbolic_octagon_translation<FT> 		Hyperbolic_translation;
	/// @}
}; /* end Periodic_4_hyperbolic_Delaunay_triangulation_traits_2 */
} /* end namespace CGAL */
