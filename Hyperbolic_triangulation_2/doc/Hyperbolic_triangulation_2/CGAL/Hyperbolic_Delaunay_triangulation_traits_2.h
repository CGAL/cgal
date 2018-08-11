// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicTriangulation2TraitsClasses

The class `Hyperbolic_Delaunay_triangulation_traits_2` is designed as one of the default traits classes for the
class `Hyperbolic_Delaunay_triangulation_2<HyperbolicDelaunayTriangulationTraits_2,TriangulationDataStructure_2>`.

\tparam K must be a model of `Kernel`.

If `K` provides exact computations with square roots, then this class automatically provides exact constructions
and predicates. The default value for `K` is `CGAL::Cartesian<CORE::Expr>`.

\sa Hyperbolic_Delaunay_triangulation_CK_traits_2

\cgalModels HyperbolicDelaunayTriangulationTraits_2
*/

template< class K >
class Hyperbolic_Delaunay_triangulation_traits_2 : public K {

public:

  /// \name Types
  /// @{
  	
  	typedef typename K::Point_2								Voronoi_point_2;
  	typedef unspecified_type								Circular_arc_2; 
    typedef typename K::Segment_2                           Euclidean_segment_2;
  	typedef boost::variant< Circular_arc_2, 
                            Euclidean_segment_2 >	        Hyperbolic_segment_2;
  /// @}
  

  /// \name Creation
  /// @{
    /*!
      Default constructor
    */
    Hyperbolic_Delaunay_triangulation_traits_2();

    /*!
      Copy constructor
    */
    Hyperbolic_Delaunay_triangulation_traits_2(const Hyperbolic_Delaunay_triangulation_traits_2 & other);
  /// @}


}; 


} // namespace CGAL
