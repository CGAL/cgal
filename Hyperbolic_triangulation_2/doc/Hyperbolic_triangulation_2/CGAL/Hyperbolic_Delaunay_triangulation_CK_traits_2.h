// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicTriangulation2TraitsClasses

The class `Hyperbolic_Delaunay_triangulation_CK_traits_2` is designed as one of the 
default models for the traits concept `HyperbolicDelaunayTriangulationTraits_2`
offered by %CGAL. 

\tparam K must be a model of `CircularKernel`.

This class provides exact constructions and predicates. The default value for `K` is 
`CGAL::Exact_circular_kernel_2`, which guarantees exact constructions of Delaunay 
triangulations and dual objects when the input points have rational coordinates.

\sa `Hyperbolic_Delaunay_triangulation_traits_2`

\cgalModels `HyperbolicDelaunayTriangulationTraits_2`
*/

template < class K >
class Hyperbolic_Delaunay_triangulation_CK_traits_2  : public K {

public:

  /// \name Types
  /// @{
    /*!
    \cgalModifBegin
    \cgalModifEnd
    */
    typedef typename K::FT                          FT;
    typedef typename K::Point_2                     Hyperbolic_point_2;
    /*!
    \cgalModifBegin
    \cgalModifEnd
    */
    typedef typename K::Circular_arc_point_2        Hyperbolic_Voronoi_point_2;
    typedef typename K::Circular_arc_2              Circular_arc_2;
    typedef typename K::Line_arc_2                  Line_arc_2; 
    typedef boost::variant<Circular_arc_2, 
                           Line_arc_2>              Hyperbolic_segment_2;
    typedef typename K::Triangle_2                  Hyperbolic_triangle_2;
  /// @}

  
  /// \name Creation
  /// @{
    /*!
      %Default constructor
    */
    Hyperbolic_Delaunay_triangulation_CK_traits_2(); 
    
    /*!
      Copy constructor
    */
    Hyperbolic_Delaunay_triangulation_CK_traits_2(const Hyperbolic_Delaunay_triangulation_CK_traits_2 & other);
  /// @}
    
};

} //namespace CGAL 


