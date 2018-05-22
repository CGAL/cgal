// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

namespace CGAL {

/*!
\ingroup PkgHyperbolicTriangulation2TraitsClasses

The class `Hyperbolic_Delaunay_triangulation_CK_traits_2` is designed as one of the default traits classes for the
class `Hyperbolic_Delaunay_triangulation_2<HyperbolicDelaunayTriangulationTraits_2,TriangulationDataStructure_2>`.

\tparam K must be a model of `CircularKernel`.

This class provides exact constructions and predicates. The default value for `K` is `CGAL::Exact_circular_kernel_2`.

\sa `Hyperbolic_Delaunay_triangulation_traits_2`

\cgalModels `HyperbolicDelaunayTriangulationTraits_2`
*/

template < class K >
class Hyperbolic_Delaunay_triangulation_CK_traits_2  : public K {

public:

  /// \name Types
  /// @{
    typedef typename K::FT                          FT;
    typedef typename K::Point_2                     Point_2;
    typedef typename K::Circular_arc_point_2        Voronoi_point_2;
    typedef typename K::Circular_arc_2              Circular_arc_2;
    typedef typename K::Line_arc_2                  Line_arc_2; 
    typedef boost::variant<Circular_arc_2, 
                           Line_arc_2>              Hyperbolic_segment_2;
    typedef typename K::Triangle_2                  Triangle_2;
  /// @}

  /// \name Predicate Types
  /// @{                   
    typedef typename K::Orientation_2               Orientation_2;
    typedef typename K::Side_of_oriented_circle_2   Side_of_oriented_circle_2;
    typedef unspecified_type                        Side_of_hyperbolic_triangle_2;
    typedef unspecified_type                        Is_hyperbolic;
  /// @}  

  /// \name Construction Types
  /// @{                   
    typedef unspecified_type                        Construct_hyperbolic_segment_2;
    typedef unspecified_type                        Construct_hyperbolic_circumcenter_2;
    typedef unspecified_type                        Construct_hyperbolic_bisector_2;
    typedef unspecified_type                        Construct_intersection_2;
  /// @}
  
  /// \name Creation
  /// @{
    /*!
      Default constructor
    */
    Hyperbolic_Delaunay_triangulation_CK_traits_2(); 
    
    /*!
      Copy constructor
    */
    Hyperbolic_Delaunay_triangulation_CK_traits_2(const Hyperbolic_Delaunay_triangulation_CK_traits_2 & other);
  /// @}
    
};

} //namespace CGAL 


