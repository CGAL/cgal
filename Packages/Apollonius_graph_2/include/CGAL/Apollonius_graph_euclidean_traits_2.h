#ifndef CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Hyperbola_2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>

#include <CGAL/Apollonius_graph_method_tags.h>


#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/predicates/Apollonius_graph_ftC2.h>
#include <CGAL/Apollonius_graph_constructions_ftC2.h>
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates/Apollonius_graph_rtH2.h>
#include <CGAL/Apollonius_graph_constructions_rtH2.h>
#endif



#include <CGAL/Apollonius_graph_euclidean_traits_2.C>



CGAL_BEGIN_NAMESPACE


//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
template < class R, class W, class Method_tag >
class Apollonius_graph_euclidean_traits_2
  : private Triangulation_euclidean_traits_2< R >
{
private:
  typedef Triangulation_euclidean_traits_2<R>    Base_traits;

public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  typedef R                                            Rep;
  typedef W                                            Weight;
  typedef typename R::Point_2                          Bare_point;
  typedef CGAL::Weighted_point<Bare_point, W>          Weighted_point;
  // I must add the following for some reason, otherwise it uses the
  // definition of Point_2 from the Triangulation_euclidean_traits_2
  typedef Weighted_point                               Point_2;
  //
  typedef Weighted_point                               Site;
  typedef typename R::Line_2                           Line_2;
  typedef typename R::Ray_2                            Ray_2;
  typedef typename R::Segment_2                        Segment_2;
  typedef Parabola_segment_2< Bare_point, W, Line_2 >  Parabola_segment_2;
  typedef Hyperbola_2< Bare_point, W >                 Hyperbola_2;
  typedef Hyperbola_ray_2< Bare_point, W >             Hyperbola_ray_2;
  typedef Hyperbola_segment_2< Bare_point, W >         Hyperbola_segment_2;

private:
  typedef Bare_point                             Pt;
  typedef Method_tag                             MTag;

public:
  // CONSTRUCTIONS
  //--------------
  // vertex and weighted point
  typedef CGAL::Construct_Apollonius_vertex_2<R>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_weighted_point_2<Pt,W,Line_2>
  /*                              */ Construct_Apollonius_weighted_point_2;

  // bisectors and subsets
  typedef CGAL::Construct_Apollonius_bisector_2<Pt,W,Line_2>
  /*                                    */ Construct_Apollonius_bisector_2;
  typedef CGAL::Construct_Apollonius_bisector_ray_2<Pt,W,Line_2,Ray_2>
  /*                                */ Construct_Apollonius_bisector_ray_2;
  typedef CGAL::Construct_Apollonius_bisector_segment_2<Pt,W,Segment_2>
  /*                            */ Construct_Apollonius_bisector_segment_2;

  // primal edges
  typedef CGAL::Construct_Apollonius_primal_ray_2<Pt,W,Line_2,Ray_2> 
  /*                                   */ Construct_Apollonius_primal_ray_2;
  typedef CGAL::Construct_Apollonius_primal_segment_2<Pt,W,Line_2,Segment_2>
  /*                               */ Construct_Apollonius_primal_segment_2;


  // PREDICATES
  //-----------
  typedef typename Base_traits::Compare_x_2             Compare_x_2;
  typedef typename Base_traits::Compare_y_2             Compare_y_2;
  typedef CGAL::Compare_weight_2<Pt,W>                  Compare_weight_2;
  typedef typename Base_traits::Orientation_2           Orientation_2;
  typedef CGAL::Is_trivial_2<Pt,W,MTag>                 Is_trivial_2;
  typedef CGAL::Oriented_side_of_bisector_2<Pt,W,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<Pt,W,MTag >           Vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<Pt,W,MTag >
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<Pt,W,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<Pt,W,MTag>         Is_degenerate_edge_2;



public:
  //-----------------------------------------------------------------------
  //                  ACCESS TO OBJECTS
  //-----------------------------------------------------------------------

  // CONSTRUCTIONS
  //--------------
  inline Construct_Apollonius_vertex_2
  construct_Apollonius_vertex_2_object() const { 
    return Construct_Apollonius_vertex_2();
  }

  inline Construct_Apollonius_weighted_point_2
  construct_Apollonius_weighted_point_2_object() const {
    return Construct_Apollonius_weighted_point_2();
  }

  inline Construct_Apollonius_bisector_2
  construct_Apollonius_bisector_2_object() const {
    return Construct_Apollonius_bisector_2();
  }

  inline Construct_Apollonius_bisector_ray_2
  construct_Apollonius_bisector_ray_2_object() const {
    return Construct_Apollonius_bisector_ray_2();
  }

  inline Construct_Apollonius_bisector_segment_2
  construct_Apollonius_bisector_segment_2_object() const { 
    return Construct_Apollonius_bisector_segment_2(); 
  }

  inline Construct_Apollonius_primal_ray_2
  construct_Apollonius_primal_ray_2_object() const {
    return Construct_Apollonius_primal_ray_2(); 
  }

  inline Construct_Apollonius_primal_segment_2
  construct_Apollonius_primal_segment_2_object() const { 
    return Construct_Apollonius_primal_segment_2();
  }

  // PREDICATES
  //-----------
  inline Compare_x_2
  compare_x_2_object() const {
    return Compare_x_2();
  }

  inline Compare_y_2
  compare_y_2_object() const {
    return Compare_y_2();
  }

  inline Compare_weight_2
  compare_weight_2_object() const {
    return Compare_weight_2();
  }

  inline Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  inline Is_trivial_2
  is_trivial_2_object() const {
    return Is_trivial_2();
  }

  inline Oriented_side_of_bisector_2
  oriented_side_of_bisector_2_object() const {
    return Oriented_side_of_bisector_2();
  }

  inline Vertex_conflict_2
  vertex_conflict_2_object() const {
    return Vertex_conflict_2();
  }

  inline Finite_edge_interior_conflict_2
  finite_edge_interior_conflict_2_object() const {
    return Finite_edge_interior_conflict_2();
  }

  inline Infinite_edge_interior_conflict_2
  infinite_edge_interior_conflict_2_object() const {
    return Infinite_edge_interior_conflict_2();
  }

  inline Is_degenerate_edge_2
  is_degenerate_edge_2_object() const {
    return Is_degenerate_test_2();
  }

};


CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H
