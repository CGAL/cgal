#ifndef CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Hyperbola_2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Filtered_predicate.h>

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


#include <CGAL/Number_type_traits.h>

#include <CGAL/Apollonius_graph_euclidean_traits_2.C>

#include <CGAL/Kernel_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
template < class R, class MTag = Ring_tag >
class Apollonius_graph_euclidean_traits_2
{
public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  
  typedef Kernel_wrapper_2<R>                           Kernel;
  typedef R                                             Rep;
  typedef MTag                                          Method_tag;
  typedef typename Kernel::RT                           Weight;
  typedef typename Kernel::Point_2                      Bare_point;
  typedef typename Kernel::Weighted_point_2             Weighted_point;
  // I must add the following for some reason, otherwise it uses the
  // definition of Point_2 from the Triangulation_euclidean_traits_2
  typedef Weighted_point                                Point_2;
  //
  typedef Weighted_point                                Site;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef Parabola_segment_2<Bare_point,Weight,Line_2>  Parabola_segment_2;
  typedef Hyperbola_2<Bare_point,Weight>                Hyperbola_2;
  typedef Hyperbola_ray_2<Bare_point,Weight>            Hyperbola_ray_2;
  typedef Hyperbola_segment_2<Bare_point,Weight>        Hyperbola_segment_2;

public:
  // CONSTRUCTIONS
  //--------------
  // vertex and weighted point
  typedef CGAL::Construct_Apollonius_vertex_2<R>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_weighted_point_2<R>
  /*                              */ Construct_Apollonius_weighted_point_2;

  // bisectors and subsets
  typedef CGAL::Construct_Apollonius_bisector_2<R>
  /*                                    */ Construct_Apollonius_bisector_2;
  typedef CGAL::Construct_Apollonius_bisector_ray_2<R>
  /*                                */ Construct_Apollonius_bisector_ray_2;
  typedef CGAL::Construct_Apollonius_bisector_segment_2<R>
  /*                            */ Construct_Apollonius_bisector_segment_2;

  // primal edges
  typedef CGAL::Construct_Apollonius_primal_ray_2<R> 
  /*                                   */ Construct_Apollonius_primal_ray_2;
  typedef CGAL::Construct_Apollonius_primal_segment_2<R>
  /*                               */ Construct_Apollonius_primal_segment_2;


  // PREDICATES
  //-----------
  typedef typename Kernel::Compare_x_2                  Compare_x_2;
  typedef typename Kernel::Compare_y_2                  Compare_y_2;
  typedef CGAL::Compare_weight_2<R>                     Compare_weight_2;
  typedef typename Kernel::Orientation_2                Orientation_2;
  typedef CGAL::Is_hidden_2<R,MTag>                     Is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<R,MTag> 
  /*                                          */ Oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<R,MTag >              Vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<R,MTag >
  /*                                      */ Finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<R,MTag>
  /*                                    */ Infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<R,MTag>            Is_degenerate_edge_2;


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

  inline Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
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
    return Is_degenerate_edge_2();
  }

};


#if 0
// I AM REMOVING THE FILTERED KERNEL STUFF SO THAT I DON'T GET ANY
// PROBLEMS WITH PARTIAL SPECIALIZATION AND ALSO BECAUSE THE CURRENT
// DEFINITION OF FILTERED KERNEL DOES NOT SUPPORT MORE THAN ONE
// TEMPLATE ARGUMENT

//-----------------------------------------------------------------------
// the Traits class for a filtered kernel
//-----------------------------------------------------------------------
template< class MTag, class _CK, class _EK, class _FK,
  class _C2E, class _C2F >
class Apollonius_graph_euclidean_traits_2<
Filtered_kernel<_CK,_EK,_FK,_C2E,_C2F>, MTag>
{
private:
  // add the wrappers; these should be removed once weighted points
  // find their way into the kernel
  typedef Kernel_wrapper_2<_CK>  CK;
  typedef Kernel_wrapper_2<_EK>  EK;
  typedef Kernel_wrapper_2<_FK>  FK;
  typedef Extended_cartesian_converter<_C2E>  C2E;
  typedef Extended_cartesian_converter<_C2F>  C2F;

public:
  //-----------------------------------------------------------------------
  //                  TYPE DEFINITIONS
  //-----------------------------------------------------------------------

  // BASIC TYPES
  //------------
  typedef Kernel_wrapper_2< Filtered_kernel<CK,EK,FK,C2E,C2F> >  Kernel;
  typedef Kernel                                        Rep;
  typedef MTag                                          Method_tag;
  typedef typename Kernel::RT                           Weight;
  typedef typename Kernel::Point_2                      Bare_point;
  typedef typename Kernel::Weighted_point_2             Weighted_point;
  // I must add the following for some reason, otherwise it uses the
  // definition of Point_2 from the Triangulation_euclidean_traits_2
  typedef Weighted_point                                Point_2;
  //
  typedef Weighted_point                                Site;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef Parabola_segment_2<Bare_point,Weight,Line_2>  Parabola_segment_2;
  typedef Hyperbola_2<Bare_point,Weight>                Hyperbola_2;
  typedef Hyperbola_ray_2<Bare_point,Weight>            Hyperbola_ray_2;
  typedef Hyperbola_segment_2<Bare_point,Weight>        Hyperbola_segment_2;

private:
  // Types for the construction kernel
  typedef typename CK::RT                                CK_weight;
  typedef typename CK::Point_2                           CK_bare_point;
  typedef CGAL::Weighted_point<CK_bare_point,CK_weight>  CK_weighted_point;

  // Types for the exact kernel
  typedef typename EK::RT                                EK_weight;
  typedef typename EK::Point_2                           EK_bare_point;
  typedef CGAL::Weighted_point<EK_bare_point,EK_weight>  EK_weighted_point;

  // Types for the filtering kernel
  typedef typename FK::RT                                FK_weight;
  typedef typename FK::Point_2                           FK_bare_point;
  typedef CGAL::Weighted_point<FK_bare_point,FK_weight>  FK_weighted_point;

public:
  // CONSTRUCTIONS
  //--------------
  // vertex and weighted point
  typedef CGAL::Construct_Apollonius_vertex_2<Rep>
  /*                                      */ Construct_Apollonius_vertex_2;
  typedef CGAL::Construct_Apollonius_weighted_point_2<Rep>
  /*                              */ Construct_Apollonius_weighted_point_2;

  // bisectors and subsets
  typedef CGAL::Construct_Apollonius_bisector_2<Rep>
  /*                                    */ Construct_Apollonius_bisector_2;
  typedef CGAL::Construct_Apollonius_bisector_ray_2<Rep>
  /*                                */ Construct_Apollonius_bisector_ray_2;
  typedef CGAL::Construct_Apollonius_bisector_segment_2<Rep>
  /*                            */ Construct_Apollonius_bisector_segment_2;

  // primal edges
  typedef CGAL::Construct_Apollonius_primal_ray_2<Rep>
  /*                                   */ Construct_Apollonius_primal_ray_2;
  typedef CGAL::Construct_Apollonius_primal_segment_2<Rep>
  /*                               */ Construct_Apollonius_primal_segment_2;


private:
  // Predicates for the construction kernel
  typedef CGAL::Compare_weight_2<EK>                   CK_compare_weight_2;
  typedef CGAL::Is_hidden_2<CK,MTag>                   CK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<CK,MTag> 
  /*                                      */ CK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<CK,MTag >            CK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<CK,MTag >
  /*                                  */ CK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<CK,MTag>
  /*                                */ CK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<CK,MTag>       CK_is_degenerate_edge_2;

  // Predicates for the exact kernel
  typedef CGAL::Compare_weight_2<EK>                   EK_compare_weight_2;
  typedef CGAL::Is_hidden_2<EK,MTag>                   EK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<EK,MTag> 
  /*                                      */ EK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<EK,MTag >            EK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<EK,MTag >
  /*                                  */ EK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<EK,MTag>
  /*                                */ EK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<EK,MTag>       EK_is_degenerate_edge_2;

  // Predicates for the filtering kernel
  typedef CGAL::Compare_weight_2<FK>                   FK_compare_weight_2;
  typedef CGAL::Is_hidden_2<FK,MTag>                   FK_is_hidden_2;
  typedef CGAL::Oriented_side_of_bisector_2<FK,MTag> 
  /*                                      */ FK_oriented_side_of_bisector_2;
  typedef CGAL::Vertex_conflict_2<FK,MTag >            FK_vertex_conflict_2;
  typedef CGAL::Finite_edge_interior_conflict_2<FK,MTag >
  /*                                  */ FK_finite_edge_interior_conflict_2;
  typedef CGAL::Infinite_edge_interior_conflict_2<FK,MTag>
  /*                                */ FK_infinite_edge_interior_conflict_2;
  typedef CGAL::Is_degenerate_edge_2<FK,MTag>       FK_is_degenerate_edge_2;

public:
  // PREDICATES
  //-----------
  typedef typename Kernel::Compare_x_2                   Compare_x_2;
  typedef typename Kernel::Compare_y_2                   Compare_y_2;
  typedef Filtered_predicate<EK_compare_weight_2,
    FK_compare_weight_2,C2E,C2F> Compare_weight_2;
  typedef typename Kernel::Orientation_2                 Orientation_2;
  typedef Filtered_predicate<EK_is_hidden_2,
    FK_is_hidden_2,C2E,C2F> Is_hidden_2;
  typedef Filtered_predicate<EK_oriented_side_of_bisector_2,
    FK_oriented_side_of_bisector_2,C2E,C2F> Oriented_side_of_bisector_2;
  typedef Filtered_predicate<EK_vertex_conflict_2,
    FK_vertex_conflict_2,C2E,C2F> Vertex_conflict_2;
  typedef Filtered_predicate<EK_finite_edge_interior_conflict_2,
    FK_finite_edge_interior_conflict_2,C2E,C2F>
  Finite_edge_interior_conflict_2;
  typedef Filtered_predicate<EK_infinite_edge_interior_conflict_2,
    FK_infinite_edge_interior_conflict_2,C2E,C2F>
  Infinite_edge_interior_conflict_2;
  typedef Filtered_predicate<EK_is_degenerate_edge_2,
    FK_is_degenerate_edge_2,C2E,C2F> Is_degenerate_edge_2;

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

  inline Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
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
    return Is_degenerate_edge_2();
  }

};
#endif


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_2_H
