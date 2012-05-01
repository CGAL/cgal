// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>
//                 Eric Berberich    <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_H

/*! \file
 * The topology traits for great spherical arcs embedded on a sphere for the
 * arrangement package.
 */

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_overlay_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_overlay_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_overlay_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_batched_pl_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_vert_decomp_sl_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>

#include <CGAL/Arr_topology_traits/Arr_spherical_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_spherical_insertion_helper.h>
#include <CGAL/Arr_topology_traits/Arr_spherical_overlay_helper.h>
#include <CGAL/Arr_topology_traits/Arr_spherical_batched_pl_helper.h>
#include <CGAL/Arr_topology_traits/Arr_spherical_vert_decomp_helper.h>
#include <CGAL/Arr_topology_traits/Arr_inc_insertion_zone_visitor.h>

#include <map>

namespace CGAL {

// Forward declaration:
template <class GeomTraits, class TopTraits> 
class Arrangement_on_surface_2;

/*! This class handles the topology for arrangements of great spherical
 * arcs on the sphere embedded on 2D parametric surdace.
 */
template <class GeomTraits, class T_Dcel = Arr_default_dcel<GeomTraits> >
class Arr_spherical_topology_traits_2 {
public:

  ///! \name The geometry-traits types.
  //@{
  typedef GeomTraits                                      Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
  //@}

  ///! \name The DCEL types.
  //@{
  typedef T_Dcel                                          Dcel;
  typedef typename Dcel::Size                             Size;
  typedef typename Dcel::Vertex                           Vertex;
  typedef typename Dcel::Halfedge                         Halfedge;
  typedef typename Dcel::Face                             Face;
  typedef typename Dcel::Outer_ccb                        Outer_ccb;
  typedef typename Dcel::Inner_ccb                        Inner_ccb;
  typedef typename Dcel::Isolated_vertex                  Isolated_vertex;
  //@}

  //! \name Arrangement types
  //!@{
  typedef Arr_spherical_topology_traits_2<Geometry_traits_2, Dcel> Self;
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>   Traits_adaptor_2;
  //!@}
  
  ///! \name The side tags
  //@{
  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;
  
  BOOST_MPL_ASSERT
  (
   (boost::mpl::or_< 
    boost::is_same< Left_side_category, Arr_oblivious_side_tag >,
    boost::is_same< Left_side_category, Arr_identified_side_tag > >)
  );
  BOOST_MPL_ASSERT
  (
   (boost::mpl::or_< 
    boost::is_same< Bottom_side_category, Arr_oblivious_side_tag >,
    boost::is_same< Bottom_side_category, Arr_contracted_side_tag > >)
  );
  BOOST_MPL_ASSERT
  (
   (boost::mpl::or_< 
    boost::is_same< Top_side_category, Arr_oblivious_side_tag >,
    boost::is_same< Top_side_category, Arr_contracted_side_tag > >)
  );
  BOOST_MPL_ASSERT
  (
   (boost::mpl::or_< 
    boost::is_same< Right_side_category, Arr_oblivious_side_tag >,
    boost::is_same< Right_side_category, Arr_identified_side_tag > >)
  );
  //@}
  
  /*! \struct
   * An auxiliary structure for rebinding the topology traits with a new 
   * geometry-traits class and a new DCEL class.
   */
  template<typename T, typename D>
  struct rebind
  {
    typedef Arr_spherical_topology_traits_2<T,D> other;
  };
  
private:

  //! A container of boundary vertices.
  struct Vertex_key_comparer {
    /*! Construct default */
    Vertex_key_comparer() : m_traits(NULL) {}
    
    /*! Construct */    
    Vertex_key_comparer(const Traits_adaptor_2 * traits) : m_traits(traits) {}
    const Traits_adaptor_2 * m_traits;
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return (m_traits->compare_y_on_boundary_2_object()(p1, p2) ==
              SMALLER);
    }
  };

  //! \todo define the key to be 'const Point_2*'.
  typedef std::map<Point_2,Vertex*,Vertex_key_comparer>  Vertex_map;
  typedef std::pair<Point_2,Vertex*>                     Vertex_value;

protected:
  // Data members:
  //! The DCEL.
  Dcel m_dcel;

  //! The spherical (bounded) face.  
  Face * m_spherical_face;

  //! The north pole
  Vertex * m_north_pole;

  //! The north pole
  Vertex * m_south_pole;

  //! The vertices on the discontinuity arc
  Vertex_map m_boundary_vertices;

  //! The geometry-traits adaptor.
  const Traits_adaptor_2 * m_traits;

  // Inidicates whether the traits object should evetually be freed.
  bool m_own_traits;
  
  // Copy constructor and assignment operator - not supported.
  Arr_spherical_topology_traits_2(const Self &);
  Self & operator=(const Self &);

public:
  ///! \name Construction methods.
  //@{

  /*! Default constructor. */
  Arr_spherical_topology_traits_2();

  /*! Constructor with a geometry-traits class.
   * \param traits the traits.
   */
  Arr_spherical_topology_traits_2(const Geometry_traits_2 * traits);

  /*! Assign the contents of another topology-traits class.
   * \param other the other spherical topology-traits.
   */
  void assign(const Self & other);
  //@}

  ///! \name Topology-traits methods.
  //@{

  /*! Obtain the DCEL (const version). */
  const Dcel & dcel() const
  {
    return (m_dcel);
  }

  /*! Obtain the DCEL (non-const version). */
  Dcel & dcel()
  {
    return m_dcel;
  }

  /*! Determine whether the DCEL reprsenets an empty structure.
   * \return true if the dcel reprsenets an empty structure; false otherwise.
   */
  bool is_empty_dcel() const
  {
    return (m_dcel.size_of_vertices() == 0);
  }

  /*! Initialize an empty DCEL structure. */
  void init_dcel();

  /*! Make the necessary updates after the DCEL structure have been updated. */
  void dcel_updated();

  /*! Determine whether the given vertex is concrete.
   * \param v the vertex.
   * \return true if v is mapped to a point on the discontinuity arc; false
   * otherwise.
   */
  bool is_concrete_vertex(const Vertex * /* v */) const
  {
    return true;
  }

  /*! Obtain the number of concrete vertices.
   * \return the number of concrete vertices.
   */
  Size number_of_concrete_vertices() const
  {
    return (m_dcel.size_of_vertices());
  }

  /*! Determine whether the given vertex is valid.
   * \param v the vertex.
   * \todo why is this needed, and where used?
   */
  bool is_valid_vertex (const Vertex * /* v */) const
  {
    return true;
  }

  /*! Obtain the number of valid vertices. */
  Size number_of_valid_vertices() const
  {
    return (m_dcel.size_of_vertices());
  }

  /*! Determine whether the given halfedge is valid. */
  bool is_valid_halfedge (const Halfedge * /* he */) const
  {
    return true;
  }

  /*! Obtain the number of valid halfedges. */
  Size number_of_valid_halfedges() const
  {
    return (m_dcel.size_of_halfedges());
  }

  /*! Determine whether the given face is valid. */
  bool is_valid_face(const Face * /* f */) const
  {
    return true;
  }
  
  /*! Obtain the number of valid faces. */
  Size number_of_valid_faces() const
  {
    return m_dcel.size_of_faces();
  }

  /*! Obtain the spherical face (const version). */
  const Face * spherical_face() const
  {
    return m_spherical_face;
  }

  /*! Obtain the spherical face (non-const version). */
  Face * spherical_face()
  {
    return m_spherical_face;
  }

  /*! Obtain the face containing the south pole (const version). */
  const Face * south_face() const
  {
    if (m_boundary_vertices.empty()) return m_spherical_face;
    typename Vertex_map::const_iterator it = m_boundary_vertices.begin();
    return _face_below_vertex_on_discontinuity(it->second);
  }

  /*! Obtain the face containing the south pole (non-const version). */
  Face * south_face()
  {
    if (m_boundary_vertices.empty()) return m_spherical_face;
    typename Vertex_map::iterator it = m_boundary_vertices.begin();
    return _face_below_vertex_on_discontinuity(it->second);
  }

  /*! Obtain the south pole (const version). */
  const Vertex * south_pole() const
  {
    return m_south_pole;
  }

  /*! Obtain the south pole (non-const version). */
  Vertex * south_pole()
  {
    return m_south_pole;
  }

  /*! Obtain the north pole (const version). */
  const Vertex * north_pole() const
  {
    return m_north_pole;
  }

  /*! Obtain the north pole (non-const version). */
  Vertex * north_pole()
  {
    return m_north_pole;
  }

  /*! Obtain a vertex on the line of discontinuity that corresponds to
   *  the given curve-end (or return NULL if no such vertex exists).
   */ 
  Vertex * discontinuity_vertex (const X_monotone_curve_2 xc, Arr_curve_end ind)
  {
    Point_2 key;
    if (ind == ARR_MIN_END)
      key = m_traits->construct_min_vertex_2_object() (xc);
    else
      key = m_traits->construct_max_vertex_2_object() (xc);

    typename Vertex_map::iterator  it = m_boundary_vertices.find(key);
    return (it != m_boundary_vertices.end()) ? it->second : NULL;
  }

  //@}
  
private:
  /// \name Auxiliary type definitions.
  //@{
  typedef Arrangement_on_surface_2<Geometry_traits_2, Self>        Arr;

  // Type definition for the constuction sweep-line visitor.
  typedef Arr_construction_subcurve<Geometry_traits_2>          CSubcurve; 
  typedef Arr_construction_event<Geometry_traits_2,CSubcurve,Arr>
                                                                CEvent;
  typedef Arr_spherical_construction_helper<Geometry_traits_2,Arr,
    CEvent,CSubcurve>                                           CHelper;

  // Type definition for the basic insertion sweep-line visitor.
  typedef Arr_basic_insertion_traits_2<Geometry_traits_2, Arr>  BInsTraits;
  typedef Arr_construction_subcurve<BInsTraits>                 BISubcurve; 
  typedef Arr_construction_event<BInsTraits,BISubcurve,Arr>     BIEvent;
  typedef Arr_spherical_insertion_helper<BInsTraits,Arr,BIEvent,BISubcurve>
                                                                BIHelper;

  // Type definition for the insertion sweep-line visitor.
  typedef Arr_insertion_traits_2<Geometry_traits_2, Arr>        InsTraits;
  typedef Arr_construction_subcurve<InsTraits>                  ISubcurve; 
  typedef Arr_construction_event<InsTraits,ISubcurve,Arr>       IEvent;
  typedef Arr_spherical_insertion_helper<InsTraits,Arr,IEvent,ISubcurve>
                                                                IHelper;

  // Type definition for the batched point-location sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              BplTraits;
  typedef Arr_spherical_batched_pl_helper<BplTraits, Arr>       BplHelper;

  // Type definition for the vertical decomposition sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              VdTraits;
  typedef Arr_spherical_vert_decomp_helper<VdTraits, Arr>       VdHelper;

  // Type definition for the overlay sweep-line visitor.
  template <class ExGeomTraits_, class ArrangementA_, class ArrangementB_>
  struct _Overlay_helper :
    public Arr_spherical_overlay_helper<ExGeomTraits_, ArrangementA_,
    ArrangementB_, Arr, Arr_construction_event<ExGeomTraits_,
    Arr_overlay_subcurve<ExGeomTraits_>, Arr>,
    Arr_overlay_subcurve<ExGeomTraits_> >
  {
    typedef Arr_spherical_overlay_helper<ExGeomTraits_, ArrangementA_,
      ArrangementB_, Arr, Arr_construction_event<ExGeomTraits_,
      Arr_overlay_subcurve<ExGeomTraits_>, Arr>,
      Arr_overlay_subcurve<ExGeomTraits_> >             Base;

    typedef typename Base::Traits_2                     Traits_2;
    typedef typename Base::Arrangement_red_2            Arrangement_red_2;
    typedef typename Base::Arrangement_blue_2           Arrangement_blue_2;
    typedef typename Base::Arrangement_2                Arrangement_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;
    typedef typename Base::Construction_helper          Construction_helper;

    _Overlay_helper(const ArrangementA_ * arrA, const ArrangementB_ * arrB) :
      Base(arrA, arrB) {}
  };
  //@}

public:
  ///! \name Visitor types.
  //@{

  typedef Arr_construction_sl_visitor<CHelper>
                            Sweep_line_construction_visitor;

  typedef Arr_insertion_sl_visitor<IHelper>
                            Sweep_line_insertion_visitor;

  typedef Sweep_line_construction_visitor
                            Sweep_line_non_intersecting_construction_visitor;
  
  typedef Arr_basic_insertion_sl_visitor<BIHelper>
                            Sweep_line_non_intersecting_insertion_visitor;

  template <class OutputIterator_>
  struct Sweep_line_batched_point_location_visitor :
    public Arr_batched_pl_sl_visitor<BplHelper, OutputIterator_>
  {
    typedef OutputIterator_                             Output_iterator;
    typedef Arr_batched_pl_sl_visitor<BplHelper,Output_iterator>  Base;

    typedef typename Base::Traits_2                     Traits_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;

    Sweep_line_batched_point_location_visitor(const Arr * arr,
                                              Output_iterator * oi) :
      Base(arr, oi)
    {}
  };

  template <class OutputIterator_>
  struct Sweep_line_vertical_decomposition_visitor :
    public Arr_vert_decomp_sl_visitor<VdHelper, OutputIterator_>
  {
    typedef OutputIterator_                             Output_iterator;
    typedef Arr_vert_decomp_sl_visitor<VdHelper,Output_iterator>  Base;

    typedef typename Base::Traits_2                     Traits_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;

    Sweep_line_vertical_decomposition_visitor(const Arr * arr,
                                              Output_iterator * oi) :
      Base(arr, oi)
    {}
  };

  template <class ArrangementA_, class ArrangementB_, class OverlayTraits_>
  struct Sweep_line_overlay_visitor :
    public Arr_overlay_sl_visitor
  <_Overlay_helper<Arr_overlay_traits_2<Geometry_traits_2,ArrangementA_,
    ArrangementB_>,
    ArrangementA_, ArrangementB_>, OverlayTraits_>
  {
    typedef ArrangementA_                               ArrangementA_2;
    typedef ArrangementB_                               ArrangementB_2;
    typedef Arr                                         Arrangement_result_2;
    typedef OverlayTraits_                              Overlay_traits;

    typedef Arr_overlay_traits_2<Geometry_traits_2,ArrangementA_2,
      ArrangementB_2>                                   Geom_ovl_traits_2;

    typedef _Overlay_helper<Geom_ovl_traits_2,ArrangementA_2,ArrangementB_2>
                                                        Ovl_helper;

    typedef Arr_overlay_sl_visitor<Ovl_helper,Overlay_traits>   Base;

    typedef typename Base::Traits_2                     Traits_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;

    Sweep_line_overlay_visitor(const ArrangementA_2 * arrA,
                               const ArrangementB_2 * arrB,
                               Arrangement_result_2 * arr_res,
                               Overlay_traits * overlay_tr) :
      Base(arrA, arrB, arr_res, overlay_tr)
    {}
  };

  typedef Arr_inc_insertion_zone_visitor<Arr> Zone_insertion_visitor;

  typedef Arr_naive_point_location<Arr>       Default_point_location_strategy;
  //@}

  ///! \name Topology-traits methods.
  //@{

  /*! Receive a notification on the creation of a new boundary vertex that
   * corresponds to the given curve end.
   * \param v The new boundary vertex.
   * \param xc The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   */
  void notify_on_boundary_vertex_creation(Vertex * v,
                                          const X_monotone_curve_2 & xc,
                                          Arr_curve_end ind,
                                          Arr_parameter_space ps_x,
                                          Arr_parameter_space ps_y);

  /*! Given two predecessor halfedges that belong to the same inner CCB of
   * a face, determine what happens when we insert an edge connecting the
   * target vertices of the two edges.
   * \param prev1 The first predecessor halfedge.
   * \param prev2 The second predecessor halfedge.
   * \param cv The curve to be inserted
   * \pre The two halfedges belong to the same inner CCB.
   * \return A pair indicating whether the insertion will cause the face
   *         to split (the first flag), and if so - whether the split face
   *         will form a hole in the original face.
   */
  std::pair<bool, bool>
  face_split_after_edge_insertion(const Halfedge * prev1,
                                  const Halfedge * prev2,
                                  const X_monotone_curve_2 & /*cv*/) const
  {
    CGAL_precondition(prev1->is_on_inner_ccb());
    CGAL_precondition(prev2->is_on_inner_ccb());
    CGAL_precondition(prev1->inner_ccb() == prev2->inner_ccb());

    // In case of a planar topology, connecting two vertices on the same
    // inner CCB closes a new face that becomes a hole in the original face:
    return (std::make_pair(true, true));
  }

  /*! Determine whether the removal of the given edge will cause the creation
   * of a hole.
   * \param he The halfedge to be removed.
   * \pre Both he and its twin lie on an outer CCB of their incident faces.
   * \return Whether a new hole will be created.
   */
  bool hole_creation_after_edge_removal(const Halfedge * he) const
  {
    CGAL_precondition(!he->is_on_inner_ccb());
    CGAL_precondition(!he->opposite()->is_on_inner_ccb());

    /* Check whether the halfedge and its twin belong to the same outer CCB
     * (and are therefore incident to the same face).
     * If they do, cut an antenna. That is, seperate a new hole from the outer
     * CCB of the face (return true).
     * Otherwise, the edge separates two faces. When removed, these two faces
     * will be merged, but no new hole will be created (return false).
     */
    return (he->outer_ccb() == he->opposite()->outer_ccb());
  }

  /*! Given two predecessor halfedges that will be used for inserting a
   * new halfedge pair (prev1 will be the predecessor of the halfedge he1,
   * and prev2 will be the predecessor of its twin he2), such that the
   * insertion will create a new perimetric face that forms a hole inside
   * an existing perimetric face, determine whether he1 will be incident to
   * this new face.
   * \param prev1 The first predecessor halfedge.
   * \param prev2 The second predecessor halfedge.
   * \param xc The x-monotone curve we use to connect prev1's target and
   *           prev2's target vertex.
   * \param try_other_way Used to temporarily solve a bug in case that if 
   *                      we go from prev1 to prev2 we get a perimetric face
   *                      but from prev2 to prev1 we don't get a perimetric 
   *                      face. The flag is used to try the other direction.
   *                      
   * \pre prev1 and prev2 belong to the same inner connected component.
   * \return true if he1 (and prev1) lies in the interior of the face we
   *         are about to create, false otherwise - in which case he2
   *         (and prev2) must be incident to this new face.
   */
  bool is_on_new_perimetric_face_boundary(const Halfedge * prev1,
                                          const Halfedge * prev2,
                                          const X_monotone_curve_2 & xc,
                                          bool try_other_way = true) const;

  /*! Determine whether the two halfedges, belonging to different outer CCBs,
   * belong to the outer boundary of the same face.
   * \param e1 The first halfedge.
   * \param e2 The second halfedge.
   * \return Whether the two halfedge belong to the outer boundary of the same
   *         face.
   */
  bool boundaries_of_same_face(const Halfedge * /* e1 */,
                               const Halfedge * /* e2 */) const
  {
    // This function is never called in case of an arrangement on a sphere:
    CGAL_error();
    return false;
  }

  /*! Determine whether a given point lies in the interior of a given face.
   * \param f The face.
   * \param p The query point.
   * \param v The vertex associated with p (if exists).
   * \param f must not be fictitious, and v must not lie at infinity.
   * \return true if p is contained in f's interior; false otherwise.
   */
  bool is_in_face(const Face * f, const Point_2 & p, const Vertex * v) const;

  /*! Compare the relative y-position of a given point and a given edge.
   * \param p The point.
   * \param he The edge (one of the pair of halfedges).
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  Comparison_result compare_y_at_x(const Point_2 & p,
                                   const Halfedge * he) const;

  /*! Determine whether a given vertex is associated with a given curve end.
   * \param v The vertex.
   * \param xc The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return Whether v represents the given curve end.
   */
  bool are_equal(const Vertex * v,
                 const X_monotone_curve_2 & xc, Arr_curve_end ind,
                 Arr_parameter_space ps_x, Arr_parameter_space ps_y) const;

  /*! Given a curve end with boundary conditions and a face that contains the
   * interior of the curve, find a place for a boundary vertex that will
   * represent the curve end along the face boundary.
   * \param f The face.
   * \param xc The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return An object that contains the curve end.
   */
  CGAL::Object place_boundary_vertex(Face * f,
                                     const X_monotone_curve_2 & xc,
                                     Arr_curve_end ind,
                                     Arr_parameter_space ps_x,
                                     Arr_parameter_space ps_y);

  /*!
   * Locate the predecessor halfedge for the given curve around a given
   * vertex with boundary conditions.
   * \param v The vertex.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y, and should be
   *      incident to the vertex v.
   * \return An object that contains the curve end.
   */
  Halfedge * locate_around_boundary_vertex(Vertex * v,
                                           const X_monotone_curve_2 & cv,
                                           Arr_curve_end ind,
                                           Arr_parameter_space ps_x,
                                           Arr_parameter_space ps_y) const;

  /*! Locate a DCEL feature that contains the given curve end.
   * \param xc The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve end is incident to the boundary.
   * \return An object that contains the curve end.
   */
  CGAL::Object locate_curve_end(const X_monotone_curve_2 & xc, Arr_curve_end ce,
                                Arr_parameter_space ps_x,
                                Arr_parameter_space ps_y);

  /*! Split a fictitious edge using the given vertex.
   * \param e The edge to split (one of the pair of halfedges).
   * \param v The split vertex.
   * \pre e is a fictitious halfedge.
   * \return A halfedge whose direction is the same as e's and whose target is
   *         the split vertex v.
   */
  Halfedge * split_fictitious_edge(Halfedge * /* e */, Vertex * /* v */)
  {
    // There are no fictitious halfedges:
    CGAL_error();
    return NULL;
  }

  /*! Determine whether the given face is unbounded.
   * \param f The face.
   * \return true if f is unbounded; false otherwise.
   */
  bool is_unbounded(const Face * /* f */) const
  {
    // All faces on a sphere are bounded:
    return false;
  }

  /*! Determine whether the given boundary vertex is redundant.
   * \param v The vertex.
   * \return Whether v is redundant, and should be erased.
   */
  bool is_redundant(const Vertex * v) const;

  /*! Erase the given redundant vertex by merging a fictitious edge.
   * The function does not free the vertex v itself.
   * \param v The vertex.
   * \pre v is a redundant vertex.
   * \return One of the pair of halfedges that form the merged edge.
   */
  Halfedge * erase_redundant_vertex(Vertex * v);

  //! reference_face (const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common 
      point.
      \return A pointer to the reference face.
  */
  const Face* reference_face() const
  {
    return spherical_face();
  }
  
  //! reference_face (non-const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common 
      point.
      \return A pointer to the reference face.
  */
  Face* reference_face()
  {
    return spherical_face();
  }
  
  //@}

protected:

  /// \name Auxiliary functions.
  //@{
  
  /*! Obtain the curve associated with a boundary vertex.
   * \param v The boundary vertex.
   * \param ind Output: ARR_MIN_END if the vertex is induced by the minimal end;
   *                    ARR_MAX_END if it is induced by the curve's maximal end.
   * \pre v is a valid boundary.
   * \return The curve that induces v.
   */
  const X_monotone_curve_2 & _curve(const Vertex * v,
                                    Arr_curve_end & ind) const;

  /*! Return the halfedge, the target vertex of which is given, that is
   * the predecessor of a halfedge, the curve of which is given, that is about
   * to be inserted into the dcel.
   */
  Halfedge *
  _locate_around_vertex_on_discontinuity(Vertex * v,
                                         const X_monotone_curve_2 & xc,
                                         Arr_curve_end ind) const;

  /*! Return the halfedge, the target vertex of which is a given pole,
   * that is the predecessor of a halfedge, the curve of which is given, that
   * is about to be inserted into the dcel.
   */
  Halfedge * _locate_around_pole(Vertex * v, const X_monotone_curve_2 & xc,
                                 Arr_curve_end ind) const;


  /*! Return the face that lies below the given vertex, which lies
   * on the line of discontinuity.
   */
  Face * _face_below_vertex_on_discontinuity(Vertex * v) const;

  //@}
};

} //namespace CGAL

#include <CGAL/Arr_topology_traits/Arr_spherical_topology_traits_2_impl.h>

#endif
