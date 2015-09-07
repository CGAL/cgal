// Copyright (c) 2006,2007,2008,2009,2010,2011,2012,2013 Tel-Aviv University (Israel).
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
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                 Efi Fogel  <efif@post.tau.ac.il>

#ifndef CGAL_ARR_UNB_PLANAR_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_UNB_PLANAR_TOPOLOGY_TRAITS_2_H

/*! \file
 * Definition of the Arr_unb_planar_topology_traits_2<GeomTraits> class.
 */

#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_topology_traits/Arr_planar_topology_traits_base_2.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_insertion_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_overlay_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_batched_pl_helper.h>
#include <CGAL/Arr_topology_traits/Arr_unb_planar_vert_decomp_helper.h>
#include <CGAL/Arr_topology_traits/Arr_inc_insertion_zone_visitor.h>
#include <CGAL/assertions.h>

namespace CGAL {

// Forward declaration:
template <typename GeomTraits_, typename TopTraits_>
class Arrangement_on_surface_2;

/*! \class Arr_unb_planar_topology_traits_2
 * A topology-traits class that encapsulates the embedding of 2D arrangements
 * of unbounded curves on the plane.
 */
template <typename GeomTraits_,
          typename Dcel_ = Arr_default_dcel<GeomTraits_> >
class Arr_unb_planar_topology_traits_2 :
  public Arr_planar_topology_traits_base_2<GeomTraits_, Dcel_>
{
private:
  typedef Arr_planar_topology_traits_base_2<GeomTraits_, Dcel_>
                                                          Base;

public:
  ///! \name The geometry-traits types.
  //@{
  typedef GeomTraits_                                     Geometry_traits_2;
  typedef typename Base::Point_2                          Point_2;
  typedef typename Base::X_monotone_curve_2               X_monotone_curve_2;
  //@}

  ///! \name The DCEL types.
  //@{
  typedef Dcel_                                           Dcel;
  typedef typename Base::Size                             Size;
  typedef typename Base::Vertex                           Vertex;
  typedef typename Base::Halfedge                         Halfedge;
  typedef typename Base::Face                             Face;
  typedef typename Base::Outer_ccb                        Outer_ccb;
  typedef typename Base::Inner_ccb                        Inner_ccb;
  typedef typename Base::Isolated_vertex                  Isolated_vertex;
  //@}

  //! \name Arrangement types
  //!@{
  typedef Arr_unb_planar_topology_traits_2<Geometry_traits_2, Dcel> Self;
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>   Traits_adaptor_2;
  //!@}

  ///! \name The side tags
  //@{
  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  BOOST_MPL_ASSERT(
      (boost::mpl::or_<
       boost::is_same< Left_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Left_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_<
       boost::is_same< Bottom_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Bottom_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_<
       boost::is_same< Top_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Top_side_category, Arr_open_side_tag > >
      )
  );
  BOOST_MPL_ASSERT(
      (boost::mpl::or_<
       boost::is_same< Right_side_category, Arr_oblivious_side_tag >,
       boost::is_same< Right_side_category, Arr_open_side_tag > >
      )
  );
  //@}

  /*! \struct
   * An auxiliary structure for rebinding the topology traits with a new
   * geometry-traits class and a new DCEL class.
   */
  template <typename T, typename D>
  struct rebind {
    typedef Arr_unb_planar_topology_traits_2<T, D> other;
  };

protected:
  // Data members:
  Vertex* v_bl;         // A fictitious vertex at (-oo,-oo).
  Vertex* v_tl;         // A fictitious vertex at (-oo,+oo).
  Vertex* v_br;         // A fictitious vertex at (-oo,+oo).
  Vertex* v_tr;         // A fictitious vertex at (+oo,+oo).
  Size n_inf_verts;     // Number of vertices at infinity.
  Face* fict_face;      // The fictitious DCEL face.

  // Copy constructor and assignment operator - not supported.
  Arr_unb_planar_topology_traits_2(const Self&);
  Self& operator=(const Self&);

public:
  ///! \name Construction methods.
  //@{

  /*! Default constructor. */
  Arr_unb_planar_topology_traits_2();

  /*! Constructor with a geometry-traits class. */
  Arr_unb_planar_topology_traits_2(const Geometry_traits_2* tr);

  /*! Assign the contents of another topology-traits class. */
  void assign(const Self& other);
  //@}

  ///! \name Accessing the DCEL and constructing iterators.
  //@{

  /*! Determine whether the DCEL reprsenets an empty structure. */
  bool is_empty_dcel() const
  {
    // An empty arrangement contains just two four vertices at infinity
    // and eight fictitious halfedges connecting them.
    return (this->m_dcel.size_of_vertices() == 4 &&
            this->m_dcel.size_of_halfedges() == 8);
  }

  /*! Check if the given vertex is concrete (associated with a point). */
  bool is_concrete_vertex(const Vertex* v) const
  {
    return (! v->has_null_point());
  }

  /*! Get the number of concrete vertices. */
  Size number_of_concrete_vertices() const
  {
    // All vertices not lying at infinity are concrete.
    return (this->m_dcel.size_of_vertices() - n_inf_verts);
  }

  /*! Check if the given vertex is valid (not a fictitious one). */
  bool is_valid_vertex(const Vertex* v) const
  {
    return (! v->has_null_point() ||
            (v != v_bl && v != v_tl && v != v_br && v != v_tr));
  }

  /*! Get the number of valid vertices. */
  Size number_of_valid_vertices() const
  {
    // All vertices, except the four fictitious one, are valid.
    return (this->m_dcel.size_of_vertices() - 4);
  }

  /*! Check if the given halfedge is valid (not a fictitious one). */
  bool is_valid_halfedge(const Halfedge* he) const
  {
    return (! he->has_null_curve());
  }

  /*! Get the number of valid halfedges. */
  Size number_of_valid_halfedges() const
  {
    // Note that we do not count fictitious halfedges (each vertex at infinity
    // induces two fictitious halfedges).
    return (this->m_dcel.size_of_halfedges() - 2*n_inf_verts);
  }

  /*! Check if the given face is valid (not a fictitious one). */
  bool is_valid_face (const Face* f) const
  {
    return (! f->is_fictitious());
  }

  /*! Get the number of valid faces. */
  Size number_of_valid_faces() const
  {
    // We do not count the ficitious DCEL face.
    return (this->m_dcel.size_of_faces() - 1);
  }
  //@}

private:
  /// \name Auxiliary type definitions.
  //@{
  typedef Arrangement_on_surface_2<Geometry_traits_2, Self>    Arr;

  // Type definition for the constuction sweep-line visitor.
  typedef Arr_construction_subcurve<Geometry_traits_2>         CSubcurve;
  typedef Arr_construction_event<Geometry_traits_2, CSubcurve, Arr>
                                                               CEvent;
  typedef Arr_unb_planar_construction_helper<Geometry_traits_2,
                                             Arr,
                                             CEvent,
                                             CSubcurve>        CHelper;

  // Type definition for the basic insertion sweep-line visitor.
  typedef Arr_basic_insertion_traits_2<Geometry_traits_2, Arr> BInsTraits;
  typedef Arr_construction_subcurve<BInsTraits>                BISubcurve;
  typedef Arr_construction_event<BInsTraits, BISubcurve, Arr>
                                                               BIEvent;
  typedef Arr_unb_planar_insertion_helper<BInsTraits, Arr, BIEvent, BISubcurve>
                                                               BIHelper;

  // Type definition for the insertion sweep-line visitor.
  typedef Arr_insertion_traits_2<Geometry_traits_2, Arr>       InsTraits;
  typedef Arr_construction_subcurve<InsTraits>                 ISubcurve;
  typedef Arr_construction_event<InsTraits, ISubcurve, Arr>
                                                               IEvent;
  typedef Arr_unb_planar_insertion_helper<InsTraits, Arr, IEvent, ISubcurve>
                                                               IHelper;

  // Type definition for the batched point-location sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>             BplTraits;
  typedef Arr_unb_planar_batched_pl_helper<BplTraits, Arr>     BplHelper;

  // Type definition for the vertical decomposition sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>             VdTraits;
  typedef Arr_unb_planar_vert_decomp_helper<VdTraits, Arr>     VdHelper;

  // Type definition for the overlay sweep-line visitor.
  template <class ExGeomTraits_, class ArrangementA_, class ArrangementB_>
  struct _Overlay_helper : public Arr_unb_planar_overlay_helper
      <ExGeomTraits_, ArrangementA_, ArrangementB_, Arr,
       Arr_construction_event<ExGeomTraits_,
                              Arr_overlay_subcurve<ExGeomTraits_>,
                              Arr>,
       Arr_overlay_subcurve<ExGeomTraits_> >
  {
    typedef Arr_unb_planar_overlay_helper
              <ExGeomTraits_, ArrangementA_, ArrangementB_, Arr,
               Arr_construction_event<ExGeomTraits_,
                                      Arr_overlay_subcurve<ExGeomTraits_>,
                                      Arr>,
               Arr_overlay_subcurve<ExGeomTraits_> >     Base;

    typedef typename Base::Traits_2                      Traits_2;
    typedef typename Base::Arrangement_red_2             Arrangement_red_2;
    typedef typename Base::Arrangement_blue_2            Arrangement_blue_2;
    typedef typename Base::Arrangement_2                 Arrangement_2;
    typedef typename Base::Event                         Event;
    typedef typename Base::Subcurve                      Subcurve;
    typedef typename Base::Construction_helper           Construction_helper;

    _Overlay_helper(const ArrangementA_* arrA, const ArrangementB_* arrB) :
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
    typedef OutputIterator_                                   Output_iterator;

    typedef Arr_batched_pl_sl_visitor<BplHelper, Output_iterator>   Base;
    typedef typename Base::Traits_2                                 Traits_2;
    typedef typename Base::Event                                    Event;
    typedef typename Base::Subcurve                                 Subcurve;

    Sweep_line_batched_point_location_visitor(const Arr* arr,
                                              Output_iterator& oi) :
      Base(arr, oi)
    {}
  };

  template <class OutputIterator_>
  struct Sweep_line_vertical_decomposition_visitor :
    public Arr_vert_decomp_sl_visitor<VdHelper, OutputIterator_>
  {
    typedef OutputIterator_                                   Output_iterator;
    typedef Arr_vert_decomp_sl_visitor<VdHelper, Output_iterator>
                                                             Base;

    typedef typename Base::Traits_2                           Traits_2;
    typedef typename Base::Event                              Event;
    typedef typename Base::Subcurve                           Subcurve;

    Sweep_line_vertical_decomposition_visitor(const Arr* arr,
                                              Output_iterator* oi) :
      Base(arr, oi) {}
  };

  template <class ArrangementA_, class ArrangementB_, class OverlayTraits_>
  struct Sweep_line_overlay_visitor :
    public Arr_overlay_sl_visitor
        <_Overlay_helper<Arr_overlay_traits_2<Geometry_traits_2,
                                              ArrangementA_,
                                              ArrangementB_>,
                         ArrangementA_,
                         ArrangementB_>,
         OverlayTraits_>
  {
    typedef ArrangementA_                            ArrangementA_2;
    typedef ArrangementB_                            ArrangementB_2;
    typedef Arr                                      Arrangement_result_2;
    typedef OverlayTraits_                           Overlay_traits;

    typedef Arr_overlay_traits_2<Geometry_traits_2,
                                 ArrangementA_2,
                                 ArrangementB_2>     Geom_ovl_traits_2;

    typedef _Overlay_helper<Geom_ovl_traits_2,
                            ArrangementA_2,
                            ArrangementB_2>          Ovl_helper;

    typedef Arr_overlay_sl_visitor<Ovl_helper,
                                   Overlay_traits>   Base;

    typedef typename Base::Traits_2                  Traits_2;
    typedef typename Base::Event                     Event;
    typedef typename Base::Subcurve                  Subcurve;

    Sweep_line_overlay_visitor(const ArrangementA_2* arrA,
                               const ArrangementB_2* arrB,
                               Arrangement_result_2* arr_res,
                               Overlay_traits* overlay_tr) :
      Base(arrA, arrB, arr_res, overlay_tr)
    {}
  };

  typedef Arr_inc_insertion_zone_visitor<Arr>
                                        Zone_insertion_visitor;

  typedef Arr_walk_along_line_point_location<Arr>
                                        Default_point_location_strategy;

  typedef Arr_walk_along_line_point_location<Arr>
                                        Default_vertical_ray_shooting_strategy;
  //@}

  ///! \name Topology-traits methods.
  //@{

  /*!
   * Initialize an empty DCEL structure.
   */
  void init_dcel();

  /*!
   * Make the necessary updates after the DCEL structure have been updated.
   */
  void dcel_updated();

  /*!
   * Check if the given vertex is associated with the given curve end.
   * \param v The vertex.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return Whether v represents the given curve end.
   */
  bool are_equal(const Vertex* v,
                 const X_monotone_curve_2& cv, Arr_curve_end ind,
                 Arr_parameter_space ps_x, Arr_parameter_space ps_y) const;

  /*!
   * Given a curve end with boundary conditions and a face that contains the
   * interior of the curve, find a place for a boundary vertex that will
   * represent the curve end along the face boundary.
   * \param f The face.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return An object that contains the curve end.
   *         In our case this object always wraps a fictitious edge.
   */
  CGAL::Object place_boundary_vertex(Face* f,
                                     const X_monotone_curve_2& cv,
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
  Halfedge*
  locate_around_boundary_vertex(Vertex* /* v */,
                                const X_monotone_curve_2& /* cv */,
                                Arr_curve_end /* ind */,
                                Arr_parameter_space /* ps_x */,
                                Arr_parameter_space /* ps_y */) const
  {
    CGAL_error();
    return (NULL);
  }

  /*!
   * Locate a DCEL feature that contains the given unbounded curve end.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve end is unbounded in either x or y.
   * \return An object that contains the curve end.
   *         In our case this object may either wrap an unbounded face,
   *         or an edge with an end-vertex at infinity (in case of an overlap).
   */
  CGAL::Object locate_curve_end(const X_monotone_curve_2& cv,
                                Arr_curve_end ind,
                                Arr_parameter_space ps_x,
                                Arr_parameter_space ps_y);

  /*!
   * Split a fictitious edge using the given vertex.
   * \param e The edge to split (one of the pair of halfedges).
   * \param v The split vertex.
   * \pre e is a fictitious halfedge.
   * \return A halfedge whose direction is the same as e's and whose target is
   *         the split vertex v.
   */
  Halfedge* split_fictitious_edge(Halfedge* e, Vertex* v);

  /*!
   * Determine whether the given face is unbounded.
   * \param f The face.
   * \return Whether f is unbounded.
   */
  bool is_unbounded(const Face* f) const;

  /*!
   * Determine whether the given boundary vertex is redundant.
   * \param v The vertex.
   * \return Whether v is redundant, and should be erased.
   */
  bool is_redundant(const Vertex* v) const;

  /*!
   * Erase the given redundant vertex by merging a fictitious edge.
   * The function does not free the vertex v itself.
   * \param v The vertex.
   * \pre v is a redundant vertex.
   * \return One of the pair of halfedges that form the merged edge.
   */
  Halfedge* erase_redundant_vertex(Vertex* v);

    //! reference_face (const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common
      point.
      \return A pointer to the reference face.
  */
  const Face* reference_face() const
  {
    CGAL_assertion(v_tr->halfedge()->direction() == ARR_LEFT_TO_RIGHT);
    return v_tr->halfedge()->outer_ccb()->face();
  }

  //! reference_face (non-const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common
      point.
      \return A pointer to the reference face.
  */
  Face* reference_face()
  {
    CGAL_assertion(v_tr->halfedge()->direction() == ARR_LEFT_TO_RIGHT);
    return v_tr->halfedge()->outer_ccb()->face();
  }

  //@}

  /// \name Additional accessors, specialized for this topology-traits class.
  //@{

  /*! This function is used by the "walk" point-location strategy. */
  const Face* initial_face() const { return fict_face; }

  /*! Get the fictitious face (const version). */
  const Face* fictitious_face() const { return fict_face; }

  /*! Get the fictitious face (non-const version). */
  Face* fictitious_face() { return fict_face; }

  /*! Get the bottom-left fictitious vertex (const version). */
  const Vertex* bottom_left_vertex() const { return (v_bl); }

  /*! Get the bottom-left fictitious vertex (non-const version). */
  Vertex* bottom_left_vertex() { return (v_bl); }

  /*! Get the top-left fictitious vertex (const version). */
  const Vertex* top_left_vertex() const { return (v_tl); }

  /*! Get the top-left fictitious vertex (non-const version). */
  Vertex* top_left_vertex() { return (v_tl); }

  /*! Get the bottom-right fictitious vertex (const version). */
  const Vertex* bottom_right_vertex() const { return (v_br); }

  /*! Get the bottom-right fictitious vertex (non-const version). */
  Vertex* bottom_right_vertex() { return (v_br); }

  /*! Get the top-right fictitious vertex (const version). */
  const Vertex* top_right_vertex() const { return (v_tr); }

  /*! Get the top-right fictitious vertex (non-const version). */
  Vertex* top_right_vertex() { return (v_tr); }
  //@}

  /// \name Additional predicates, specialized for this topology-traits class.
  //@{

  /*!
   * Compare the given vertex (which may lie at infinity) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the comparison of the x-coordinates of p and v.
   */
  virtual Comparison_result compare_x(const Point_2& p,
                                      const Vertex* v) const;

  /*!
   * Compare the given vertex (which may lie at infinity) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the xy-lexicographic comparison of p and v.
   */
  virtual Comparison_result compare_xy(const Point_2& p,
                                       const Vertex* v) const;

  /*!
   * Compare the relative y-position of the given point and the given edge
   * (which may be fictitious).
   * \param p The point.
   * \param he The edge (one of the pair of halfedges).
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  virtual Comparison_result compare_y_at_x(const Point_2& p,
                                           const Halfedge* he) const;
  //@}

protected:

  /// \name Auxiliary functions.
  //@{

  /*!
   * Get the curve associated with a boundary vertex.
   * \param v The vertex as infinity.
   * \param ind Output: ARR_MIN_END if the vertex is induced by the minimal end;
   *                    ARR_MAX_END if it is induced by the curve's maximal end.
   * \pre v is a valid (not fictitious) boundary.
   * \return The curve that induces v, or NULL if v has no incident curves yet.
   */
  const X_monotone_curve_2* _curve(const Vertex* v, Arr_curve_end& ind) const;

  /*!
   * Check whether the given infinite curve end lies on the given fictitious
   * halfedge.
   * \param cv The curve.
   * \param ind Whether we refer to the minimal or maximal end of cv.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \param he The fictitious halfedge.
   * \param eq_source Output: Whether the curve coincides with he's source.
   * \param eq_target Output: Whether the curve coincides with he's target.
   * \return Whether the curve end lies on the fictitious halfedge.
   */
  bool _is_on_fictitious_edge(const X_monotone_curve_2& cv, Arr_curve_end ind,
                               Arr_parameter_space ps_x,
                               Arr_parameter_space ps_y,
                               const Halfedge* he,
                               bool& eq_source, bool& eq_target);
  //@}
};

} //namespace CGAL

#include <CGAL/Arr_topology_traits/Arr_unb_planar_topology_traits_2_impl.h>

#endif
