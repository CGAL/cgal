// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s) : Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efif@post.tau.ac.il>
//             Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_bounded_planar_topology_traits_2<GeomTraits> class.
 */

#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_topology_traits/Arr_planar_topology_traits_base_2.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_insertion_helper.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_overlay_helper.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_batched_pl_helper.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_vert_decomp_helper.h>
#include <CGAL/Arr_topology_traits/Arr_inc_insertion_zone_visitor.h>
#include <CGAL/Surface_sweep_2/No_overlap_subcurve.h>
#include <CGAL/Surface_sweep_2/No_overlap_event.h>
#include <CGAL/use.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

// Forward declaration:
template <typename GeometryTraits_2, typename TopTraits_>
class Arrangement_on_surface_2;

/*! \class Arr_bounded_planar_topology_traits_2
 *
 * A topology-traits class that encapsulates the embedding of 2D arrangements
 * of bounded curves on the plane.
 */
template <typename GeometryTraits_2,
          typename Dcel_ = Arr_default_dcel<GeometryTraits_2> >
class Arr_bounded_planar_topology_traits_2 :
  public Arr_planar_topology_traits_base_2<GeometryTraits_2, Dcel_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Dcel_                                         Dcel;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_planar_topology_traits_base_2<Gt2, Dcel_> Base;

public:
  ///! \name The geometry-traits types.
  //@{
  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  //@}

  ///! \name The DCEL types.
  //@{
  typedef typename Base::Size                            Size;
  typedef typename Base::Vertex                          Vertex;
  typedef typename Base::Halfedge                        Halfedge;
  typedef typename Base::Face                            Face;
  typedef typename Base::Outer_ccb                       Outer_ccb;
  typedef typename Base::Inner_ccb                       Inner_ccb;
  typedef typename Base::Isolated_vertex                 Isolated_vertex;
  //@}

  //! \name Arrangement types
  //!@{
  typedef Arr_bounded_planar_topology_traits_2<Gt2, Dcel> Self;
  typedef Arr_traits_basic_adaptor_2<Gt2>                 Traits_adaptor_2;
  //!@}

  ///! \name The side tags
  //@{
  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  BOOST_MPL_ASSERT
  ((boost::is_same< Left_side_category, Arr_oblivious_side_tag >));
  BOOST_MPL_ASSERT
  ((boost::is_same< Bottom_side_category, Arr_oblivious_side_tag >));
  BOOST_MPL_ASSERT
  ((boost::is_same< Top_side_category, Arr_oblivious_side_tag >));
  BOOST_MPL_ASSERT
  ((boost::is_same< Right_side_category, Arr_oblivious_side_tag >));
  //@}

  /*! \struct
   * An auxiliary structure for rebinding the topology traits with a new
   * geometry-traits class and a new DCEL class.
   */
  template <typename T, typename D>
  struct rebind {
    typedef Arr_bounded_planar_topology_traits_2<T, D> other;
  };

protected:
  // Data members:
  Face* unb_face;     // The unbounded face.

  // Copy constructor and assignment operator - not supported.
  Arr_bounded_planar_topology_traits_2(const Self&);
  Self& operator=(const Self&);

public:
  ///! \name Construction methods.
  //@{

  /*! Default constructor. */
  Arr_bounded_planar_topology_traits_2() :
    Base(),
    unb_face(NULL)
  {}

  /*! Constructor from a geometry-traits object. */
  Arr_bounded_planar_topology_traits_2(const Gt2* traits) :
    Base(traits),
    unb_face(NULL)
  {}

  /*! Assign the contents of another topology-traits class. */
  void assign(const Self& other);
  //@}

  ///! \name Accessing the DCEL and constructing iterators.
  //@{

  /*! Determine whether the DCEL reprsenets an empty structure. */
  bool is_empty_dcel() const
  {
    // An empty bounded arrangement has no edges or vertices.
    return (this->m_dcel.size_of_vertices() == 0 &&
            this->m_dcel.size_of_halfedges() == 0);
  }

  /*! Check if the given vertex is concrete (associated with a point). */
  inline bool is_concrete_vertex(const Vertex*) const { return true; }

  /*! Get the number of concrete vertices. */
  Size number_of_concrete_vertices() const
  {
    // All vertices are concrete.
    return (this->m_dcel.size_of_vertices());
  }

  /*! Check if the given vertex is valid (not a fictitious one). */
  inline bool is_valid_vertex(const Vertex*) const { return true; }

  /*! Get the number of valid vertices. */
  Size number_of_valid_vertices() const
  {
    // All vertices are valid.
    return (this->m_dcel.size_of_vertices());
  }

  /*! Check if the given halfedge is valid (not a fictitious one). */
  inline bool is_valid_halfedge(const Halfedge*) const { return true; }

  /*! Get the number of valid halfedges. */
  Size number_of_valid_halfedges() const
  {
    // All halfedges are valid.
    return (this->m_dcel.size_of_halfedges());
  }

  /*! Check if the given face is valid (not a fictitious one). */
  inline bool is_valid_face (const Face*) const { return true; }

  /*! Get the number of valid faces. */
  Size number_of_valid_faces() const
  {
    // All faces are valid.
    return (this->m_dcel.size_of_faces());
  }
  //@}

private:
  /// \name Auxiliary type definitions.
  //@{
  typedef Arrangement_on_surface_2<Gt2, Self>                   Arr;

  // Type definition for the construction sweep-line visitor.
  typedef Arr_construction_event<Gt2, Arr>                      CEvent;
  typedef Arr_construction_subcurve<Gt2, CEvent>                CSubcurve;
  typedef Arr_bounded_planar_construction_helper<Gt2, Arr, CEvent, CSubcurve>
                                                                CHelper;

  // Type definition for the no-intersection construction sweep-line visitor.
  typedef Arr_construction_event<Gt2, Arr, Ss2::No_overlap_event_base,
                                 Ss2::No_overlap_subcurve>      NxCEvent;
  typedef Arr_construction_subcurve<Gt2, NxCEvent, Ss2::No_overlap_subcurve>
                                                                NxCSubcurve;
  typedef Arr_bounded_planar_construction_helper<Gt2, Arr, NxCEvent,
                                                 NxCSubcurve>   NxCHelper;

  // Type definition for the insertion sweep-line visitor.
  typedef Arr_insertion_traits_2<Gt2, Arr>                      ITraits;
  typedef Arr_construction_event<ITraits, Arr>                  IEvent;
  typedef Arr_construction_subcurve<ITraits, IEvent>            ISubcurve;
  typedef Arr_bounded_planar_insertion_helper<ITraits, Arr, IEvent, ISubcurve>
                                                                IHelper;

  // Type definition for the no-intersection insertion sweep-line visitor.
  // The second parameter of Arr_construction_subcurve is the base class of
  // Arr_construction_subcurve. By default Arr_construction_subcurve derives
  // from Surface_sweep_subcurve, which is suitable for general curves
  // including overlapping curves. Here we bypass Surface_sweep_subcurve and
  // force Arr_construction_subcurve to derive directly from
  // No_overlap_subcurve (which is the base class of Surface_sweep_subcurve).
  typedef Arr_basic_insertion_traits_2<Gt2, Arr>                NxITraits;
  typedef Arr_construction_event<NxITraits, Arr, Ss2::No_overlap_event_base,
                                 Ss2::No_overlap_subcurve>      NxIEvent;
  typedef Arr_construction_subcurve<NxITraits, NxIEvent,
                                    Ss2::No_overlap_subcurve>   NxISubcurve;
  typedef Arr_bounded_planar_insertion_helper<NxITraits, Arr, NxIEvent,
                                              NxISubcurve>      NxIHelper;

  // Type definition for the batched point-location sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              BplTraits;
  typedef Ss2::No_overlap_event<BplTraits>                      BplEvent;
  typedef Ss2::No_overlap_subcurve<BplTraits, BplEvent>         BplSubcurve;
  typedef Arr_bounded_planar_batched_pl_helper<BplTraits, Arr, BplEvent,
                                               BplSubcurve>     BplHelper;

  // Type definition for the vertical decomposition sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              VdTraits;
  typedef Ss2::No_overlap_event<VdTraits>                       VdEvent;
  typedef Ss2::No_overlap_subcurve<VdTraits, VdEvent>           VdSubcurve;
  typedef Arr_bounded_planar_vert_decomp_helper<VdTraits, Arr, VdEvent,
                                                VdSubcurve>     VdHelper;

  // Type definition for the overlay sweep-line visitor.
  template <typename ExtendedGeometryTraits, typename ArrangementA,
            typename ArrangementB>
  struct Overlay_helper : public Arr_bounded_planar_overlay_helper<
    ExtendedGeometryTraits, ArrangementA, ArrangementB, Arr,
    Arr_overlay_event<ExtendedGeometryTraits, Arr>,
    Arr_overlay_subcurve<ExtendedGeometryTraits, Arr_overlay_event<
                                                   ExtendedGeometryTraits,
                                                   Arr> > >
  {
    typedef ExtendedGeometryTraits                      Ex_geom_traits;
    typedef ArrangementA                                Arrangement_a;
    typedef ArrangementB                                Arrangement_b;

    typedef Arr_overlay_event<Ex_geom_traits, Arr>      Ex_event;
    typedef Arr_overlay_subcurve<Ex_geom_traits, Ex_event>
                                                        Ex_subcurve;
    typedef Arr_bounded_planar_overlay_helper<Ex_geom_traits, Arrangement_a,
                                              Arrangement_b, Arr, Ex_event,
                                              Ex_subcurve>
                                                        Base;

    typedef typename Base::Geometry_traits_2            Geometry_traits_2;
    typedef typename Base::Arrangement_red_2            Arrangement_red_2;
    typedef typename Base::Arrangement_blue_2           Arrangement_blue_2;
    typedef typename Base::Arrangement_2                Arrangement_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;
    typedef typename Base::Construction_helper          Construction_helper;

    Overlay_helper(const Arrangement_a* arr_a, const Arrangement_b* arr_b) :
      Base(arr_a, arr_b)
    {}
  };
  //@}

public:
  ///! \name Visitor types.
  //@{

  typedef Arr_construction_sl_visitor<CHelper>
    Surface_sweep_construction_visitor;

  typedef Arr_insertion_sl_visitor<IHelper>
    Surface_sweep_insertion_visitor;

  typedef Arr_construction_sl_visitor<NxCHelper>
    Surface_sweep_no_intersection_construction_visitor;

  typedef Arr_no_intersection_insertion_sl_visitor<NxIHelper>
    Surface_sweep_no_intersection_insertion_visitor;

  /*! \class Surface_sweep_batched_point_location_visitor
   */
  template <typename OutputIterator>
  struct Surface_sweep_batched_point_location_visitor :
    public Arr_batched_pl_sl_visitor<
      BplHelper, OutputIterator,
      Surface_sweep_batched_point_location_visitor<OutputIterator> >
  {
    typedef OutputIterator                              Output_iterator;

    typedef Surface_sweep_batched_point_location_visitor<Output_iterator>
                                                        Self;
    typedef Arr_batched_pl_sl_visitor<BplHelper, Output_iterator, Self>
                                                        Base;

    typedef typename BplHelper::Geometry_traits_2       Geometry_traits_2;
    typedef typename BplHelper::Event                   Event;
    typedef typename BplHelper::Subcurve                Subcurve;

    Surface_sweep_batched_point_location_visitor(const Arr* arr,
                                                 Output_iterator& oi) :
      Base(arr, oi)
    {}
  };

  /*! \class Surface_sweep_vertical_decomposition_visitor
   */
  template <typename OutputIterator>
  struct Surface_sweep_vertical_decomposition_visitor :
    public Arr_vert_decomp_sl_visitor<
      VdHelper, OutputIterator,
      Surface_sweep_vertical_decomposition_visitor<OutputIterator> >
  {
    typedef OutputIterator                              Output_iterator;
    typedef Surface_sweep_vertical_decomposition_visitor<Output_iterator>
                                                        Self;
    typedef Arr_vert_decomp_sl_visitor<VdHelper, Output_iterator, Self>
                                                        Base;

    typedef typename Base::Geometry_traits_2            Geometry_traits_2;
    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;

    Surface_sweep_vertical_decomposition_visitor(const Arr* arr,
                                                 Output_iterator* oi) :
      Base(arr, oi)
    {}
  };

  /*! \class Surface_sweep_overlay_visitor
   */
  template <typename ArrangementA, typename ArrangementB,
            typename OverlayTraits>
  struct Surface_sweep_overlay_visitor :
    public Arr_overlay_sl_visitor<Overlay_helper<Arr_overlay_traits_2<
                                                   Arr_traits_basic_adaptor_2<
                                                     Geometry_traits_2>,
                                                   ArrangementA, ArrangementB>,
                                                 ArrangementA, ArrangementB>,
                                  OverlayTraits,
                                  Surface_sweep_overlay_visitor<ArrangementA,
                                                                ArrangementB,
                                                                OverlayTraits> >
  {
    typedef ArrangementA                                Arrangement_a;
    typedef ArrangementB                                Arrangement_b;
    typedef OverlayTraits                               Overlay_traits;
    typedef Arr                                         Arrangement_result_2;

    typedef Arr_overlay_traits_2<Arr_traits_basic_adaptor_2<Geometry_traits_2>,
                                 Arrangement_a, Arrangement_b>
                                                        Geom_ovl_traits_2;

    typedef Overlay_helper<Geom_ovl_traits_2, Arrangement_a, Arrangement_b>
                                                        Ovl_helper;

    typedef Surface_sweep_overlay_visitor<Arrangement_a, Arrangement_b,
                                          Overlay_traits>
                                                        Self;
    typedef Arr_overlay_sl_visitor<Ovl_helper, Overlay_traits, Self>
                                                        Base;

    typedef typename Base::Event                        Event;
    typedef typename Base::Subcurve                     Subcurve;

    Surface_sweep_overlay_visitor(const Arrangement_a* arr_a,
                                  const Arrangement_b* arr_b,
                                  Arrangement_result_2* arr_res,
                                  Overlay_traits* overlay_tr) :
      Base(arr_a, arr_b, arr_res, overlay_tr)
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

  /*! Initialize an empty DCEL structure.
   */
  void init_dcel();

  /*! Make the necessary updates after the DCEL structure have been updated.
   */
  void dcel_updated();

  /*! Check if the given vertex is associated with the given curve end.
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
                 Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
  {
    CGAL_USE(ps_x);
    CGAL_USE(ps_y);

    CGAL_assertion((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR));

    if (ind == ARR_MIN_END) {
      // Compare v with the left endpoint of cv.
      return (this->m_geom_traits->equal_2_object()
              (this->m_geom_traits->construct_min_vertex_2_object()(cv),
               v->point()));
    }
    else {
      // Compare v with the right endpoint of cv.
      return (this->m_geom_traits->equal_2_object()
              (this->m_geom_traits->construct_max_vertex_2_object()(cv),
               v->point()));
    }
  }

  /*! Given a curve end with boundary conditions and a face that contains the
   * interior of the curve, find a place for a boundary vertex that will
   * represent the curve end along the face boundary.
   * \param f The face.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return An object that contains the curve end.
   */
  CGAL::Object place_boundary_vertex(Face*,
                                     const X_monotone_curve_2&,
                                     Arr_curve_end,
                                     Arr_parameter_space /* ps_x */,
                                     Arr_parameter_space /* ps_y */)
  {
    // This function should never be called:
    CGAL_error();
    return Object();
  }

  /*! Locate the predecessor halfedge for the given curve around a given
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
  locate_around_boundary_vertex(Vertex*,
                                const X_monotone_curve_2&,
                                Arr_curve_end,
                                Arr_parameter_space /* ps_x */,
                                Arr_parameter_space /* ps_y */) const
  {
    CGAL_error();
    return NULL;
  }

  /*! Locate a DCEL feature that contains the given curve end.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve end is incident to the boundary.
   * \return An object that contains the curve end.
   */
  CGAL::Object locate_curve_end(const X_monotone_curve_2&,
                                Arr_curve_end,
                                Arr_parameter_space /* ps_x */,
                                Arr_parameter_space /* ps_y */)
  {
    // This function should never be called:
    CGAL_error();
    return Object();
  }

  /*! Split a fictitious edge using the given vertex.
   * \param e The edge to split (one of the pair of halfedges).
   * \param v The split vertex.
   * \pre e is a fictitious halfedge.
   * \return A halfedge whose direction is the same as e's and whose target is
   *         the split vertex v.
   */
  Halfedge* split_fictitious_edge (Halfedge*, Vertex*)
  {
    // This function should never be called:
    CGAL_error();
    return NULL;
  }

  /*! Determine whether the given face is unbounded.
   * \param f The face.
   * \return Whether f is unbounded.
   * There is only one unbounded face in the arrangement:
   */
  bool is_unbounded(const Face* f) const { return (f == unb_face); }

  /*! Determine whether the given boundary vertex is redundant.
   * \param v The vertex.
   * \return Whether v is redundant, and should be erased.
   * There are no redundant vertices.
   */
  bool is_redundant(const Vertex*) const { return false; }

  /*! Erase the given redundant vertex by merging a fictitious edge.
   * The function does not free the vertex v itself.
   * \param v The vertex.
   * \pre v is a redundant vertex.
   * \return One of the pair of halfedges that form the merged edge.
   */
  Halfedge* erase_redundant_vertex(Vertex*)
  {
    // This function should never be called:
    CGAL_error();
    return NULL;
  }

    //! reference_face (const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common
      point.
      \return A pointer to the reference face.
  */
  const Face* reference_face() const { return unbounded_face(); }

  //! reference_face (non-const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common
      point.
      \return A pointer to the reference face.
  */
  Face* reference_face() { return unbounded_face(); }

  //@}

  /// \name Additional accessors, specialized for this topology-traits class.
  //@{

  /*! This function is used by the "walk" point-location strategy. */
  const Face* initial_face() const { return (unb_face); }

  /*! Get the unbounded face (const version). */
  const Face* unbounded_face() const { return (unb_face); }

  /*! Get the unbounded face (non-const version). */
  Face* unbounded_face() { return (unb_face); }
  //@}

  /// \name Additional predicates, specialized for this topology-traits class.
  //@{

  /*! Compare the given vertex and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the comparison of the x-coordinates of p and v.
   */
  virtual Comparison_result compare_x(const Point_2& p, const Vertex* v) const
  { return (this->m_geom_traits->compare_x_2_object()(p, v->point())); }

  /*! Compare the given vertex and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the xy-lexicographic comparison of p and v.
   */
  virtual Comparison_result compare_xy(const Point_2& p, const Vertex* v) const
  { return (this->m_geom_traits->compare_xy_2_object()(p, v->point())); }

  /*! Compare the relative y-position of the given point and the given edge
   * (which may be fictitious).
   * \param p The point.
   * \param he The edge (one of the pair of halfedges).
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  virtual Comparison_result compare_y_at_x(const Point_2& p,
                                           const Halfedge* he) const
  { return (this->m_geom_traits->compare_y_at_x_2_object()(p, he->curve())); }
  //@}
};

} // namespace CGAL

#include <CGAL/Arr_topology_traits/Arr_bounded_planar_topology_traits_2_impl.h>

#endif
