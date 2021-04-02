// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein        <wein@post.tau.ac.il>
//            Efi Fogel       <efif@post.tau.ac.il>
//            Eric Berberich  <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_BOUNDED_PLANAR_TOPOLOGY_TRAITS_2_H

#include <boost/variant.hpp>

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

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
#include <CGAL/use.h>

namespace CGAL {

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
  typedef Arr_traits_basic_adaptor_2<Gt2>                 Gt_adaptor_2;
  //!@}

  ///! \name The side tags
  //@{
  typedef typename Gt_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Gt_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Gt_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Gt_adaptor_2::Right_side_category  Right_side_category;

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
    unb_face(nullptr)
  {}

  /*! Constructor from a geometry-traits object. */
  Arr_bounded_planar_topology_traits_2(const Gt2* traits) :
    Base(traits),
    unb_face(nullptr)
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

public:
  // The following definitions of helper functions use a C++11 features called
  // alias template. This feature can be mimiked by non C++11 support simply by
  // defining a type that inherits from the template we need to alias.  However,
  // the non-C++11 code requires the (re)definition of all constructors of the
  // derived class.
  // Type definition for the construction surface-sweep visitor.
  template <typename Evt, typename Crv>
  using Construction_helper =
    Arr_bounded_planar_construction_helper<Gt2, Arr, Evt, Crv>;

  // Type definition for the no-intersection construction surface-sweep visitor.
  template <typename Evt, typename Crv>
  using No_intersection_construction_helper =
    Arr_bounded_planar_construction_helper<Gt2, Arr, Evt, Crv>;

  // Type definition for the insertion surface-sweep visitor.
  typedef Arr_insertion_traits_2<Gt2, Arr>                      I_traits;
  template <typename Evt, typename Crv>
  using Insertion_helper =
    Arr_bounded_planar_insertion_helper<I_traits, Arr, Evt, Crv>;

  // Type definition for the no-intersection insertion surface-sweep visitor.
  typedef Arr_basic_insertion_traits_2<Gt2, Arr>                Nxi_traits;
  template <typename Evt, typename Crv>
  using No_intersection_insertion_helper =
    Arr_bounded_planar_insertion_helper<Nxi_traits, Arr, Evt, Crv>;

  // Type definition for the batched point-location surface-sweep visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              Bpl_traits;
  template <typename Evt, typename Crv>
  using Batched_point_location_helper =
    Arr_bounded_planar_batched_pl_helper<Bpl_traits, Arr, Evt, Crv>;

  // Type definition for the vertical decomposition sweep-line visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              Vd_traits;
  template <typename Evt, typename Crv>
  using Vertical_decomposition_helper =
    Arr_bounded_planar_vert_decomp_helper<Vd_traits, Arr, Evt, Crv>;

  // Type definition for the overlay surface-sweep visitor.
  template <typename Gt, typename Evt, typename Crv,
            typename ArrA, typename ArrB>
  using Overlay_helper =
    Arr_bounded_planar_overlay_helper<Gt, ArrA, ArrB, Arr, Evt, Crv>;
  //@}

public:
  ///! \name Visitor types.
  //@{

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
   * \return An object that wraps the curve end.
   */
  boost::optional<boost::variant<Vertex*, Halfedge*> >
  place_boundary_vertex(Face*,
                        const X_monotone_curve_2&,
                        Arr_curve_end,
                        Arr_parameter_space /* ps_x */,
                        Arr_parameter_space /* ps_y */)
  {
    // This function should never be called:
    CGAL_error();
    return boost::none;
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
    return nullptr;
  }

  /*! Locate a DCEL feature that contains the given curve end.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve end is incident to the boundary.
   * \return An object that contains the curve end.
   */
  boost::variant<Vertex*, Halfedge*, Face*>
  locate_curve_end(const X_monotone_curve_2&,
                   Arr_curve_end,
                   Arr_parameter_space /* ps_x */,
                   Arr_parameter_space /* ps_y */)
  {
    typedef boost::variant<Vertex*, Halfedge*, Face*>   Result;
    // This function should never be called:
    CGAL_error();
    Vertex* v(nullptr);
    return Result(v);
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
    return nullptr;
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
    return nullptr;
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

#include <CGAL/enable_warnings.h>

#endif
