// Copyright (c) 2006,2007,2008,2009,2010,2011,2012,2013 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Eric Berberich    <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_H
#define CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 *
 * The topology traits for great spherical arcs embedded on a sphere for the
 * arrangement package.
 */

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Surface_sweep_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Surface_sweep_2/Arr_insertion_traits_2.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_traits_2.h>
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
template <typename GeomTraits, typename TopTraits>
class Arrangement_on_surface_2;

/*! This class handles the topology for arrangements of great spherical
 * arcs on the sphere embedded on 2D parametric surdace.
 */
template <typename GeometryTraits_2,
          typename Dcel_ = Arr_default_dcel<GeometryTraits_2> >
class Arr_spherical_topology_traits_2 {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Dcel_                                         Dcel;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  ///! \name The geometry-traits types.
  //@{
  typedef typename Gt2::Point_2                         Point_2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  //@}

  ///! \name The DCEL types.
  //@{
  typedef typename Dcel::Size                           Size;
  typedef typename Dcel::Vertex                         Vertex;
  typedef typename Dcel::Halfedge                       Halfedge;
  typedef typename Dcel::Face                           Face;
  typedef typename Dcel::Outer_ccb                      Outer_ccb;
  typedef typename Dcel::Inner_ccb                      Inner_ccb;
  typedef typename Dcel::Isolated_vertex                Isolated_vertex;
  //@}

  //! \name Arrangement types
  //!@{
  typedef Arr_spherical_topology_traits_2<Gt2, Dcel>    Self;
  typedef Arr_traits_basic_adaptor_2<Gt2>               Gt_adaptor_2;
  //!@}

  ///! \name The side tags
  //@{
  typedef typename Gt_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Gt_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Gt_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Gt_adaptor_2::Right_side_category  Right_side_category;

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
  template <typename T, typename D>
  struct rebind {
    typedef Arr_spherical_topology_traits_2<T, D> other;
  };

private:
  //! A container of boundary vertices.
  struct Vertex_key_comparer {
    /*! Construct default */
    Vertex_key_comparer() : m_geom_traits(nullptr) {}

    /*! Construct */
    Vertex_key_comparer(const Gt_adaptor_2* geom_traits) :
      m_geom_traits(geom_traits)
    {}

    const Gt_adaptor_2* m_geom_traits;

    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      return (m_geom_traits->compare_y_on_boundary_2_object()(p1, p2) ==
              SMALLER);
    }
  };

  //! \todo define the key to be 'const Point_2*'.
  typedef std::map<Point_2, Vertex*, Vertex_key_comparer> Vertex_map;
  typedef std::pair<Point_2, Vertex*>                     Vertex_value;

protected:
  // Data members:
  //! The DCEL.
  Dcel m_dcel;

  //! The spherical (bounded) face.
  Face* m_spherical_face;

  //! The north pole
  Vertex* m_north_pole;

  //! The north pole
  Vertex* m_south_pole;

  //! The vertices on the discontinuity arc
  Vertex_map m_boundary_vertices;

  //! The geometry-traits adaptor.
  const Gt_adaptor_2* m_geom_traits;

  //! Inidicates whether the traits object should evetually be freed.
  bool m_own_geom_traits;

  // Copy constructor and assignment operator - not supported.
  Arr_spherical_topology_traits_2(const Self&);
  Self& operator=(const Self&);

public:
  ///! \name Construction methods.
  //@{

  /*! Default constructor. */
  Arr_spherical_topology_traits_2();

  /*! Constructor from a geometry-traits object.
   * \param traits the traits.
   */
  Arr_spherical_topology_traits_2(const Gt2* traits);

  /*! Destructor */
  ~Arr_spherical_topology_traits_2();

  /*! Assign the contents of another topology-traits class.
   * \param other the other spherical topology-traits.
   */
  void assign(const Self& other);
  //@}

  ///! \name Topology-traits methods.
  //@{

  /*! Obtain the DCEL (const version). */
  const Dcel& dcel() const { return (m_dcel); }

  /*! Obtain the DCEL (non-const version). */
  Dcel& dcel() { return m_dcel; }

  /*! Determine whether the DCEL reprsenets an empty structure.
   * \return true if the dcel reprsenets an empty structure; false otherwise.
   */
  bool is_empty_dcel() const { return (m_dcel.size_of_vertices() == 0); }

  /*! Initialize an empty DCEL structure. */
  void init_dcel();

  /*! Make the necessary updates after the DCEL structure have been updated. */
  void dcel_updated();

  /*! Determine whether the given vertex is concrete.
   * \param v the vertex.
   * \return true if v is mapped to a point on the discontinuity arc; false
   * otherwise.
   */
  bool is_concrete_vertex(const Vertex* /* v */) const { return true; }

  /*! Obtain the number of concrete vertices.
   * \return the number of concrete vertices.
   */
  Size number_of_concrete_vertices() const
  { return (m_dcel.size_of_vertices()); }

  /*! Determine whether the given vertex is valid.
   * \param v the vertex.
   * \todo why is this needed, and where used?
   */
  bool is_valid_vertex (const Vertex* /* v */) const { return true; }

  /*! Obtain the number of valid vertices. */
  Size number_of_valid_vertices() const { return (m_dcel.size_of_vertices()); }

  /*! Determine whether the given halfedge is valid. */
  bool is_valid_halfedge (const Halfedge* /* he */) const { return true; }

  /*! Obtain the number of valid halfedges. */
  Size number_of_valid_halfedges() const
  { return (m_dcel.size_of_halfedges()); }

  /*! Determine whether the given face is valid. */
  bool is_valid_face(const Face* /* f */) const { return true; }

  /*! Obtain the number of valid faces. */
  Size number_of_valid_faces() const { return m_dcel.size_of_faces(); }

  /*! Obtain the spherical face (const version). */
  const Face* spherical_face() const { return m_spherical_face; }

  /*! Obtain the spherical face (non-const version). */
  Face* spherical_face() { return m_spherical_face; }

  /*! Obtain the face containing the south pole (const version). */
  const Face* south_face() const
  {
    if (m_boundary_vertices.empty()) return m_spherical_face;
    typename Vertex_map::const_iterator it = m_boundary_vertices.begin();
    return _face_below_vertex_on_discontinuity(it->second);
  }

  /*! Obtain the face containing the south pole (non-const version). */
  Face* south_face()
  {
    if (m_boundary_vertices.empty()) return m_spherical_face;
    typename Vertex_map::iterator it = m_boundary_vertices.begin();
    return _face_below_vertex_on_discontinuity(it->second);
  }

  /*! Obtain the south pole (const version). */
  const Vertex* south_pole() const { return m_south_pole; }

  /*! Obtain the south pole (non-const version). */
  Vertex* south_pole() { return m_south_pole; }

  /*! Obtain the north pole (const version). */
  const Vertex* north_pole() const { return m_north_pole; }

  /*! Obtain the north pole (non-const version). */
  Vertex* north_pole() { return m_north_pole; }

  /*! Obtain a vertex on the line of discontinuity that corresponds to
   *  the given curve-end (or return nullptr if no such vertex exists).
   */
  Vertex* discontinuity_vertex(const X_monotone_curve_2 xc, Arr_curve_end ind)
  {
    Point_2 key = (ind == ARR_MIN_END) ?
      m_geom_traits->construct_min_vertex_2_object()(xc) :
      m_geom_traits->construct_max_vertex_2_object()(xc);
    typename Vertex_map::iterator it = m_boundary_vertices.find(key);
    return (it != m_boundary_vertices.end()) ? it->second : nullptr;
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
  // derived class. The non-C++11 code follows the commented out C++11 code.
  // When we move to C++11 we can use the more elgant code.
  // Type definition for the construction surface-sweep visitor.
  template <typename Evt, typename Crv>
  struct Construction_helper :
    public Arr_spherical_construction_helper<Gt2, Arr, Evt, Crv>
  {
    typedef Arr_spherical_construction_helper<Gt2, Arr, Evt, Crv>
                                                                Base;
    Construction_helper(Arr* arr) : Base(arr) {}
  };

  // Type definition for the no-intersection construction surface-sweep visitor.
  template <typename Evt, typename Crv>
  struct No_intersection_construction_helper :
    public Arr_spherical_construction_helper<Gt2, Arr, Evt, Crv>
  {
    typedef Arr_spherical_construction_helper<Gt2, Arr, Evt, Crv>
                                                                Base;
    No_intersection_construction_helper(Arr* arr) : Base(arr) {}
  };

  // Type definition for the insertion surface-sweep visitor.
  typedef Arr_insertion_traits_2<Gt2, Arr>                      I_traits;
  template <typename Evt, typename Crv>
  struct Insertion_helper :
    public Arr_spherical_insertion_helper<I_traits, Arr, Evt, Crv>
  {
    typedef Arr_spherical_insertion_helper<I_traits, Arr, Evt, Crv>
                                                                Base;
    Insertion_helper(Arr* arr) : Base(arr) {}
  };

  // Type definition for the no-intersection insertion surface-sweep visitor.
  typedef Arr_basic_insertion_traits_2<Gt2, Arr>                Nxi_traits;
  template <typename Evt, typename Crv>
  struct No_intersection_insertion_helper :
    public Arr_spherical_insertion_helper<Nxi_traits, Arr, Evt, Crv>
  {
    typedef Arr_spherical_insertion_helper<Nxi_traits, Arr, Evt, Crv>
                                                                Base;
    No_intersection_insertion_helper(Arr* arr) : Base(arr) {}
  };

  // Type definition for the batched point-location surface-sweep visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              Bpl_traits;
  template <typename Evt, typename Crv>
  struct Batched_point_location_helper :
    public Arr_spherical_batched_pl_helper<Bpl_traits, Arr, Evt, Crv>
  {
    typedef Arr_spherical_batched_pl_helper<Bpl_traits, Arr, Evt, Crv>
                                                                Base;
    Batched_point_location_helper(const Arr* arr) : Base(arr) {}
  };

  // Type definition for the vertical decomposition surface-sweep visitor.
  typedef Arr_batched_point_location_traits_2<Arr>              Vd_traits;
  template <typename Evt, typename Crv>
  struct Vertical_decomposition_helper :
    public Arr_spherical_vert_decomp_helper<Vd_traits, Arr, Evt, Crv>
  {
    typedef Arr_spherical_vert_decomp_helper<Vd_traits, Arr, Evt, Crv>
                                                                Base;
    Vertical_decomposition_helper(const Arr* arr) : Base(arr) {}
  };

  // Type definition for the overlay surface-sweep visitor.
  template <typename Gt, typename Evt, typename Crv,
            typename ArrA, typename ArrB>
  struct Overlay_helper :
    public Arr_spherical_overlay_helper<Gt, ArrA, ArrB, Arr, Evt, Crv>
  {
    typedef Arr_spherical_overlay_helper<Gt, ArrA, ArrB, Arr, Evt, Crv>
                                                                Base;
    Overlay_helper(const ArrA* arr_a, const ArrB* arr_b) : Base(arr_a, arr_b) {}
  };
  //@}

public:
  ///! \name Visitor types.
  //@{

  typedef Arr_inc_insertion_zone_visitor<Arr> Zone_insertion_visitor;

  typedef Arr_naive_point_location<Arr> Default_point_location_strategy;
  typedef Arr_naive_point_location<Arr> Default_vertical_ray_shooting_strategy;
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
  void notify_on_boundary_vertex_creation(Vertex* v,
                                          const X_monotone_curve_2& xc,
                                          Arr_curve_end ind,
                                          Arr_parameter_space ps_x,
                                          Arr_parameter_space ps_y);

  /*! Determines whether the function should decide on swapping the predecssor
   * halfedges that imply two ccb (and whose signs are given here).
   * If true, swap_predecessors will be correctly set. If false,
   * generic way of searching for lexicographically minimal point and checking
   * its incident halfedges will do the job to decide on swapping
   * \param signs1 signs of first implied ccb in x- and y-direction
   * \param signs2 signs of second implied ccb in x- and y-direction
   * \param swap_predecessors Output swap predeccesors or not;
   *        set correctly only if true is returned
   */
  bool let_me_decide_the_outer_ccb(std::pair< CGAL::Sign, CGAL::Sign> signs1,
                                   std::pair< CGAL::Sign, CGAL::Sign> signs2,
                                   bool& swap_predecessors) const;


  /*! Given signs of two ccbs that show up when splitting upon insertion of
   * curve into two, determine what happens to the face(s).
   * \param signs1 signs in x and y of the first implied ccb
   * \param signs2 signs in x and y of the secondd implied ccb
   * \return A pair indicating whether the insertion will cause the face
   *         to split (the first flag), and if so - whether the split face
   *         will form a hole in the original face.
   */
  std::pair<bool, bool>
  face_split_after_edge_insertion(std::pair< CGAL::Sign,
                                             CGAL::Sign > /* signs1 */,
                                  std::pair< CGAL::Sign,
                                             CGAL::Sign > /* signs2 */) const
  {
    // In case of a spherical topology, connecting two vertices on the same
    // inner CCB closes a new face that becomes a hole in the original face:
    return (std::make_pair(true, true));
  }

  /*! Determine whether a given point lies in the interior of a given face.
   * \param f The face.
   * \param p The query point.
   * \param v The vertex associated with p (if exists).
   * \param f must not be fictitious, and v must not lie at infinity.
   * \return true if p is contained in f's interior; false otherwise.
   */
  bool is_in_face(const Face* f, const Point_2& p, const Vertex* v) const;

  /*! Compare the relative y-position of a given point and a given edge.
   * \param p The point.
   * \param he The edge (one of the pair of halfedges).
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  Comparison_result compare_y_at_x(const Point_2& p,
                                   const Halfedge* he) const;

  /*! Determine whether a given vertex is associated with a given curve end.
   * \param v The vertex.
   * \param xc The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \pre The curve has a boundary condition in either x or y.
   * \return Whether v represents the given curve end.
   */
  bool are_equal(const Vertex* v,
                 const X_monotone_curve_2& xc, Arr_curve_end ind,
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
  boost::optional<boost::variant<Vertex*, Halfedge*> >
  place_boundary_vertex(Face* f,
                        const X_monotone_curve_2& xc,
                        Arr_curve_end ind,
                        Arr_parameter_space ps_x,
                        Arr_parameter_space ps_y);

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
  Halfedge* locate_around_boundary_vertex(Vertex* v,
                                          const X_monotone_curve_2& cv,
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
  boost::variant<Vertex*, Halfedge*, Face*>
  locate_curve_end(const X_monotone_curve_2& xc, Arr_curve_end ce,
                   Arr_parameter_space ps_x,
                   Arr_parameter_space ps_y);

  /*! Split a fictitious edge using the given vertex.
   * \param e The edge to split (one of the pair of halfedges).
   * \param v The split vertex.
   * \pre e is a fictitious halfedge.
   * \return A halfedge whose direction is the same as e's and whose target is
   *         the split vertex v.
   */
  Halfedge* split_fictitious_edge(Halfedge* /* e */, Vertex* /* v */)
  {
    // There are no fictitious halfedges:
    CGAL_error();
    return nullptr;
  }

  /*! Determine whether the given face is unbounded.
   * \param f The face.
   * \return true if f is unbounded; false otherwise.
   * All faces on a sphere are bounded:
   */
  bool is_unbounded(const Face* /* f */) const { return false; }

  /*! Determine whether the given boundary vertex is redundant.
   * \param v The vertex.
   * \return Whether v is redundant, and should be erased.
   */
  bool is_redundant(const Vertex* v) const;

  /*! Erase the given redundant vertex by merging a fictitious edge.
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
  const Face* reference_face() const { return spherical_face(); }

  //! reference_face (non-const version).
  /*! The function returns a reference face of the arrangement.
      All reference faces of arrangements of the same type have a common
      point.
      \return A pointer to the reference face.
  */
  Face* reference_face() { return spherical_face(); }
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
  const X_monotone_curve_2& _curve(const Vertex* v, Arr_curve_end& ind) const;

  /*! Return the halfedge, the target vertex of which is given, that is
   * the predecessor of a halfedge, the curve of which is given, that is about
   * to be inserted into the dcel.
   */
  Halfedge* _locate_around_vertex_on_discontinuity(Vertex* v,
                                                   const X_monotone_curve_2& xc,
                                                   Arr_curve_end ind) const;

  /*! Return the halfedge, the target vertex of which is a given pole,
   * that is the predecessor of a halfedge, the curve of which is given, that
   * is about to be inserted into the dcel.
   */
  Halfedge* _locate_around_pole(Vertex* v, const X_monotone_curve_2& xc,
                                Arr_curve_end ind) const;


  /*! Return the face that lies below the given vertex, which lies
   * on the line of discontinuity.
   */
  Face* _face_below_vertex_on_discontinuity(Vertex* v) const;

  //@}
};

} // namespace CGAL

#include <CGAL/Arr_topology_traits/Arr_spherical_topology_traits_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif
