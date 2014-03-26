#ifndef CGAL_OVERLAY_TEST_H
#define CGAL_OVERLAY_TEST_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/IO/Arr_text_formatter.h>

#include "utils.h"
#include "IO_base_test.h"

/*! Overlay test */
template <typename T_Geom_traits, typename T_Topol_traits>
class Overlay_test :
  public IO_base_test<typename T_Geom_traits::Base_traits_2>
{
public:
  typedef T_Geom_traits                                         Geom_traits;
  typedef T_Topol_traits                                        Topol_traits;
  typedef IO_base_test<typename Geom_traits::Base_traits_2>     Base;

  typedef typename Geom_traits::Base_traits_2           Base_geom_traits;
  typedef typename Base_geom_traits::Point_2            Base_point_2;
  typedef typename Base_geom_traits::Curve_2            Base_curve_2;
  typedef typename Base_geom_traits::X_monotone_curve_2 Base_xcurve_2;

  typedef typename Geom_traits::Point_2                 Point_2;
  typedef typename Geom_traits::Curve_2                 Curve_2;
  typedef typename Geom_traits::X_monotone_curve_2      Xcurve_2;

  typedef typename std::vector<Curve_2>                 Curve_container;
  typedef typename std::vector<Point_2>                 Point_container;
  typedef typename std::vector<Xcurve_2>                Xcurve_container;

  typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                        Arrangement;

  typedef CGAL::Arr_extended_dcel_text_formatter<Arrangement>
                                                        Formatter;

  // Handles
  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;

  // Const handles
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;

  // Iterators
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement::Face_iterator           Face_iterator;
  typedef typename Arrangement::Outer_ccb_iterator      Outer_ccb_iterator;
  typedef typename Arrangement::Inner_ccb_iterator      Inner_ccb_iterator;

  // Const iterators
  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Arrangement::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement::Outer_ccb_const_iterator
    Outer_ccb_const_iterator;
  typedef typename Arrangement::Inner_ccb_const_iterator
    Inner_ccb_const_iterator;

  // Circulators
  typedef typename Arrangement::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Arrangement::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;

  // Const Circulators
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

private:
  /*! The geometry traits */
  const Geom_traits& m_geom_traits;

  /*! Verbosity */
  unsigned int m_verbose_level;

  typename Arrangement::Size m_num_vertices;
  typename Arrangement::Size m_num_edges;
  typename Arrangement::Size m_num_faces;

  /*! The input data file of points*/
  std::string m_filename;

  // Arrangement 1
  Arrangement m_arr1;
  Xcurve_container m_xcurves1;
  Point_container m_isolated_points1;

  // Arrangement 2
  Arrangement m_arr2;
  Xcurve_container m_xcurves2;
  Point_container m_isolated_points2;

  // Arrangement expected
  Arrangement m_arr;

public:
  /*! Constructor */
  Overlay_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Overlay_test() { clear(); }

  void set_verbose_level(unsigned int verbose_level)
  { m_verbose_level = verbose_level; }

  void set_filename(const char* filename) { m_filename.assign(filename); }

  /*! Initialize the test */
  virtual bool init();

  /*! Perform the test */
  virtual bool perform();

  /*! Clear the data structures */
  virtual void clear();

protected:
  //! Overlay traits
  class Overlay_traits {
  private:
    /*! Verbosity */
    unsigned int m_verbose_level;

    std::size_t count_outer(Face_const_handle f) const
    {
      std::size_t cnt = 0;
      Outer_ccb_const_iterator ocit;
      for (ocit = f->outer_ccbs_begin(); ocit != f->outer_ccbs_end(); ++ocit) {
        Ccb_halfedge_const_circulator curr = *ocit;
        do ++cnt;
        while (++curr != *ocit);
      }
      return cnt;
    }

    std::size_t count_outer(Face_handle f) const
    {
      std::size_t cnt = 0;
      Outer_ccb_iterator ocit;
      for (ocit = f->outer_ccbs_begin(); ocit != f->outer_ccbs_end(); ++ocit) {
        Ccb_halfedge_circulator curr = *ocit;
        do ++cnt;
        while (++curr != *ocit);
      }
      return cnt;
    }

    std::size_t count_inner(Face_const_handle f) const
    {
      std::size_t cnt = 0;
      Inner_ccb_const_iterator icit;
      for (icit = f->inner_ccbs_begin(); icit != f->inner_ccbs_end(); ++icit) {
        Ccb_halfedge_const_circulator curr = *icit;
        do ++cnt;
        while (++curr != *icit);
      }
      return cnt;
    }

    std::size_t count_inner(Face_handle f) const
    {
      std::size_t cnt = 0;
      Inner_ccb_iterator icit;
      for (icit = f->inner_ccbs_begin(); icit != f->inner_ccbs_end(); ++icit) {
        Ccb_halfedge_circulator curr = *icit;
        do ++cnt;
        while (++curr != *icit);
      }
      return cnt;
    }

  public:
    /*! Destructor. */
    virtual ~Overlay_traits() {}

    /*! Constructor */
    Overlay_traits(unsigned int verbose_level)
    { m_verbose_level = verbose_level; }

    /*! Create a vertex v that corresponds to the coinciding vertices v1 and v2.
     */
    virtual void create_vertex(Vertex_const_handle v1, Vertex_const_handle v2,
                               Vertex_handle v) const
    {
#ifdef CGAL_OVERLAY_TRAITS_VERBOSE
      std::cout << "  v1: " << v1->point() << ", " << v1->data() << std::endl;
      std::cout << "  v2: " << v2->point() << ", " << v2->data() << std::endl;
#endif
      v->set_data(v1->data() + v2->data());
#ifdef CGAL_OVERLAY_TRAITS_VERBOSE
      std::cout << "  v: " << v->point() << ", " << v->data() << std::endl;
      std::cout << std::endl;
#endif
    }

    /*! Create a vertex v that mathces v1, which lies of the edge e2. */
    virtual void create_vertex(Vertex_const_handle  v1,
                               Halfedge_const_handle e2,
                               Vertex_handle v) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  v1: " << v1->point() << ", " << v1->data() << std::endl;
        std::cout << "  e2: " << e2->source()->point() << "=>"
                  << e2->target()->point()
                  << ", " << e2->data() << std::endl;
      }
      v->set_data(v1->data() + e2->data());
      if (m_verbose_level > 2) {
        std::cout << "  v: " << v->point() << ", " << v->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create a vertex v that mathces v1, contained in the face f2. */
    virtual void create_vertex(Vertex_const_handle v1, Face_const_handle f2,
                               Vertex_handle v) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  v1: " << v1->point() << ", " << v1->data() << std::endl;
        std::cout << "  f2: " << "Inner(" << count_inner(f2) << ")"
                  << ", Outer(" << count_outer(f2) << ")"
                  << ", " << f2->data() << std::endl;
      }
      v->set_data(v1->data() + f2->data());
      if (m_verbose_level > 2) {
        std::cout << "  V: " << v->point() << ", " << v->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create a vertex v that mathces v2, which lies of the edge e1. */
    virtual void create_vertex(Halfedge_const_handle e1, Vertex_const_handle v2,
                               Vertex_handle v) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  e1: " << e1->source()->point() << "=>"
                  << e1->target()->point()
                  << ", " << e1->data() << std::endl;
        std::cout << "  v2: " << v2->point() << ", " << v2->data() << std::endl;
      }
      v->set_data(e1->data() + v2->data());
      if (m_verbose_level > 2) {
        std::cout << "  v: " << v->point() << ", " << v->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create a vertex v that mathces v2, contained in the face f1. */
    virtual void create_vertex(Face_const_handle f1, Vertex_const_handle v2,
                               Vertex_handle v) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  f1: " << "Inner(" << count_inner(f1) << ")"
                  << ", Outer(" << count_outer(f1) << ")"
                  << ", " << f1->data() << std::endl;
        std::cout << "  v2: " << v2->point() << ", " << v2->data() << std::endl;
      }
      v->set_data(f1->data() + v2->data());
      if (m_verbose_level > 2) {
        std::cout << "  v: " << v->point() << ", " << v->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create a vertex v that mathces the intersection of the edges e1 and e2.
     */
    virtual void create_vertex(Halfedge_const_handle e1,
                               Halfedge_const_handle e2,
                               Vertex_handle v) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  e1: " << e1->source()->point() << "=>"
                  << e1->target()->point()
                  << ", " << e1->data() << std::endl;
        std::cout << "  e2: " << e2->source()->point()  << "=>"
                  << e2->target()->point()
                  << ", " << e2->data() << std::endl;
      }
      v->set_data(e1->data() + e2->data());
      if (m_verbose_level > 2) {
        std::cout << "  v: " << v->point() << ", " << v->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create an edge e that matches the overlap between e1 and e2. */
    virtual void create_edge(Halfedge_const_handle e1, Halfedge_const_handle e2,
                             Halfedge_handle e) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  e1: " << e1->source()->point() << "=>"
                  << e1->target()->point()
                  << ", " << e1->data() << std::endl;
        std::cout << "  e2: " << e2->source()->point() << "=>"
                  << e2->target()->point()
                  << ", " << e2->data() << std::endl;
      }
      e->set_data(e1->data() + e2->data());
      e->twin()->set_data(e1->data() + e2->data());
      if (m_verbose_level > 2) {
        std::cout << "  e: " << e->source()->point() << "=>"
                  << e->target()->point()
                  << ", " << e->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create an edge e that matches the edge e1, contained in the face f2. */
    virtual void create_edge(Halfedge_const_handle e1, Face_const_handle f2,
                             Halfedge_handle e) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  e1: " << e1->source()->point() << "=>"
                  << e1->target()->point()
                  << ", " << e1->data()
                  << std::endl;
        std::cout << "  f2: " << "Inner(" << count_inner(f2) << ")"
                  << ", Outer(" << count_outer(f2) << ")"
                  << ", " << f2->data() << std::endl;
      }
      e->set_data(e1->data() + f2->data());
      e->twin()->set_data(e1->data() + f2->data());
      if (m_verbose_level > 2) {
        std::cout << "  e: " << e->source()->point() << "=>"
                  << e->target()->point()
                  << ", " << e->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create an edge e that matches the edge e2, contained in the face f1. */
    virtual void create_edge(Face_const_handle f1, Halfedge_const_handle e2,
                             Halfedge_handle e) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  f1: " << "Inner(" << count_inner(f1) << ")"
                  << ", Outer(" << count_outer(f1) << ")"
                  << ", " << f1->data() << std::endl;
        std::cout << "  e2: " << e2->source()->point() << "=>"
                  << e2->target()->point()
                  << ", " << e2->data() << std::endl;
      }
      e->set_data(f1->data() + e2->data());
      e->twin()->set_data(f1->data() + e2->data());
      if (m_verbose_level > 2) {
        std::cout << "  e: " << e->source()->point() << "=>"
                  << e->target()->point()
                  << ", " << e->data() << std::endl;
        std::cout << std::endl;
      }
    }

    /*! Create a face f that matches the overlapping region between f1 and f2. */
    virtual void create_face(Face_const_handle f1, Face_const_handle f2,
                             Face_handle f) const
    {
      if (m_verbose_level > 2) {
        std::cout << "  f1: " << "Inner(" << count_inner(f1) << ")"
                  << ", Outer(" << count_outer(f1) << ")"
                  << ", " << f1->data() << std::endl;
        std::cout << "  f2: " << "Inner(" << count_inner(f2) << ")"
                  << ", Outer(" << count_outer(f2) << ")"
                  << ", " << f2->data() << std::endl;
      }
      f->set_data(f1->data() + f2->data());
      if (m_verbose_level > 2) {
        std::cout << "  f: " << "Inner(" << count_inner(f) << ")"
                  << ", Outer(" << count_outer(f) << ")"
                  << ", " << f->data() << std::endl;
        std::cout << std::endl;
      }
    }
  };

  // A functor that tests whether two points are equal
  class Point_equal {
  public:
    Point_equal(const Geom_traits& traits) : m_traits(traits) {}
    bool operator()(const Point_2& p1, const Point_2& p2)
    { return (m_traits.equal_2_object()(p1, p2)); }

  private:
    const Geom_traits& m_traits;
  };

  // Maps
  struct Less_than_handle {
    template <typename Type>
    bool operator()(Type s1, Type s2) const { return (&(*s1) < &(*s2)); }
  };

  typedef std::map<Vertex_const_handle, Vertex_const_handle,
                            Less_than_handle>           Vertex_map;
  typedef std::map<Halfedge_const_handle, Halfedge_const_handle,
                            Less_than_handle>           Halfedge_map;
  typedef std::map<Face_const_handle, Face_const_handle,
                            Less_than_handle>           Face_map;

  // A functor that tests whether two c-monotone curves are equal
  class Curve_equal {
  public:
    Curve_equal(const Geom_traits& traits) : m_traits(traits) {}
    bool operator()(const Xcurve_2& c1, const Xcurve_2& c2)
    { return (m_traits.equal_2_object()(c1, c2) && (c1.data() == c2.data())); }

  private:
    const Geom_traits& m_traits;
  };

  // Read input data for a single arrangement
  template <typename XcurveOutputIterator, typename PointOutputIterator>
  bool read_arr(std::istream& in, Arrangement& arr,
                XcurveOutputIterator xcurves,
                PointOutputIterator isolated_points);

  // Construct an arrangement
  template <typename Curve_iterator, typename Point_iterator>
  void construct_arr(Arrangement& arr,
                     Curve_iterator xcurves_begin, Curve_iterator xcurves_end,
                     Point_iterator points_begin, Point_iterator points_end);

  // Initialize an arrangement
  bool init_arr(Arrangement& arr);

  // A predicate that verifies the results
  bool is_interior(Vertex_const_handle vh);
  // bool are_same_results();

  // Checks whether two ccb's are equivalent
  bool equivalent_ccb(Ccb_halfedge_const_circulator ccb1,
                      Ccb_halfedge_const_circulator ccb2,
                      const Halfedge_map& halfedge_map);

  // Check whether two faces are equivalent
  bool equivalent_face(Face_const_handle fit1, Face_const_handle fit2,
                       const Halfedge_map& halfedge_map);

  // Check whether two arrangement are equivalent
  bool equivalent_arr(const Arrangement& arr1, const Arrangement& arr2);
};

/*! Constructor */
template <typename T_Geom_traits, typename T_Topol_traits>
Overlay_test<T_Geom_traits, T_Topol_traits>::
Overlay_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_verbose_level(0)
{}

/*! Clear the data structures */
template<class T_Geom_traits, typename T_Topol_traits>
void Overlay_test<T_Geom_traits, T_Topol_traits>::clear()
{
  m_arr1.clear();
  m_xcurves1.clear();
  m_isolated_points1.clear();

  m_arr2.clear();
  m_xcurves2.clear();
  m_isolated_points2.clear();

  m_arr.clear();
  m_filename.clear();
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::
is_interior(Vertex_const_handle vh)
{
  return ((vh->parameter_space_in_x() == CGAL::ARR_INTERIOR) &&
          (vh->parameter_space_in_y() == CGAL::ARR_INTERIOR));
}

#if 0
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::are_same_results()
{
  if (m_verbose_level > 1) {
    std::cout << "# vertices, edge, faces obtained: ("
              << m_arr->number_of_vertices() << ","
              << m_arr->number_of_edges() << ","
              << m_arr->number_of_faces() << ")"
              << ", expected: ("
              << m_num_vertices << ","
              << m_num_edges << ","
              << m_num_faces << ")" << std::endl;
  }

  if (m_arr->number_of_vertices() != m_num_vertices) return false;
  if (m_arr->number_of_edges() != m_num_edges) return false;
  if (m_arr->number_of_faces() != m_num_faces) return false;

  Point_container points_res(m_num_vertices);
  typename Point_container::iterator pit = points_res.begin();
  Vertex_const_iterator vit;
  for (vit = m_arr->vertices_begin(); vit != m_arr->vertices_end(); ++vit) {
    if (is_interior(vit))
      *pit++ = vit->point();
  }
  Point_compare<Geom_traits> pt_compare(m_geom_traits);
  std::sort(points_res.begin(), pit, pt_compare);

  if (m_verbose_level > 2) {
    std::copy(points_res.begin(), pit,
              std::ostream_iterator<Point_2>(std::cout, "\n"));
  }

  Point_equal point_eq(m_geom_traits);
  if (! std::equal(points_res.begin(), pit, m_points.begin(), point_eq))
    return false;

  std::vector<Xcurve_2> curves_res(m_arr->number_of_edges());
  typename Xcurve_container::iterator xcit = curves_res.begin();

  Edge_const_iterator eit;
  for (eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit) {
    if (is_interior(eit->source()) && is_interior(eit->target()))
      *xcit++ = eit->curve();
  }
  Curve_compare<Geom_traits> curve_compare(m_geom_traits);
  std::sort(curves_res.begin(), xcit, curve_compare);

  if (m_verbose_level > 2) {
    std::copy(curves_res.begin(), xcit,
              std::ostream_iterator<Xcurve_2>(std::cout, "\n"));
  }

  Curve_equal curve_eq(m_geom_traits);
  if (! std::equal(curves_res.begin(), xcit, m_xcurves.begin(), curve_eq))
    return false;

  return true;
}
#endif

// Read input data for a single arrangement
template <typename T_Geom_traits, typename T_Topol_traits>
template <typename XcurveOutputIterator, typename PointOutputIterator>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::
read_arr(std::istream& in, Arrangement& /* arr */,
         XcurveOutputIterator xcurves, PointOutputIterator isolated_points)
{
  unsigned int i;

  unsigned int num_of_curves;
  in >> num_of_curves;
  for (i = 0; i < num_of_curves; ++i) {
    Base_xcurve_2 base_xcv;
    if (!this->read_xcurve(in, base_xcv)) return false;
    Xcurve_2 xcv(base_xcv, 1);
    *xcurves++ = xcv;
  }

  unsigned int num_of_isolated_points;
  in >> num_of_isolated_points;
  for (i = 0; i < num_of_isolated_points; ++i) {
    Point_2 point;
    if (!this->read_point(in, point)) return false;
    *isolated_points++ = point;
  }

  return true;
}

template <typename T_Geom_traits, typename T_Topol_traits>
template <typename Curve_iterator, typename Point_iterator>
void Overlay_test<T_Geom_traits, T_Topol_traits>::
construct_arr(Arrangement& arr,
              Curve_iterator xcurves_begin, Curve_iterator xcurves_end,
              Point_iterator points_begin, Point_iterator points_end)
{
  typedef T_Geom_traits                 Geom_traits;

#if 1
  // Insert the curves incrementally.<
  Curve_iterator cit;
  for (cit = xcurves_begin; cit != xcurves_end; ++cit) {
    if (m_verbose_level > 2) std::cout << "inserting " << *cit << " ... ";
    std::cout.flush();
    CGAL::insert(arr, *cit);
    if (m_verbose_level > 2) std::cout << "inserted" << std::endl;
  }
#else
  // Insert the curves aggregately.
  if (m_verbose_level > 2) std::cout << "inserting x-monotone curves"
                                     << " ... ";
  std::cout.flush();
  CGAL::insert(arr, xcurves_begin, xcurves_end);
  if (m_verbose_level > 2) std::cout << "inserted" << std::endl;
#endif

  // Insert the isolated points.
  if (m_verbose_level > 2) std::cout << "inserting isolated vertices"
                                     << " ... ";
  Point_iterator pit;
  for (pit = points_begin; pit != points_end; ++pit) {
    Point_2 point(*pit);
    CGAL::insert_point(arr, point);
  }
  if (m_verbose_level > 2) std::cout << "inserted" << std::endl;

  #if 1
  // Merge mergeable vertices
  Vertex_iterator vit;
  Vertex_iterator vit_next;
  const Geom_traits* traits = arr.geometry_traits();
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); vit = vit_next) {
    vit_next = vit;
    ++vit_next;
    if (vit->degree() != 2) continue;
    Halfedge_around_vertex_circulator eit = vit->incident_halfedges();
    if (traits->are_mergeable_2_object()(eit->curve(), eit->next()->curve())) {
      // std::cout << "merging: " << vit->point() << std::endl;
      Xcurve_2 xcv;
      traits->merge_2_object()(eit->curve(), eit->next()->curve(), xcv);
      arr.merge_edge(eit, eit->next(), xcv);
    }
  }
  #endif
}

// Read input data for a single arrangement
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::init_arr(Arrangement& arr)
{
  // Initialize the data of the halfedges of the arrangement.
  Halfedge_iterator heit;
  for (heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit)
    heit->set_data(heit->curve().data() *
                   ((heit->direction() == CGAL::ARR_LEFT_TO_RIGHT) ? 1 : 2));

  // Initialize the data of the faces of the arrangement.
  Face_iterator fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    unsigned int count = 0;

    // Outer ccb
    Outer_ccb_iterator ocit;
    for (ocit = fit->outer_ccbs_begin(); ocit != fit->outer_ccbs_end(); ++ocit)
    {
      Ccb_halfedge_circulator curr = *ocit;
      do count += curr->data() * 2;
      while (++curr != *ocit);
    }

    // Inner ccbs
    Inner_ccb_iterator icit;
    for (icit = fit->inner_ccbs_begin(); icit != fit->inner_ccbs_end(); ++icit)
    {
      Ccb_halfedge_circulator curr = *icit;
      do count += curr->data();
      while (++curr != *icit);
    }

    fit->set_data(count);
  }

  // Initialize the data of the vertices of the arrangement.
  Vertex_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    unsigned int count = 0;
    if (vit->is_isolated()) count = vit->face()->data();
    else {
      Halfedge_around_vertex_const_circulator curr = vit->incident_halfedges();
      do count += curr->data();
      while (++curr != vit->incident_halfedges());
    }
    vit->set_data(count);
  }

  if (m_verbose_level > 0) std::cout << "Arrangement Input: " << std::endl;

  if (m_verbose_level > 2) {
    std::cout << "Vertex Data: " << std::endl;
    Vertex_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
      std::cout << vit->point() << " " << vit->data() << std::endl;
  }

  if (m_verbose_level > 1) {
    std::cout << "Halfedge Data: " << std::endl;
    Halfedge_iterator heit;
    for (heit = arr.halfedges_begin(); heit != arr.halfedges_end(); ++heit)
      std::cout << heit->source()->point() << " "
                << heit->target()->point() << " " << heit->data()
                << std::endl;
  }

  if (m_verbose_level > 0) {
    std::cout << "Face Data: " << std::endl;
    Face_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
      std::cout << fit->data() << std::endl;
  }

  return true;
}

// Initialize the test
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::init()
{
  std::ifstream p_stream(m_filename.c_str());
  if (!p_stream.is_open()) {
    std::cerr << "Cannot open file " << m_filename << "!" << std::endl;
    return false;
  }

  // 1st arrangement.
  read_arr(p_stream, m_arr1,
           std::back_inserter(m_xcurves1),
           std::back_inserter(m_isolated_points1));
  construct_arr(m_arr1, m_xcurves1.begin(), m_xcurves1.end(),
                m_isolated_points1.begin(), m_isolated_points1.end());
  init_arr(m_arr1);

  // 2nd arrangement.
  read_arr(p_stream, m_arr2,
           std::back_inserter(m_xcurves2),
           std::back_inserter(m_isolated_points2));
  construct_arr(m_arr2, m_xcurves2.begin(), m_xcurves2.end(),
                m_isolated_points2.begin(), m_isolated_points2.end());
  init_arr(m_arr2);

  // Consume eol
  int c;
  while ((c = p_stream.get()) != '\n');

  // Expected arrangement.
  Formatter formatter;
  CGAL::read(m_arr, p_stream, formatter);

  p_stream.close();

  return true;
}

// Check whether two ccb's are equivalent
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::
equivalent_ccb(Ccb_halfedge_const_circulator ccb1,
               Ccb_halfedge_const_circulator ccb2,
               const Halfedge_map& halfedge_map)
{
  // Find a matching starting point.
  typename Halfedge_map::const_iterator it = halfedge_map.find(ccb1);
  if (it == halfedge_map.end()) return false;
  Ccb_halfedge_const_circulator curr2 = ccb2;
  do if ((*it).second == curr2) break;
  while (++curr2 != ccb2);
  if ((*it).second != curr2) return false;
  ccb2 = curr2;

  // Match the rest
  Ccb_halfedge_const_circulator curr1 = ccb1;
  do {
    typename Halfedge_map::const_iterator it = halfedge_map.find(curr1);
    if ((it == halfedge_map.end()) || ((*it).second != curr2++)) return false;
  } while (++curr1 != ccb1);
  if (curr2 != ccb2) return false;

  return true;
}

// Check whether two faces are equivalent
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::
equivalent_face(Face_const_handle fit1, Face_const_handle fit2,
                const Halfedge_map& halfedge_map)
{
  if (m_verbose_level > 0) {
    std::cout << "face 1: "
              << fit1->number_of_outer_ccbs() << ","
              << fit1->number_of_inner_ccbs() << ","
              << fit1->data() << std::endl;
    std::cout << "face 2: "
              << fit2->number_of_outer_ccbs() << ","
              << fit2->number_of_inner_ccbs() << ","
              << fit2->data() << std::endl;
  }

  if (fit1->number_of_outer_ccbs() != fit2->number_of_outer_ccbs())
    return false;

  if (fit1->number_of_inner_ccbs() != fit2->number_of_inner_ccbs())
    return false;

  // Outer ccb
  Outer_ccb_const_iterator ocit1;
  Outer_ccb_const_iterator ocit2;
  for (ocit1 = fit1->outer_ccbs_begin(); ocit1 != fit1->outer_ccbs_end();
       ++ocit1)
  {
    Ccb_halfedge_const_circulator ccb1 = *ocit1;
    bool found = false;
    for (ocit2 = fit2->outer_ccbs_begin(); ocit2 != fit2->outer_ccbs_end();
         ++ocit2)
    {
      Ccb_halfedge_const_circulator ccb2 = *ocit2;
      if (equivalent_ccb(ccb1, ccb2, halfedge_map)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  // Inner ccb
  Outer_ccb_const_iterator icit1;
  Outer_ccb_const_iterator icit2;
  for (icit1 = fit1->inner_ccbs_begin(); icit1 != fit1->inner_ccbs_end();
       ++icit1)
  {
    Ccb_halfedge_const_circulator ccb1 = *icit1;
    bool found = false;
    for (icit2 = fit2->inner_ccbs_begin(); icit2 != fit2->inner_ccbs_end();
         ++icit2)
    {
      Ccb_halfedge_const_circulator ccb2 = *icit2;
      if (equivalent_ccb(ccb1, ccb2, halfedge_map)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  return true;
}

// Check whether two arrangement are equivalent
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::
equivalent_arr(const Arrangement& arr1, const Arrangement& arr2)
{
  if (arr1.number_of_vertices() != arr2.number_of_vertices()) return false;
  if (arr1.number_of_halfedges() != arr2.number_of_halfedges()) return false;
  if (arr1.number_of_faces() != arr2.number_of_faces()) return false;

  const Geom_traits* traits = arr1.geometry_traits();
  typename Geom_traits::Equal_2 equal = traits->equal_2_object();

  // Compare the vertices
  const Vertex_const_iterator invalid_vit;
  Vertex_map vertex_map;
  Vertex_const_iterator vit1;
  for (vit1 = arr1.vertices_begin(); vit1 != arr1.vertices_end(); ++vit1) {
    const Point_2& p1 = vit1->point();
    Vertex_const_iterator vit2;
    bool found = false;
    for (vit2 = arr2.vertices_begin(); vit2 != arr2.vertices_end(); ++vit2) {
      const Point_2& p2 = vit2->point();
      if (equal(p1, p2)) {
        if (vertex_map[vit1] != invalid_vit) {
          std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                    << ") has been mapped already!"
                    << std::endl;
          return false;
        }
        vertex_map[vit1] = vit2;
        found = true;

        if (vit1->data() != vit2->data()) {
          std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                    << ") data does not match (" << vit2->data()
                    << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The vertex ((" << p1 << "), " << vit1->data()
                << ") was not found!" << std::endl;
       return false;
    }
  }

  // Compare the halfedges.
  const Halfedge_const_iterator invalid_heit;
  Halfedge_map halfedge_map;
  Halfedge_const_iterator heit1;
  for (heit1 = arr1.halfedges_begin(); heit1 != arr1.halfedges_end(); ++heit1) {
    const Xcurve_2& xcv1 = heit1->curve();
    Halfedge_const_iterator heit2;
    bool found = false;
    for (heit2 = arr2.halfedges_begin(); heit2 != arr2.halfedges_end(); ++heit2)
    {
      const Xcurve_2& xcv2 = heit2->curve();
      if ((vertex_map[heit1->source()] == heit2->source()) &&
          (vertex_map[heit1->target()] == heit2->target()) &&
          equal(xcv1, xcv2))
      {
        if (halfedge_map[heit1] != invalid_heit) {
          std::cerr << "The halfedge ((" << heit1->source()->point()
                    << " => " << heit1->target()->point() << "), "
                    << heit1->data() << ") has been mapped already!"
                    << std::endl;
          return false;
        }
        halfedge_map[heit1] = heit2;
        found = true;

        if (heit1->data() != heit2->data()) {
          std::cerr << "The halfedge ((" << heit1->source()->point()
                    << " => " << heit1->target()->point() << "), "
                    << heit1->data() << ") data does not match ("
                    << heit2->data() << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The halfedge ((" << heit1->source()->point()
                << " => " << heit1->target()->point() << "), "
                << heit1->data() << ") was not found!"
                << std::endl;
      return false;
    }
  }

  // Compare the faces.
  const Face_const_iterator invalid_fit;
  Face_map face_map;
  Face_const_iterator fit1;
  for (fit1 = arr1.faces_begin(); fit1 != arr1.faces_end(); ++fit1) {
    Face_const_iterator fit2;
    bool found = false;
    for (fit2 = arr2.faces_begin(); fit2 != arr2.faces_end(); ++fit2) {
      if (equivalent_face(fit1, fit2, halfedge_map)) {
        if (face_map[fit1] != invalid_fit) {
          std::cerr << "The face (" << fit1->data()
                    << ") has been mapped already!" << std::endl;
          return false;
        }
        face_map[fit1] = fit2;
        found = true;

        if (fit1->data() != fit2->data()) {
          std::cerr << "The face (" << fit1->data()
                    << ") data does not match ("
                    << fit2->data() << ")!" << std::endl;
          return false;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "The face (" << fit1->data()
                << ") was not found!" << std::endl;
      return false;
    }
  }

  return true;
}

// Perform the test
template <typename T_Geom_traits, typename T_Topol_traits>
bool Overlay_test<T_Geom_traits, T_Topol_traits>::perform()
{
  // Overlay the input arrangements:
  Arrangement arr;
  Overlay_traits overlay_traits(m_verbose_level);
  // Formatter formatter;
  // CGAL::write(m_arr2, std::cout, formatter);

  CGAL::overlay(m_arr1, m_arr2, arr, overlay_traits);

  // Generate the output for debugging purposes
  // Formatter formatter;
  // CGAL::write(arr, std::cout, formatter);

  // Verify the resulting arrangement:
  if (!equivalent_arr(arr, m_arr)) {
    arr.clear();
    return false;
  }

  arr.clear();
  return true;
}

#endif
