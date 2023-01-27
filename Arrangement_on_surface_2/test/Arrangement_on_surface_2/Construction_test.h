#ifndef CGAL_CONSTRUCTION_TEST_H
#define CGAL_CONSTRUCTION_TEST_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


#include <CGAL/Timer.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arrangement_on_surface_2.h>

#include "utils.h"
#include "IO_base_test.h"

#include <CGAL/disable_warnings.h>

/*! Construction test */
template <typename T_Geom_traits, typename T_Topol_traits>
class Construction_test :
  public IO_base_test<typename T_Geom_traits::Base_traits_2>
{
public:
  typedef T_Geom_traits                                         Geom_traits;
  typedef T_Topol_traits                                        Topol_traits;
  typedef IO_base_test<typename Geom_traits::Base_traits_2>     Base;

  typedef typename Geom_traits::Base_traits_2           Base_geom_traits;
  typedef typename Base_geom_traits::Point_2            Base_point_2;
  typedef typename Base_geom_traits::Curve_2            Base_curve_2;
  typedef typename Base_geom_traits::X_monotone_curve_2 Base_x_monotone_curve_2;

  typedef typename Geom_traits::Point_2                 Point_2;
  typedef typename Geom_traits::Curve_2                 Curve_2;
  typedef typename Geom_traits::X_monotone_curve_2      X_monotone_curve_2;

  typedef typename std::vector<Curve_2>                 Curve_container;
  typedef typename std::vector<Point_2>                 Point_container;
  typedef typename std::vector<X_monotone_curve_2>      Xcurve_container;

  typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                        Arrangement;
  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;

  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;

  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement::Edge_iterator           Edge_iterator;
  typedef typename Arrangement::Face_iterator           Face_iterator;

  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;

private:
  /*! The geometry traits */
  const Geom_traits& m_geom_traits;

  /*! The arrangement */
  Arrangement* m_arr;

  /*! Verbosity */
  int m_verbose_level;

  typename Arrangement::Size m_num_vertices;
  typename Arrangement::Size m_num_edges;
  typename Arrangement::Size m_num_faces;

  /*! The input data file of points*/
  std::string m_filename;

  Curve_container m_curves;
  Point_container m_isolated_points;
  Point_container m_points;
  Xcurve_container m_xcurves;

  /*! The number of inner and outer ccbs of faces. */
  typedef std::pair<typename Arrangement::Size, typename Arrangement::Size>
                                                        Face_ccbs;
  typedef std::vector<Face_ccbs>                        Face_ccbs_vector;
  typedef typename Face_ccbs_vector::const_iterator     Face_ccbs_const_iter;
  typedef typename Face_ccbs_vector::iterator           Face_ccbs_iter;
  Face_ccbs_vector m_faces;

public:
  /*! Constructor */
  Construction_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Construction_test() { clear(); }

  void set_verbose_level(int verbose_level) { m_verbose_level = verbose_level; }
  void set_filename(const char* filename) { m_filename.assign(filename); }

  bool allocate_arrangement();

  void deallocate_arrangement();

  /*! Initialize the test */
  virtual bool init();

  /*! Perform the test */
  virtual bool perform();

  /*! Clear the data structures */
  virtual void clear();

protected:
  // A functor that tests whether two points are equal
  class Point_equal {
  public:
    Point_equal(const Geom_traits& traits) : m_traits(traits) {}
    bool operator()(const Point_2& p1, const Point_2& p2)
    { return (m_traits.equal_2_object()(p1, p2)); }

  private:
    const Geom_traits& m_traits;
  };

  // A functor that tests whether two c-monotone curves are equal
  class Curve_equal {
  public:
    Curve_equal(const Geom_traits& traits) : m_traits(traits) {}
    bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
    { return (m_traits.equal_2_object()(c1, c2) && (c1.data() == c2.data())); }

  private:
    const Geom_traits& m_traits;
  };

  bool test1();
  bool test2();
  bool test3();
  bool test4();
  bool test5();
  bool test6();
  bool test7();
  bool test8();
  bool test9();

  // A predicate that verifies the results
  bool is_interior(Vertex_const_handle vh);
  bool are_same_results();

  typedef typename Geom_traits::Left_side_category      Left_side_category;
  typedef typename Geom_traits::Bottom_side_category    Bottom_side_category;
  typedef typename Geom_traits::Top_side_category       Top_side_category;
  typedef typename Geom_traits::Right_side_category     Right_side_category;
  typedef typename
  CGAL::Arr_all_sides_not_open_category<Left_side_category,
                                        Bottom_side_category,
                                        Top_side_category,
                                        Right_side_category>::result
    All_sides_not_open_category;

  bool is_open() const
  { return is_open(All_sides_not_open_category()); }

  bool is_open(CGAL::Arr_all_sides_not_open_tag) const { return false; }

  bool is_open(CGAL::Arr_not_all_sides_not_open_tag) const { return true; }
};

/*! Constructor */
template <typename T_Geom_traits, typename T_Topol_traits>
Construction_test<T_Geom_traits, T_Topol_traits>::
Construction_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_arr(nullptr),
  m_verbose_level(0)
{}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::allocate_arrangement()
{
  if (!(m_arr = new Arrangement(&m_geom_traits))) return false;
  return true;
}

template <typename T_Geom_traits, typename T_Topol_traits>
void Construction_test<T_Geom_traits, T_Topol_traits>::deallocate_arrangement()
{
  if (m_arr) {
    delete m_arr;
    m_arr = nullptr;
  }
}

/*! Clear the data structures */
template<class T_Geom_traits, typename T_Topol_traits>
void Construction_test<T_Geom_traits, T_Topol_traits>::clear()
{
  deallocate_arrangement();

  m_curves.clear();
  m_isolated_points.clear();
  m_filename.clear();
  m_points.clear();
  m_faces.clear();
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::
is_interior(Vertex_const_handle vh)
{
  return ((vh->parameter_space_in_x() == CGAL::ARR_INTERIOR) &&
          (vh->parameter_space_in_y() == CGAL::ARR_INTERIOR));
}

template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::are_same_results()
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

  if ((m_arr->number_of_vertices() != m_num_vertices) ||
      (m_arr->number_of_edges() != m_num_edges) ||
      (m_arr->number_of_faces() != m_num_faces)) {
    std::cout << "# vertices, edge, faces obtained: ("
              << m_arr->number_of_vertices() << ","
              << m_arr->number_of_edges() << ","
              << m_arr->number_of_faces() << ")"
              << ", expected: ("
              << m_num_vertices << ","
              << m_num_edges << ","
              << m_num_faces << ")" << std::endl;
    return false;
  }

  Point_container points_res(m_num_vertices);
  typename Point_container::iterator pit = points_res.begin();
  Vertex_const_iterator vit;
  for (vit = m_arr->vertices_begin(); vit != m_arr->vertices_end(); ++vit) {
    if (! vit->is_at_open_boundary()) *pit++ = vit->point();
  }
  Point_compare<Geom_traits> pt_compare(m_geom_traits);
  std::sort(points_res.begin(), pit, pt_compare);

  if (m_verbose_level > 2) {
    std::copy(points_res.begin(), points_res.end(),
              std::ostream_iterator<Point_2>(std::cout, "\n"));
  }

  Point_equal point_eq(m_geom_traits);
  if (! std::equal(points_res.begin(), pit, m_points.begin(), point_eq)) {
    std::cout << "Expected: " << std::endl;
    std::copy(m_points.begin(), m_points.end(),
              std::ostream_iterator<Point_2>(std::cout, "\n"));
    std::cout << "Obtained: " << std::endl;
    std::copy(points_res.begin(), points_res.end(),
              std::ostream_iterator<Point_2>(std::cout, "\n"));
    return false;
  }

  std::vector<X_monotone_curve_2> curves_res(m_arr->number_of_edges());
  typename Xcurve_container::iterator xcit = curves_res.begin();

  Edge_const_iterator eit;
  for (eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit)
    *xcit++ = eit->curve();

  Curve_compare<Geom_traits> curve_compare(m_geom_traits);
  std::sort(curves_res.begin(), xcit, curve_compare);

  if (m_verbose_level > 2) {
    for (typename Xcurve_container::iterator it = curves_res.begin();
         it != curves_res.end(); ++it)
      std::cout << *it << " " << it->data() << std::endl;
  }

  Curve_equal curve_eq(m_geom_traits);
  if (! std::equal(curves_res.begin(), xcit, m_xcurves.begin(), curve_eq)) {
    std::cout << "Expected: " << std::endl;
    for (typename Xcurve_container::iterator it = m_xcurves.begin();
         it != m_xcurves.end(); ++it)
      std::cout << *it << " " << it->data() << std::endl;
    std::cout << "Obtained: " << std::endl;
    for (typename Xcurve_container::iterator it = curves_res.begin();
         it != curves_res.end(); ++it)
      std::cout << *it << " " << it->data() << std::endl;
    return false;
  }

  if (m_arr->number_of_faces() == 1) {
    Face_const_iterator fit = m_arr->faces_begin();
    if (m_verbose_level > 1)
      std::cout << "Face: # inner " << fit->number_of_inner_ccbs()
                  << ", # outer: " << fit->number_of_outer_ccbs()
                  << std::endl;
    if (is_open()) {
      if (fit->number_of_outer_ccbs() != 1) return false;
    }
    else {
      if (fit->number_of_outer_ccbs() != 0) return false;
    }
  }
  else {
    if (m_faces.empty()) {
      m_faces.resize(m_num_faces);
      Face_ccbs_iter cit = m_faces.begin();
      Face_const_iterator fit;
      for (fit = m_arr->faces_begin(); fit != m_arr->faces_end(); ++fit) {
        if (m_verbose_level > 1)
          std::cout << "Face: # inner " << fit->number_of_inner_ccbs()
                    << ", # outer: " << fit->number_of_outer_ccbs()
                    << std::endl;
        *cit++ = Face_ccbs(fit->number_of_inner_ccbs(),
                           fit->number_of_outer_ccbs());
      }
    }
    else {
      Face_ccbs_vector faces = m_faces;
      Face_const_iterator fit;
      for (fit = m_arr->faces_begin(); fit != m_arr->faces_end(); ++fit) {
        if (m_verbose_level > 1)
          std::cout << "Face: # inner " << fit->number_of_inner_ccbs()
                    << ", # outer: " << fit->number_of_outer_ccbs()
                    << std::endl;
        Face_ccbs face_ccbs(fit->number_of_inner_ccbs(),
                            fit->number_of_outer_ccbs());
        Face_ccbs_iter cit = std::find(faces.begin(), faces.end(), face_ccbs);
        if (cit == faces.end()) return false;
        *cit = Face_ccbs(0, 0);
      }
      faces.clear();
    }
  }

  return true;
}

// Perform the test
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::init()
{
  unsigned int i;

  // Read the input curves and isolated vertices.
  std::ifstream p_stream(m_filename.c_str());
  if (!p_stream.is_open()) {
    std::cerr << "Cannot open file " << m_filename << "!" << std::endl;
    return false;
  }

  unsigned int num_of_curves;
  p_stream >> num_of_curves;
  for (i = 0; i < num_of_curves; ++i) {
    Base_curve_2 base_cv;
    this->read_curve(p_stream, base_cv);
    m_curves.push_back(Curve_2(base_cv, 1));
  }

  unsigned int num_of_isolated_points;
  p_stream >> num_of_isolated_points;
  for (i = 0; i < num_of_isolated_points; ++i) {
    Point_2 point;
    this->read_point(p_stream, point);
    m_isolated_points.push_back(point);
  }

  // Read the points and curves that correspond to the arrangement vertices
  // and edges, respectively.
  p_stream >> m_num_vertices >> m_num_edges >> m_num_faces;

  m_points.resize(m_num_vertices);
  for (i = 0; i < m_num_vertices; ++i) {
    Point_2 point;
    this->read_point(p_stream, point);
    m_points[i] = point;
  }

  m_xcurves.resize(m_num_edges);
  for (i = 0; i < m_num_edges; ++i) {
    Base_x_monotone_curve_2 base_xcv;
    this->read_xcurve(p_stream, base_xcv);
    unsigned int k;
    p_stream >> k;
    m_xcurves[i] = X_monotone_curve_2(base_xcv, k);
  }
  p_stream.close();

  if (! allocate_arrangement()) return false;

  return true;
}

// Test incremental construction.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test1()
{
  typename Curve_container::const_iterator cit;
  for (cit = m_curves.begin(); cit != m_curves.end(); ++cit)
    CGAL::insert(*m_arr, *cit);
  typename Point_container::const_iterator pit;
  for (pit = m_isolated_points.begin(); pit != m_isolated_points.end(); ++pit)
    CGAL::insert_point(*m_arr, *pit);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cerr << "ERROR : (1) The incremental insertion test failed (invalid)."
              << std::endl;
    return false;
  }
#endif
  if (! are_same_results()) {
    std::cerr << "ERROR : (1) The incremental insertion test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(1) Passed incremental insertion." << std::endl;
  m_arr->clear();
  return true;
}

// Test aggregate construction.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test2()
{
  CGAL::insert(*m_arr, m_curves.begin(), m_curves.end());
  // When creating insert_points, this call should be fixed to insert_points.
  typename Point_container::const_iterator pit;
  for (pit = m_isolated_points.begin(); pit != m_isolated_points.end(); ++pit)
    CGAL::insert_point(*m_arr, *pit);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cerr << "ERROR : (2) The aggregated construction test failed (invalid)."
              << std::endl;
    return false;
  }
#endif
  if (! are_same_results()) {
    std::cerr << "ERROR : (2) The aggregated construction test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(2) Passed aggregated construction." << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of half of the curves aggregatley followed by the
// insertion of the rest aggregatley (test the insertion visitor).
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test3()
{
  size_t n = m_curves.size() / 2;
  CGAL::insert(*m_arr, m_curves.begin(), m_curves.begin() + n);
  CGAL::insert(*m_arr, m_curves.begin() + n, m_curves.end());
  // When creating insert_points, this call should be fixed to insert_points.
  typename Point_container::const_iterator pit;
  for (pit = m_isolated_points.begin(); pit != m_isolated_points.end(); ++pit)
  CGAL::insert_point(*m_arr, *pit);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cerr << "ERROR : (3) The aggregated insertion test failed (invalid)."
              << std::endl;
    return false;
  }
#endif
  if (! are_same_results()) {
    std::cerr << "ERROR : (3) The aggregated insertion test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(3) Passed aggregated insertion." << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of the interior-disjoint subcurves incrementally.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test4()
{
  size_t i;
  for (i = 0; i < m_xcurves.size(); ++i)
    CGAL::insert(*m_arr, m_xcurves[i]);
  for (i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cerr << "ERROR : (4) The incremental x-monotone test failed (invalid)."
              << std::endl;
    return false;
  }
#endif
  if (! are_same_results()) {
    std::cerr << "ERROR : (4) The incremental x-monotone test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(4) Passed incremental x-monotone insertion." << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of the interior-disjoint subcurves aggregatley.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test5()
{
  CGAL::insert(*m_arr, m_xcurves.begin(), m_xcurves.end());
  for (size_t i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cerr << "ERROR : (5) The aggregated x-monotone construction test failed (invalid)."
              << std::endl;
  }
#endif
  if (! are_same_results()) {
    std::cerr << "ERROR : (5) The aggregated x-monotone construction test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(5) Passed aggregated x-monotone construction." << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of half of the interior-disjoint subcurves aggregatley
// followed by the insertion of the rest aggregatley with
// insert_x_monotone_curves(test the addition visitor).
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test6()
{
  CGAL::insert(*m_arr, m_xcurves.begin(),
               m_xcurves.begin() + (m_num_edges/2));
  CGAL::insert(*m_arr, m_xcurves.begin() + (m_num_edges/2),
               m_xcurves.end());
  for (size_t i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cout << "ERROR : (6) The aggregated x-monotone insertion test failed (invalid)."
              << std::endl;
  }
#endif
  if (! are_same_results()) {
    std::cout << "ERROR : (6) The aggregated x-monotone insertion test failed."
              << std::endl;
     return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(6) Passed aggregated x-monotone insertion." << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of the disjoint subcurves incrementally with
// insert_non_intersecting_curve.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test7()
{
  size_t i;
  for (i = 0; i < m_xcurves.size(); ++i)
    CGAL::insert_non_intersecting_curve(*m_arr, m_xcurves[i]);
  for (i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cout << "ERROR : (7) The incremental non-intersecting test failed (invalid)."
              << std::endl;
  }
#endif
  if (! are_same_results()) {
    std::cout << "ERROR : (7) The incremental non-intersecting test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(7) Passed incremental non-intersecting insertion."
              << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of the interior-disjoint subcurves aggregatley with
// insert_non_intersecting_curves.
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test8()
{
  CGAL::insert_non_intersecting_curves(*m_arr,
                                       m_xcurves.begin(), m_xcurves.end());
  for (size_t i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cout
      << "ERROR : (8) The aggregated non-intersecting construction test failed (invalid)."
      << std::endl;
  }
#endif
  if (! are_same_results()) {
    std::cout
      << "ERROR : (8) The aggregated non-intersecting construction test failed."
      << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(8) Passed aggregated non-intersecting construction."
              << std::endl;
  m_arr->clear();
  return true;
}

// Test the insertion of half of the disjoint subcurves aggregatley followed
// by the insertion of the rest aggregatley with
// insert_non_intersecting_curves (test the addition visitor).
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::test9()
{
  size_t n = m_xcurves.size() / 2;
  CGAL::insert_non_intersecting_curves(*m_arr,
                                       m_xcurves.begin(), m_xcurves.begin() + n);
  CGAL::insert_non_intersecting_curves(*m_arr,
                                       m_xcurves.begin() + n, m_xcurves.end());
  for (size_t i = 0; i < m_isolated_points.size(); ++i)
    CGAL::insert_point(*m_arr, m_isolated_points[i]);
#if TEST_TOPOL_TRAITS != SPHERICAL_TOPOL_TRAITS
  if (! CGAL::is_valid(*m_arr)) {
    std::cout <<"ERROR : (9) The aggregated non-intersecting insertion test failed (invalid)."
              << std::endl;
  }
#endif
  if (! are_same_results()) {
    std::cout <<"ERROR : (9) The aggregated non-intersecting insertion test failed."
              << std::endl;
    return false;
  }

  if (m_verbose_level > 0)
    std::cout << "(9) Passed aggregated non-intersecting insertion."
              << std::endl;
  m_arr->clear();
  return true;
}

// Perform the test
template <typename T_Geom_traits, typename T_Topol_traits>
bool Construction_test<T_Geom_traits, T_Topol_traits>::perform()
{
  if (! test1()) return false;
  if (! test2()) return false;
  if (! test3()) return false;
  if (! test4()) return false;
  if (! test5()) return false;
  if (! test6()) return false;
  if (! test7()) return false;
  if (! test8()) return false;
  if (! test9()) return false;
  return true;
}

#include <CGAL/enable_warnings.h>

#endif
