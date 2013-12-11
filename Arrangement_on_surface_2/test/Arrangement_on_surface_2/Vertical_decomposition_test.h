#ifndef CGAL_BATCHED_POINT_LOCATION_TEST_H
#define CGAL_BATCHED_POINT_LOCATION_TEST_H

#include <vector>
#include <list>
#include <utility>

#include <CGAL/basic.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Arr_point_location_result.h>

#include "IO_test.h"

/*! Point location test */
template <typename GeomTraits_T, typename TopolTraits_T>
class Vertical_decomposition_test : public IO_test<GeomTraits_T> {
public:
  typedef GeomTraits_T                                    Geom_traits;
  typedef TopolTraits_T                                   Topol_traits;

private:
  typedef IO_test<Geom_traits>                            Base;

public:
 typedef typename Base::Point_2                           Point_2;
  typedef typename Base::X_monotone_curve_2               X_monotone_curve_2;
  typedef typename Base::Curve_2                          Curve_2;

  typedef typename Base::Points_vector                    Points_vector;
  typedef typename Base::Xcurves_vector                   Xcurves_vector;
  typedef typename Base::Curves_vector                    Curves_vector;

  typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                          Arrangement;

  typedef typename Arrangement::Vertex_handle             Vertex_handle;
  typedef typename Arrangement::Halfedge_handle           Halfedge_handle;
  typedef typename Arrangement::Face_handle               Face_handle;

  typedef typename Arrangement::Vertex_const_handle       Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle         Face_const_handle;

  typedef typename Arrangement::Edge_const_iterator       Edge_const_iterator;
  typedef typename Arrangement::Vertex_const_iterator     Vertex_const_iterator;

  typedef typename std::pair<CGAL::Object, CGAL::Object>  Object_pair;
  typedef typename std::pair<Vertex_const_handle, Object_pair>
                                                          Vert_decomp_entry;
  typedef typename std::list<Vert_decomp_entry>           Vert_decomp_list;

  typedef CGAL::Arr_point_location_result<Arrangement>    Point_location_result;
  typedef typename Point_location_result::Type            Result_type;

protected:
  /*! The geometry traits. */
  const Geom_traits& m_geom_traits;

  /*! The arrangement. */
  Arrangement m_arr;

  /*! Verbosity */
  size_t m_verbose_level;

  /*! Verify the results.
   */
  template <typename InputIterator>
  bool verify(InputIterator begin, InputIterator end);

  /*! Compare the results.
   */
  bool compare(const Result_type& expected, const CGAL::Object& actual);

  /*! print the results.
   */
  void print(const Vert_decomp_entry& result);

public:
  /*! Constructor from a geometry traits object.
   */
  Vertical_decomposition_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Vertical_decomposition_test() { clear(); }

  /*! Perform the test.
   * \return true upon success and false otherwise.
   */
  virtual bool perform();

  /*! Clear the data structure. */
  virtual void clear();

  /*! Initialize the data structure.
   * \return true upon success and false otherwise.
   */
  virtual bool init();

  /*! Set the verbosity level.
   */
  void set_verbose_level(size_t verbose_level);
};

/*!
 * Constructor from a geometry traits object.
 */
template <typename GeomTraits_T, typename TopolTraits_T>
Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::
Vertical_decomposition_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_verbose_level(0)
{}

//! \brief sets the verbosity level.
template <typename GeomTraits_T, typename TopolTraits_T>
void Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::
set_verbose_level(size_t verbose_level)
{ m_verbose_level = verbose_level; }

/*! Clear the data structures */
template <typename GeomTraits_T, typename TopolTraits_T>
void Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::clear()
{
  m_arr.clear();
  Base::clear();
}

template <typename GeomTraits_T, typename TopolTraits_T>
bool Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::init()
{
  // Initialize the input.
  if (!Base::init()) return false;

  // Insert all into the arrangement
  CGAL::insert(m_arr, this->m_xcurves.begin(), this->m_xcurves.end());
  // insert(*m_arr, m_points.begin(), m_points.end());
  CGAL::insert(m_arr, this->m_curves.begin(), this->m_curves.end());

  // Print the size of the arrangement.
  if (m_verbose_level > 1)
    std::cout << "V = " << m_arr.number_of_vertices()
              << ",  E = " << m_arr.number_of_edges()
              << ",  F = " << m_arr.number_of_faces() << std::endl;

  return true;
}

//! \brief performs the test.
template <typename GeomTraits_T, typename TopolTraits_T>
bool Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::perform()
{
  // Apply vertical decomposition.
  Vert_decomp_list results;
  CGAL::decompose(m_arr, std::back_inserter(results));

  // Verify the results.
  return verify(results.begin(), results.end());
}

//! \brief verifies the results.
template <typename GeomTraits_T, typename TopolTraits_T>
template <typename InputIterator>
bool Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::
verify(InputIterator begin, InputIterator end)
{
  typedef TopolTraits_T TopolTraits;

  InputIterator it;
  if (m_verbose_level > 1) for (it = begin; it != end; ++it) print(*it);

  // Compare the results.
  typename TopolTraits::Default_vertical_ray_shooting_strategy vs(m_arr);
  for (it = begin; it != end; ++it) {
    Vertex_const_handle vh = it->first;
    const Object_pair& res = it->second;
    const CGAL::Object& obj_below_actual = res.first;
    const CGAL::Object& obj_above_actual = res.second;

    Result_type obj_below_expected = vs.ray_shoot_down(vh->point());
    Result_type obj_above_expected = vs.ray_shoot_up(vh->point());

    if (!compare(obj_below_expected, obj_below_actual)) return false;
    if (!compare(obj_above_expected, obj_above_actual)) return false;
  }

  return true;
}

template <typename GeomTraits_T, typename TopolTraits_T>
bool Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::
compare(const Result_type& expected, const CGAL::Object& actual)
{
  // Assign object to a fase.
  const Face_const_handle* fh_expected =
    boost::get<Face_const_handle>(&(expected));
  if (fh_expected) {
    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a face."  << std::endl;
      std::cout << "Actual: a vertex." << std::endl;
      return false;
    }
    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a face."  << std::endl;
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    return true;
  }

  // Assign object to a halfedge.
  const Halfedge_const_handle* hh_expected =
    boost::get<Halfedge_const_handle>(&(expected));
  if (hh_expected) {
    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      if (*hh_expected == hh_actual) return true;

      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
                << std::endl;
      std::cout << "Actual: a different halfedge." << hh_actual->curve()
                << std::endl;
      return false;
    }

    std::cout << "Error: vertical decomposition!" << std::endl;
    std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
              << std::endl;

    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      std::cout << "Actual: a vertex, " << vh_actual->point() << std::endl;
      return false;
    }

    Face_const_handle fh_actual;
    if (CGAL::assign(fh_actual, actual)) {
      std::cout << "Actual: a face." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a vertex.
  const Vertex_const_handle* vh_expected =
    boost::get<Vertex_const_handle>(&(expected));
  if (vh_expected) {
    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      if (*vh_expected == vh_actual) return true;

      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a vertex, " << (*vh_expected)->point()
                << std::endl;
      std::cout << "Actual: a different vertex, " << vh_actual->point()
                << std::endl;
      return false;
    }

    std::cout << "Error: vertical decomposition!" << std::endl;
    std::cout << "Expected: a vertex, " << (*vh_expected)->point() << std::endl;

    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      std::cout << "Actual: a halfedge, " << hh_actual->curve() << std::endl;
      return false;
    }

    Face_const_handle fh_actual;
    if (CGAL::assign(fh_actual, actual)) {
      std::cout << "Actual: a face." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }
  std::cout << "Error: Unknown!" << std::endl;
  return false;
}

//! \brief prints the results.
template <typename GeomTraits_T, typename TopolTraits_T>
void Vertical_decomposition_test<GeomTraits_T, TopolTraits_T>::
print(const Vert_decomp_entry& result)
{
  // Print the result.
  Vertex_const_handle vh;
  Halfedge_const_handle hh;
  Face_const_handle fh;

  const Object_pair& res = result.second;
  std::cout << "Vertex (" << result.first->point() << ") : ";

  std::cout << "  feature below: ";
  if (CGAL::assign(hh, res.first)) std::cout << '[' << hh->curve() << ']';
  else if (CGAL::assign(vh, res.first)) std::cout << '(' << vh->point() << ')';
  else if (CGAL::assign(fh, res.first)) std::cout << "NONE";
  else std::cout << "EMPTY";

  std::cout << "  feature above: ";
  if (CGAL::assign(hh, res.second)) std::cout << '[' << hh->curve() << ']';
  else if (CGAL::assign(vh, res.second)) std::cout << '(' << vh->point() << ')';
  else if (CGAL::assign(fh, res.second)) std::cout << "NONE";
  else std::cout << "EMPTY";

  std::cout << std::endl;
}

#endif
