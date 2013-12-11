#ifndef CGAL_BATCHED_POINT_LOCATION_TEST_H
#define CGAL_BATCHED_POINT_LOCATION_TEST_H

#include <vector>
#include <list>
#include <utility>

#include <CGAL/basic.h>
#include <CGAL/Arr_batched_point_location.h>
#include <CGAL/Arr_point_location_result.h>

#include "IO_test.h"

/*! Point location test */
template <typename GeomTraits_T, typename TopolTraits_T>
class Batched_point_location_test : public IO_test<GeomTraits_T> {
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

  typedef typename Points_vector::iterator                Point_iterator;

  typedef typename boost::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>      Cell_handle;

  typedef CGAL::Arr_point_location_result<Arrangement>    Point_location_result;
  typedef typename Point_location_result::Type            Result_type;
  typedef std::pair<Point_2, Result_type>                 Query_result;

protected:
  /*! The geometry traits. */
  const Geom_traits& m_geom_traits;

  /*! The arrangement. */
  Arrangement m_arr;

  /*! The input data file of the query points. */
  std::string m_filename_queries;

  /*! The query points. */
  Points_vector m_query_points;

  /*! Verbosity */
  size_t m_verbose_level;

  /*! Verify the results.
   */
  template <typename InputIterator>
  bool verify(InputIterator begin, InputIterator end);

  /*! print the results.
   */
  void print(const std::pair<Point_2, Cell_handle>& res);

  /*! print the results.
   */
  void print(const std::pair<Point_2, CGAL::Object>& res);

  /*! Compare the results.
   */
  bool compare(const Cell_handle& expected, const Cell_handle& actual);

  /*! Compare the results.
   */
  bool compare(const CGAL::Object& expected, const CGAL::Object& actual);

public:
  /*! Constructor from a geometry traits object.
   */
  Batched_point_location_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Batched_point_location_test() { clear(); }

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

  /*! Set the file names.
   */
  void set_filenames(const char* points_filename, const char* xcurves_filename,
                     const char* curves_filename, const char* queries_filename);

  /*! Set the verbosity level.
   */
  void set_verbose_level(size_t verbose_level);
};

/*!
 * Constructor from a geometry traits object.
 */
template <typename GeomTraits_T, typename TopolTraits_T>
Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
Batched_point_location_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_verbose_level(0)
{}

//! \brief sets the file names.
template <typename GeomTraits_T, typename TopolTraits_T>
void Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
set_filenames(const char* points_filename,
              const char* xcurves_filename,
              const char* curves_filename,
              const char* queries_filename)
{
  Base::set_filenames(points_filename, xcurves_filename, curves_filename);
  m_filename_queries.assign(queries_filename);
}

//! \brief sets the verbosity level.
template <typename GeomTraits_T, typename TopolTraits_T>
void Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
set_verbose_level(size_t verbose_level)
{ m_verbose_level = verbose_level; }

/*! Clear the data structures */
template <typename GeomTraits_T, typename TopolTraits_T>
void Batched_point_location_test<GeomTraits_T, TopolTraits_T>::clear()
{
  m_arr.clear();
  Base::clear();
  m_query_points.clear();
  m_filename_queries.clear();
}

template <typename GeomTraits_T, typename TopolTraits_T>
bool Batched_point_location_test<GeomTraits_T, TopolTraits_T>::init()
{
  // Initialize the input.
  if (!Base::init()) return false;

  // Read the query points
  if (!this->read_points(m_filename_queries.c_str(), m_query_points))
    return false;

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
bool Batched_point_location_test<GeomTraits_T, TopolTraits_T>::perform()
{
  // Apply batched point location.
  std::list<Query_result>  results;
  locate(m_arr, m_query_points.begin(), m_query_points.end(),
         std::back_inserter(results));

  // Verify the results.
  return verify(results.begin(), results.end());
}

//! \brief verifies the results.
template <typename GeomTraits_T, typename TopolTraits_T>
template <typename InputIterator>
bool Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
verify(InputIterator begin, InputIterator end)
{
  typedef TopolTraits_T TopolTraits;

  typename TopolTraits::Default_point_location_strategy pl(m_arr);
  for (InputIterator it = begin; it != end; ++it) {

    if (m_verbose_level > 1)
      print(*it);

    // Perform (single) point location.
    Result_type obj = pl.locate(it->first);

    // Compare the results.
    if (!compare(obj, it->second)) return false;
  }
  return true;
}

//! \brief compares the results.
template <typename GeomTraits_T, typename TopolTraits_T>
bool Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
compare(const Cell_handle& expected, const Cell_handle& actual)
{
  // Assign object to a face
  const Face_const_handle* fh_expected =
    boost::get<Face_const_handle>(&(expected));
  if (fh_expected) {
    const Face_const_handle* fh_actual =
      boost::get<Face_const_handle>(&(actual));
    if (fh_actual) {
      if ((*fh_actual) == (*fh_expected)) return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a face." << std::endl;
      std::cout << "Actual: a different face." << std::endl;
      return false;
    }
    std::cout << "Error: batched point location!" << std::endl;
    std::cout << "Expected a face."  << std::endl;
    const Halfedge_const_handle* hh_actual =
      boost::get<Halfedge_const_handle>(&(actual));
    if (hh_actual) {
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    const Vertex_const_handle* vh_actual =
      boost::get<Vertex_const_handle>(&(actual));
    if (vh_actual) {
      std::cout << "Actual: a vertex." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a halfedge
  const Halfedge_const_handle* hh_expected =
    boost::get<Halfedge_const_handle>(&(expected));
  if (hh_expected) {
    const Halfedge_const_handle* hh_actual =
      boost::get<Halfedge_const_handle>(&(actual));
    if (hh_actual) {
      if (((*hh_actual) == (*hh_expected)) ||
          ((*hh_actual)->twin() == (*hh_expected)))
        return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
                << std::endl;
      std::cout << "Actual: a different halfedge, " << (*hh_actual)->curve()
                << std::endl;
      return false;
    }

    std::cout << "Error: batched point location!" << std::endl;
    std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
              << std::endl;
    const Face_const_handle* fh_actual =
      boost::get<Face_const_handle>(&(actual));
    if (fh_actual) {
      std::cout << "Actual: a face." << std::endl;
      return false;
    }
    const Vertex_const_handle* vh_actual =
      boost::get<Vertex_const_handle>(&(actual));
    if (vh_actual) {
      std::cout << "Actual: a vertex." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a vertex
  const Vertex_const_handle* vh_expected =
    boost::get<Vertex_const_handle>(&(expected));
  if (vh_expected) {
    const Vertex_const_handle* vh_actual =
      boost::get<Vertex_const_handle>(&(actual));
    if (vh_actual) {
      if ((*vh_actual) == (*vh_expected)) return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a vertex, "<< (*vh_expected)->point()
                << std::endl;
      std::cout << "Actual: a different vertex, " << (*vh_actual)->point()
                << std::endl;
      return false;
    }
    std::cout << "Error: batched point location!";
    std::cout << "Expected: a vertex, "<< (*vh_expected)->point() << std::endl;
    const Face_const_handle* fh_actual =
      boost::get<Face_const_handle>(&(actual));
    if (fh_actual) {
      std::cout << "Actual: a face" << std::endl;
      return false;
    }
    const Halfedge_const_handle* hh_actual =
      boost::get<Halfedge_const_handle>(&(actual));
    if (hh_actual) {
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknown object." << std::endl;
    return false;
  }

  std::cout << "Error: Unknown!" << std::endl;
  return false;
}

//! \brief compares the results.
template <typename GeomTraits_T, typename TopolTraits_T>
bool Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
compare(const CGAL::Object& expected, const CGAL::Object& actual)
{
  Face_const_handle fh_expected;
  if (CGAL::assign(fh_expected, expected)) {
    Face_const_handle fh_actual;
    if (CGAL::assign(fh_actual, actual)) {
      if (fh_actual == fh_expected) return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a face" << std::endl;
      std::cout << "Actual: a different face." << std::endl;
      return false;
    }
    std::cout << "Error: batched point location!" << std::endl;
    std::cout << "Expected: a face" << std::endl;
    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      std::cout << "Actual: a vertex." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a halfedge
  Halfedge_const_handle hh_expected;
  if (CGAL::assign(hh_expected, expected)) {
    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      if ((hh_actual == hh_expected) || (hh_actual->twin() == hh_expected))
        return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a halfedge, " << hh_expected->curve()
                << std::endl;
      std::cout << "Actual: a different halfedge, " << hh_actual->curve()
                << std::endl;
      return false;
    }
    std::cout << "Error: batched point location!" << std::endl;
    std::cout << "Expected: a halfedge, " << hh_expected->curve() << std::endl;
    Face_const_handle fh_actual;
    if (CGAL::assign(fh_actual, actual)) {
      std::cout << "Actual: a face" << std::endl;
      return false;;
    }
    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      std::cout << "Actual: a vertex." << std::endl;
      return false;;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a vertex
  Vertex_const_handle vh_expected;
  if (CGAL::assign(vh_expected, expected)) {
    Vertex_const_handle vh_actual;
    if (CGAL::assign(vh_actual, actual)) {
      if (vh_actual == vh_expected) return true;

      std::cout << "Error: batched point location!" << std::endl;
      std::cout << "Expected: a vertex, "<< vh_expected->point() << std::endl;
      std::cout << "Actual: a different vertex: "<< vh_actual->point()
                << std::endl;
      return false;
    }
    std::cout << "Error: batched point location!" << std::endl;
    std::cout << "Expected: a vertex, "<< vh_expected->point() << std::endl;
    Face_const_handle fh_actual;
    if (CGAL::assign(fh_actual, actual)) {
      std::cout << "Actual: a face." << std::endl;
      return false;
    }
    Halfedge_const_handle hh_actual;
    if (CGAL::assign(hh_actual, actual)) {
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknown object." << std::endl;
    return false;
  }

  std::cout << "Error: Unknown!" << std::endl;
  return false;
}

//! \brief prints the results.
template <typename GeomTraits_T, typename TopolTraits_T>
void Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
print(const std::pair<Point_2, Cell_handle>& res)
{
  // Print the results.
  std::cout << "The point (" << res.first << ") is located ";
  if (const Face_const_handle* f =
      boost::get<Face_const_handle>(&(res.second)))       // inside a face
    std::cout << "inside "
              << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
              << " face." << std::endl;
  else if (const Halfedge_const_handle* e =
           boost::get<Halfedge_const_handle>(&(res.second))) // on an edge
    std::cout << "on an edge: " << (*e)->curve() << std::endl;
  else if (const Vertex_const_handle* v =
           boost::get<Vertex_const_handle>(&(res.second)))  // on a vertex
    std::cout << "on "
              << (((*v)->is_isolated()) ? "an isolated" : "a")
              << " vertex: " << (*v)->point() << std::endl;
}

//! \brief prints the results.
template <typename GeomTraits_T, typename TopolTraits_T>
void Batched_point_location_test<GeomTraits_T, TopolTraits_T>::
print(const std::pair<Point_2, CGAL::Object>& res)
{
  // Print the results.
  std::cout << "The point (" << res.first << ") is located ";
  Face_const_handle fh;
  if (CGAL::assign(fh, res.second))
    std::cout << "inside "
              << ((fh->is_unbounded()) ? "the unbounded" : "a bounded")
              << " face." << std::endl;
}

#endif
