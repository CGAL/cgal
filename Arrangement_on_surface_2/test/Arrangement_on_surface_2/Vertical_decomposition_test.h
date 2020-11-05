#ifndef CGAL_BATCHED_POINT_LOCATION_TEST_H
#define CGAL_BATCHED_POINT_LOCATION_TEST_H

#include <vector>
#include <list>
#include <utility>


#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Arr_point_location_result.h>

#include "IO_test.h"

/*! Point location test */
template <typename GeomTraits_2, typename TopolTraits>
class Vertical_decomposition_test : public IO_test<GeomTraits_2> {
public:
  typedef GeomTraits_2                                    Geom_traits;
  typedef TopolTraits                                   Topol_traits;

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

  typedef boost::variant<Vertex_const_handle, Halfedge_const_handle,
                         Face_const_handle>               Cell_type;
  typedef boost::optional<Cell_type>                      Vert_type;
  typedef typename std::pair<Vert_type, Vert_type>        Vert_pair;
  typedef typename std::pair<Vertex_const_handle, Vert_pair>
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
  bool compare(const Result_type& expected, Vert_type actual);

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
template <typename GeomTraits_2, typename TopolTraits>
Vertical_decomposition_test<GeomTraits_2, TopolTraits>::
Vertical_decomposition_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_verbose_level(0)
{}

//! \brief sets the verbosity level.
template <typename GeomTraits_2, typename TopolTraits>
void Vertical_decomposition_test<GeomTraits_2, TopolTraits>::
set_verbose_level(size_t verbose_level)
{ m_verbose_level = verbose_level; }

/*! Clear the data structures */
template <typename GeomTraits_2, typename TopolTraits>
void Vertical_decomposition_test<GeomTraits_2, TopolTraits>::clear()
{
  m_arr.clear();
  Base::clear();
}

template <typename GeomTraits_2, typename TopolTraits>
bool Vertical_decomposition_test<GeomTraits_2, TopolTraits>::init()
{
  // Initialize the input.
  if (! Base::init()) return false;

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
template <typename GeomTraits_2, typename TopolTraits>
bool Vertical_decomposition_test<GeomTraits_2, TopolTraits>::perform()
{
  // Apply vertical decomposition.
  Vert_decomp_list results;
  CGAL::decompose(m_arr, std::back_inserter(results));

  // Verify the results.
  return verify(results.begin(), results.end());
}

//! \brief verifies the results.
template <typename GeomTraits_2, typename TopolTraits>
template <typename InputIterator>
bool Vertical_decomposition_test<GeomTraits_2, TopolTraits>::
verify(InputIterator begin, InputIterator end)
{
  if (m_verbose_level > 1) for (auto it = begin; it != end; ++it) print(*it);

  // Compare the results.
  typename TopolTraits::Default_vertical_ray_shooting_strategy vs(m_arr);
  for (auto it = begin; it != end; ++it) {
    Vertex_const_handle vh = it->first;
    const auto& res = it->second;
    auto obj_below_actual = res.first;
    auto obj_above_actual = res.second;

    Result_type obj_below_expected = vs.ray_shoot_down(vh->point());
    Result_type obj_above_expected = vs.ray_shoot_up(vh->point());

    if (! compare(obj_below_expected, obj_below_actual)) return false;
    if (! compare(obj_above_expected, obj_above_actual)) return false;
  }

  return true;
}

template <typename GeomTraits_2, typename TopolTraits>
bool Vertical_decomposition_test<GeomTraits_2, TopolTraits>::
compare(const Result_type& expected, Vert_type actual)
{
  // This test does not test the case where the result is empty!
  if (! actual) return false;

  auto obj = *actual;
  // Assign object to a fase.
  if (const auto* fh_expected = boost::get<Face_const_handle>(&(expected))) {
    if (boost::get<Vertex_const_handle>(&obj)) {
      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a face."  << std::endl;
      std::cout << "Actual: a vertex." << std::endl;
      return false;
    }
    if (boost::get<Halfedge_const_handle>(&obj)) {
      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a face."  << std::endl;
      std::cout << "Actual: a halfedge." << std::endl;
      return false;
    }
    return true;
  }

  // Assign object to a halfedge.
  const auto* hh_expected = boost::get<Halfedge_const_handle>(&(expected));
  if (hh_expected) {
    if (const auto* hh_actual = boost::get<Halfedge_const_handle>(&obj)) {
      if (*hh_expected == *hh_actual) return true;

      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
                << std::endl;
      std::cout << "Actual: a different halfedge." << (*hh_actual)->curve()
                << std::endl;
      return false;
    }

    std::cout << "Error: vertical decomposition!" << std::endl;
    std::cout << "Expected: a halfedge, " << (*hh_expected)->curve()
              << std::endl;

    if (const auto* vh_actual = boost::get<Vertex_const_handle>(&obj)) {
      std::cout << "Actual: a vertex, " << (*vh_actual)->point() << std::endl;
      return false;
    }

    Face_const_handle fh_actual;
    if (boost::get<Face_const_handle>(&obj)) {
      std::cout << "Actual: a face." << std::endl;
      return false;
    }
    std::cout << "Actual: an unknowen object." << std::endl;
    return false;
  }

  // Assign object to a vertex.
  const auto* vh_expected = boost::get<Vertex_const_handle>(&(expected));
  if (vh_expected) {
    if (const auto* vh_actual = boost::get<Vertex_const_handle>(&obj)) {
      if (*vh_expected == *vh_actual) return true;

      std::cout << "Error: vertical decomposition!" << std::endl;
      std::cout << "Expected: a vertex, " << (*vh_expected)->point()
                << std::endl;
      std::cout << "Actual: a different vertex, " << (*vh_actual)->point()
                << std::endl;
      return false;
    }

    std::cout << "Error: vertical decomposition!" << std::endl;
    std::cout << "Expected: a vertex, " << (*vh_expected)->point() << std::endl;

    if (const auto* hh_actual = boost::get<Halfedge_const_handle>(&obj)) {
      std::cout << "Actual: a halfedge, " << (*hh_actual)->curve() << std::endl;
      return false;
    }

    if (boost::get<Face_const_handle>(&obj)) {
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
template <typename GeomTraits_2, typename TopolTraits>
void Vertical_decomposition_test<GeomTraits_2, TopolTraits>::
print(const Vert_decomp_entry& result)
{
  const auto& res = result.second;
  auto obj_below = res.first;
  auto obj_above = res.second;

  std::cout << "Vertex (" << result.first->point() << ") : ";

  assert(obj_below);
  auto obj = *obj_below;
  std::cout << "  feature below: ";
  if (const auto* hh  = boost::get<Halfedge_const_handle>(&obj))
    std::cout << '[' << (*hh)->curve() << ']';
  else if (const auto* vh = boost::get<Vertex_const_handle>(&obj))
    std::cout << '(' << (*vh)->point() << ')';
  else if (const auto* fh = boost::get<Face_const_handle>(&obj))
    std::cout << "NONE";
  else std::cout << "EMPTY";

  assert(obj_above);
  obj = *obj_above;
  std::cout << "  feature above: ";
  if (const auto* hh  = boost::get<Halfedge_const_handle>(&obj))
    std::cout << '[' << (*hh)->curve() << ']';
  else if (const auto* vh = boost::get<Vertex_const_handle>(&obj))
    std::cout << '(' << (*vh)->point() << ')';
  else if (const auto* vh = boost::get<Face_const_handle>(&obj))
    std::cout << "NONE";
  else std::cout << "EMPTY";

  std::cout << std::endl;
}

#endif
