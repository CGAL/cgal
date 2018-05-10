#ifndef CGAL_POINT_LOCATION_TEST_H
#define CGAL_POINT_LOCATION_TEST_H

#include <array>

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_specified_points_generator.h>
#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arr_triangulation_point_location.h>

#include <CGAL/disable_warnings.h>

#include "IO_test.h"

/*! Point location test */
template <typename GeomTraits, typename TopolTraits>
class Point_location_test : public IO_test<GeomTraits> {
public:
  typedef GeomTraits                                    Geom_traits;
  typedef TopolTraits                                   Topol_traits;
  typedef IO_test<Geom_traits>                          Base;

  typedef typename Base::Point_2                        Point_2;
  typedef typename Base::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Base::Curve_2                        Curve_2;

  typedef typename Base::Points_vector                  Points_vector;
  typedef typename Base::Xcurves_vector                 Xcurves_vector;
  typedef typename Base::Curves_vector                  Curves_vector;

  typedef CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>
                                                        Arrangement;

  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Face_handle             Face_handle;

  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;

  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;

  typedef typename Points_vector::iterator              Point_iterator;

  typedef std::vector<CGAL::Object>                     Objects_vector;
  typedef Objects_vector::iterator                      Object_iterator;

  typedef typename boost::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>    Cell_handle;
  typedef std::vector<Cell_handle>                      Variants_vector;
  typedef typename Variants_vector::iterator            Variant_iterator;

protected:
  typedef typename CGAL::Arr_naive_point_location<Arrangement>
                                                    Naive_pl;
  typedef typename CGAL::Arr_simple_point_location<Arrangement>
                                                    Simple_pl;
  typedef typename CGAL::Arr_walk_along_line_point_location<Arrangement>
                                                    Walk_pl;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement>
                                                    Lm_pl;
  typedef typename CGAL::Arr_random_landmarks_generator<Arrangement>
                                                    Random_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Random_lm_generator>
                                                    Lm_random_pl;
  typedef typename CGAL::Arr_grid_landmarks_generator<Arrangement>
                                                    Grid_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Grid_lm_generator>
                                                    Lm_grid_pl;
  typedef typename CGAL::Arr_halton_landmarks_generator<Arrangement>
                                                    Halton_lm_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Halton_lm_generator>
                                                    Lm_halton_pl;
  typedef typename CGAL::Arr_middle_edges_landmarks_generator<Arrangement>
                                                    Middle_edges_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Middle_edges_generator>
    Lm_middle_edges_pl;
  typedef typename CGAL::Arr_landmarks_specified_points_generator<Arrangement>
    Specified_points_generator;
  typedef typename CGAL::Arr_landmarks_point_location<Arrangement,
                                                      Specified_points_generator>
    Lm_specified_points_pl;

  typedef CGAL::Arr_triangulation_point_location<Arrangement>
    Triangulation_pl;

  typedef typename CGAL::Arr_trapezoid_ric_point_location<Arrangement>
    Trapezoid_ric_pl;

  // ===> Add new point location type here <===

protected:
  /*! The geometry traits */
  const Geom_traits& m_geom_traits;

  /*! The arrangement */
  Arrangement* m_arr;

  /*! The input data file of the query points*/
  std::string m_filename_queries;

  /*! The query points */
  Points_vector m_query_points;

  enum Pl_strategy {
    NAIVE_PL = 0,
    SIMPLE_PL,
    WALK_PL,
    LM_PL,
    LM_RANDOM_PL,
    LM_GRID_PL,
    LM_HALTON_PL,
    LM_MIDDLE_EDGES_PL,
    LM_SPECIFIED_POINTS_PL,
    TRIANGULATION_PL,
    TRAPEZOID_RIC_PL,
    TRAPEZOID_RIC_NO_GUARANTEE_PL,
    NUM_PL_STRATEGIES
  };

#define NAIVE_PL_ACTIVE 1
#define SIMPLE_PL_ACTIVE 1
#define WALK_PL_ACTIVE 1
#define LM_PL_ACTIVE 1
#define LM_RANDOM_PL_ACTIVE 1
#define LM_GRID_PL_ACTIVE 1
#define LM_HALTON_PL_ACTIVE 1
#define LM_MIDDLE_EDGES_PL_ACTIVE 1
#define LM_SPECIFIED_POINTS_PL_ACTIVE 1
#define TRIANGULATION_PL_ACTIVE 1
#define TRAPEZOID_RIC_PL_ACTIVE 1
#define TRAPEZOID_RIC_NO_GUARANTEE_PL_ACTIVE 1
#define NUM_PL_STRATEGIE_ACTIVE 1

  struct Locator {
    typedef boost::variant<Naive_pl*,
                           Simple_pl*,
                           Walk_pl*,
#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)
                           Lm_pl*,
                           Lm_random_pl*,
                           Lm_grid_pl*,
                           Lm_halton_pl*,
                           Lm_middle_edges_pl*,
                           Lm_specified_points_pl*,
#endif
#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
                           Triangulation_pl*,
#endif
                           Trapezoid_ric_pl*>       Pl_variant;
    Pl_variant m_variant;

    //! The name of the locator.
    std::string m_name;

    bool m_active;
  };

  std::array<Locator, NUM_PL_STRATEGIES> m_locators;

  // Landmark point-location generators
  Random_lm_generator* m_random_g;
  Grid_lm_generator* m_grid_g;
  Halton_lm_generator* m_halton_g;
  Middle_edges_generator* m_middle_edges_g;
  Specified_points_generator* m_specified_points_g;

  /*! Verify the results (old version).
   */
  int verify(Objects_vector objs[NUM_PL_STRATEGIES],
             size_t size, size_t pls_num);

  /*! Verify the results (new version).
   */
  int verify(Variants_vector objs[NUM_PL_STRATEGIES],
             size_t size, size_t pls_num);

public:
  /*! Constructor */
  Point_location_test(const Geom_traits& geom_traits);

  /*! Destructor */
  virtual ~Point_location_test()
  {
    deallocate_arrangement();
    clear();
    deallocate_pl_strategies();
  }

  void set_filenames(const char* points_filename, const char* xcurves_filename,
                     const char* curves_filename, const char* queries_filename);

  //! Initialize the data structures.
  virtual bool init();

  //! Perform the test.
  virtual bool perform();

  //! Clear the data structures.
  virtual void clear();

  bool allocate_arrangement();

  void deallocate_arrangement();

  bool construct_arrangement();

  void clear_arrangement();

  bool allocate_pl_strategies();

  bool construct_pl_strategies();

  bool attach_pl_strategies();

  void deallocate_pl_strategies();

  template <typename Point_location,
            typename InputIterator, typename OutputIterator>
  void query(Point_location& pl, const char* type,
             InputIterator begin, InputIterator end, OutputIterator oi)
  {
    CGAL::Timer timer;
    timer.reset(); timer.start();
    for (auto piter = begin; piter != end; ++piter) *oi++ = pl.locate(*piter);
    timer.stop();
    std::cout << type << " location took " << timer.time() << std::endl;
  }

private:
  //! Initialize point location.
  template <typename Strategy, Pl_strategy id>
  void init_pl(const std::string& name)
  {
    m_locators[id].m_variant = static_cast<Strategy*>(nullptr);
    m_locators[id].m_name = name;
    m_locators[id].m_active = true;
  }

  //! Allocate point location.
  template <typename Strategy, Pl_strategy id>
  bool allocate_pl()
  {
    auto* locator = new Strategy();
    if (! locator) return false;
    m_locators[id].m_variant = locator;
    return true;
  }

  //! Construct point location.
  template <typename Strategy, Pl_strategy id>
  bool construct_pl()
  {
    auto* locator = new Strategy(*m_arr);
    if (! locator) return false;
    m_locators[id].m_variant = locator;
    return true;
  }

  //! Construct landmark point-location with a generator.
  template <typename Strategy, Pl_strategy id, typename Generator>
  bool construct_pl(Generator* generator)
  {
    auto* locator = new Strategy(*m_arr, generator);
    if (! locator) return false;
    m_locators[id].m_variant = locator;
    return true;
  }

  //! Delete point location.
  template <typename Strategy, Pl_strategy id>
  void deallocate_pl()
  {
    auto* strategy = boost::get<Strategy*>(m_locators[id].m_variant);
    if (strategy) {
      delete strategy;
      m_locators[id].m_variant = static_cast<Strategy*>(nullptr);
    }
  }

  //! Attach point location.
  template <typename Strategy, Pl_strategy id>
  void attach_pl()
  {
    if (! m_locators[id].m_active) return;
    auto* strategy = boost::get<Strategy*>(m_locators[id].m_variant);
    strategy->attach(*m_arr);
  }

  //! Attach landmark point location with a generator.
  template <typename Strategy, Pl_strategy id, typename Generator>
  void attach_pl(Generator* generator)
  {
    if (! m_locators[id].m_active) return;
    auto* strategy = boost::get<Strategy*>(m_locators[id].m_variant);
    strategy->attach(*m_arr, generator);
  }

  //! Query using point location.
  template <typename Strategy, Pl_strategy id, typename T>
  void query_pl(T& objs)
  {
    if (! m_locators[id].m_active) return;
    const auto& name = m_locators[id].m_name;
    auto* strategy = boost::get<Strategy*>(m_locators[id].m_variant);
    query(*strategy, name.c_str(), m_query_points.begin(), m_query_points.end(),
          std::back_inserter(objs));
  }

  //! Measure the time consumption of an operation
  template <Pl_strategy id, typename Timer, typename UnaryOperation>
  void measure(Timer& timer, UnaryOperation op)
  {
    if (! m_locators[id].m_active) return;
    const auto& name = m_locators[id].m_name;
    timer.reset();
    timer.start();
    op();
    timer.stop();
    std::cout << name.c_str() << " took " << timer.time() << std::endl;
  }
};

//! Constructor.
template <typename GeomTraits, typename TopolTraits>
Point_location_test<GeomTraits, TopolTraits>::
Point_location_test(const Geom_traits& geom_traits) :
  Base(geom_traits),
  m_geom_traits(geom_traits),
  m_arr(nullptr),
  m_random_g(nullptr),
  m_grid_g(nullptr),
  m_halton_g(nullptr),
  m_middle_edges_g(nullptr),
  m_specified_points_g(nullptr)
{
  init_pl<Naive_pl, NAIVE_PL>("Naive");
  init_pl<Simple_pl, SIMPLE_PL>("Simple");
  init_pl<Walk_pl, WALK_PL>("Walk");

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  init_pl<Lm_pl, LM_PL>("Landmarks (vertices)");
  init_pl<Lm_random_pl, LM_RANDOM_PL>("Landmarks random");
  init_pl<Lm_grid_pl, LM_GRID_PL>("Landmarks grid");
  init_pl<Lm_halton_pl, LM_HALTON_PL>("Landmarks Halton");
  init_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>("Landmarks middle edges");
  init_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>("Landmarks specified points");
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  init_pl<Triangulation_pl, TRIANGULATION_PL>("Triangulation");
#endif

  init_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>("Trapezoidal RIC");
  init_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>("Trapezoidal RIC without guarantees");
}

//! Set the file names.
template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::
set_filenames(const char* points_filename,
              const char* xcurves_filename,
              const char* curves_filename,
              const char* queries_filename)
{
  Base::set_filenames(points_filename, xcurves_filename, curves_filename);
  m_filename_queries.assign(queries_filename);
}

//! Initialize the data structures.
template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::init()
{
  if (!Base::init()) return false;

  // Read the query points
  if (!this->read_points(m_filename_queries.c_str(), m_query_points))
    return false;

  return true;
}

//! Clear the data structures.
template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::clear()
{
  Base::clear();
  m_query_points.clear();
  m_filename_queries.clear();
}

//! Clear the data structures.
template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::deallocate_pl_strategies()
{
  deallocate_pl<Naive_pl, NAIVE_PL>();
  deallocate_pl<Simple_pl, SIMPLE_PL>();
  deallocate_pl<Walk_pl, WALK_PL>();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  deallocate_pl<Lm_pl, LM_PL>();
  deallocate_pl<Lm_random_pl, LM_RANDOM_PL>();
  deallocate_pl<Lm_grid_pl, LM_GRID_PL>();
  deallocate_pl<Lm_halton_pl, LM_HALTON_PL>();
  deallocate_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>();
  deallocate_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>();

  // Free Generators
  if (m_random_g) {
    delete m_random_g;
    m_random_g = nullptr;
  }
  if (m_grid_g) {
    delete m_grid_g;
    m_grid_g = nullptr;
  }
  if (m_halton_g) {
    delete m_halton_g;
    m_halton_g = nullptr;
  }
  if (m_middle_edges_g) {
    delete m_middle_edges_g;
    m_middle_edges_g = nullptr;
  }
  if (m_specified_points_g) {
    delete m_specified_points_g;
    m_specified_points_g = nullptr;
  }
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  deallocate_pl<Triangulation_pl, TRIANGULATION_PL>();
#endif

  deallocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>();
  deallocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>();
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::allocate_arrangement()
{
  if (!(m_arr = new Arrangement(&m_geom_traits))) return false;
  return true;
}

template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::
deallocate_arrangement()
{
  if (m_arr) {
    delete m_arr;
    m_arr = nullptr;
  }
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::
construct_arrangement()
{
  // Insert all into the arrangement
  CGAL::insert(*m_arr, this->m_xcurves.begin(), this->m_xcurves.end());
  // insert(*m_arr, m_points.begin(), m_points.end());
  CGAL::insert(*m_arr, this->m_curves.begin(), this->m_curves.end());

  // Print the size of the arrangement.
  std::cout << "V = " << m_arr->number_of_vertices()
            << ",  E = " << m_arr->number_of_edges()
            << ",  F = " << m_arr->number_of_faces() << std::endl;

  return true;
}

template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::clear_arrangement()
{
  if (m_arr) m_arr->clear();
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::allocate_pl_strategies()
{
  // Allocate all point location strategies.
  if (! allocate_pl<Naive_pl, NAIVE_PL>()) return false;
  if (! allocate_pl<Simple_pl, SIMPLE_PL>()) return false;
  if (! allocate_pl<Walk_pl, WALK_PL>()) return false;

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  if (! allocate_pl<Lm_pl, LM_PL>()) return false;
  if (! allocate_pl<Lm_random_pl, LM_RANDOM_PL>()) return false;
  if (! allocate_pl<Lm_grid_pl, LM_GRID_PL>()) return false;
  if (! allocate_pl<Lm_halton_pl, LM_HALTON_PL>()) return false;
  if (! allocate_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>()) return false;
  if (! allocate_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>())
    return false;

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  if (! allocate_pl<Triangulation_pl, TRIANGULATION_PL>()) return false;
#endif

  if (! allocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>()) return false;
  if (! allocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>())
    return false;

  return true;
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::construct_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  construct_pl<Naive_pl, NAIVE_PL>();
  construct_pl<Simple_pl, SIMPLE_PL>();
  construct_pl<Walk_pl, WALK_PL>();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  measure<LM_PL>(timer, [&](){ construct_pl<Lm_pl, LM_PL>(); });
  // auto op = std::bind(&Point_location_test::construct_pl<Lm_pl, LM_PL>, this);
  // measure<LM_PL>(timer, op);

  measure<LM_RANDOM_PL>(timer, [&](){
      m_random_g = new Random_lm_generator(*m_arr);
      construct_pl<Lm_random_pl, LM_RANDOM_PL>(m_random_g);
    });

  measure<LM_GRID_PL>(timer, [&](){
      m_grid_g = new Grid_lm_generator(*m_arr);
      construct_pl<Lm_grid_pl, LM_GRID_PL>(m_grid_g);
    });

  measure<LM_HALTON_PL>(timer, [&](){
      m_halton_g = new Halton_lm_generator(*m_arr);
      construct_pl<Lm_halton_pl, LM_HALTON_PL>(m_halton_g);
    });

  measure<LM_MIDDLE_EDGES_PL>(timer, [&](){
      m_middle_edges_g = new Middle_edges_generator(*m_arr);
      construct_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(m_middle_edges_g);
    });

  measure<LM_SPECIFIED_POINTS_PL>(timer, [&](){
      m_specified_points_g = new Specified_points_generator(*m_arr);
      construct_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>
        (m_specified_points_g);
    });

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  measure<TRIANGULATION_PL>(timer, [&](){
      construct_pl<Triangulation_pl, TRIANGULATION_PL>();
    });
#endif

  measure<TRAPEZOID_RIC_PL>(timer, [&](){
      construct_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>();
    });

  measure<TRAPEZOID_RIC_NO_GUARANTEE_PL>(timer, [&](){
      construct_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>();
    });

  std::cout << std::endl;

  return true;
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::attach_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  attach_pl<Naive_pl, NAIVE_PL>();
  attach_pl<Simple_pl, SIMPLE_PL>();
  attach_pl<Walk_pl, WALK_PL>();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  measure<LM_PL>(timer, [&](){ attach_pl<Lm_pl, LM_PL>(); });

  measure<LM_RANDOM_PL>(timer, [&](){
      m_random_g = new Random_lm_generator(*m_arr);
      attach_pl<Lm_random_pl, LM_RANDOM_PL>(m_random_g);
    });

  measure<LM_GRID_PL>(timer, [&](){
      m_grid_g = new Grid_lm_generator(*m_arr);
      attach_pl<Lm_grid_pl, LM_GRID_PL>(m_grid_g);
    });

  measure<LM_HALTON_PL>(timer, [&](){
      m_halton_g = new Halton_lm_generator(*m_arr);
      attach_pl<Lm_halton_pl, LM_HALTON_PL>(m_halton_g);
    });

  measure<LM_MIDDLE_EDGES_PL>(timer, [&](){
      m_middle_edges_g = new Middle_edges_generator(*m_arr);
      attach_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(m_middle_edges_g);
    });

  measure<LM_SPECIFIED_POINTS_PL>(timer, [&](){
      m_specified_points_g = new Specified_points_generator(*m_arr);
      attach_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>
        (m_specified_points_g);
    });
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  measure<TRIANGULATION_PL>(timer, [&](){
      attach_pl<Triangulation_pl, TRIANGULATION_PL>();
    });
#endif

  measure<TRAPEZOID_RIC_PL>(timer, [&](){
      attach_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>();
    });

  measure<TRAPEZOID_RIC_NO_GUARANTEE_PL>(timer, [&](){
      auto& var = m_locators[TRAPEZOID_RIC_NO_GUARANTEE_PL].m_variant;
      auto* strategy = boost::get<Trapezoid_ric_pl*>(var);
      strategy->with_guarantees(false);
      attach_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>();
    });

  std::cout << std::endl;

  return true;
}

// Perform the test
template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::perform()
{
#if ((CGAL_ARR_POINT_LOCATION_VERSION < 2) || \
     defined(CGAL_ARR_POINT_LOCATION_CONVERSION))
  Objects_vector objs[NUM_PL_STRATEGIES];
#else
  Variants_vector objs[NUM_PL_STRATEGIES];
#endif

  // Locate the points in the list using all point location strategies.

  // std::cout << "Time in seconds" << std::endl;
  std::cout << std::endl;

  size_t pl_index = 0;

  query_pl<Naive_pl, NAIVE_PL>(objs[pl_index++]);
  query_pl<Simple_pl, SIMPLE_PL>(objs[pl_index++]);
  query_pl<Walk_pl, WALK_PL>(objs[pl_index++]);

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  query_pl<Lm_pl, LM_PL>(objs[pl_index++]);
  query_pl<Lm_random_pl, LM_RANDOM_PL>(objs[pl_index++]);
  query_pl<Lm_grid_pl, LM_GRID_PL>(objs[pl_index++]);
  query_pl<Lm_halton_pl, LM_HALTON_PL>(objs[pl_index++]);
  query_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(objs[pl_index++]);
  query_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>(objs[pl_index++]);

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  query_pl<Triangulation_pl, TRIANGULATION_PL>(objs[pl_index++]);
#endif

  query_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>(objs[pl_index++]);
  query_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>(objs[pl_index++]);

  std::cout << std::endl;

  // Number of point location strategies used.
  size_t pls_num = pl_index;
  std::cout << "Number of strategies is " << pls_num << std::endl;

  // End Location

  int result = 0;

  // get size of objects
  size_t size = objs[0].size();
  std::cout << "size is " << size << std::endl;

  for (pl_index = 0; pl_index < pls_num; ++pl_index) {
    if (size != objs[pl_index].size()) {
      std::cout << "Error: size of pl number " << pl_index << " is "
                << objs[pl_index].size() << std::endl;
      result += -1;
      break;
    }
  }

  result += verify(objs, size, pls_num);
  return (result == 0);
}

//! Verify the results
template <typename GeomTraits, typename TopolTraits>
int Point_location_test<GeomTraits, TopolTraits>::
verify(Objects_vector objs[NUM_PL_STRATEGIES], size_t size, size_t pls_num)
{
  Vertex_const_handle vh_ref, vh_cur;
  Halfedge_const_handle hh_ref, hh_cur;
  Face_const_handle fh_ref, fh_cur;

  int result = 0;

  // Assign and check results
  size_t qi; //qi is the query point index
  for (qi = 0; qi < size; ++qi) {
    // Assign object to a face
    if (CGAL::assign(fh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          if (fh_cur != fh_ref) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expecte: a face." << std::endl;
            std::cout << "Actual: a different face" << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expecte: a face." << std::endl;
        result += -1;
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      //if (fh_ref->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
      continue;
    }

    // Assign object to a halfedge
    if (CGAL::assign (hh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          if ((hh_cur != hh_ref) && (hh_cur->twin() != hh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a halfedge, " << hh_ref->curve()
                      << std::endl;
            std::cout << "Actual: a different halfedge: " << hh_cur->curve()
                      << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a halfedge, " << hh_ref->curve() << std::endl;
        result += -1;
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      continue;
    }

    // Assign object to a vertex
    if (CGAL::assign(vh_ref, objs[0][qi])) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        if (CGAL::assign(vh_cur, objs[pl][qi])) {
          if (vh_cur != vh_ref) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a vertex, "<< vh_ref->point() << std::endl;
            std::cout << "Actual: a different vertex, "<< vh_cur->point()
                      << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a vertex, "<< vh_ref->point() << std::endl;
        result += -1;
        if (CGAL::assign(fh_cur, objs[pl][qi])) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        if (CGAL::assign(hh_cur, objs[pl][qi])) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknown object." << std::endl;
      }
      continue;
    }

    std::cout << "Error: Unknown!" << std::endl;
    result += -1;
  }
  return result;
}

// Verify the results
template <typename GeomTraits, typename TopolTraits>
int Point_location_test<GeomTraits, TopolTraits>::
verify(Variants_vector objs[NUM_PL_STRATEGIES], size_t size, size_t pls_num)
{
  int result = 0;

  // Assign and check results
  for (size_t qi = 0; qi < size; ++qi) {
    // Assign object to a face
    auto* fh_ref = boost::get<Face_const_handle>(&(objs[0][qi]));
    if (fh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	auto* fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          if ((*fh_cur) != (*fh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a face." << std::endl;
            std::cout << "Actual: a different face." << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a face." << std::endl;
        result += -1;
	auto* hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
	if (hh_cur) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
	auto* vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      //if ((*fh_ref)->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
      continue;
    }

    // Assign object to a halfedge
    auto* hh_ref = boost::get<Halfedge_const_handle>(&(objs[0][qi]));
    if (hh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	auto* hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
        if (hh_cur) {
          if (((*hh_cur) != (*hh_ref)) && ((*hh_cur)->twin() != (*hh_ref))) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a halfedge, " << (*hh_ref)->curve()
                      << std::endl;
            std::cout << "Actual: a different halfedge, " << (*hh_cur)->curve()
                      << std::endl;
            result += -1;
          }
          continue;
        }
        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a halfedge, " << (*hh_ref)->curve()
                  << std::endl;
        result += -1;
	auto* fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
	auto* vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      continue;
    }

    // Assign object to a vertex
    auto* vh_ref = boost::get<Vertex_const_handle>(&(objs[0][qi]));
    if (vh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
	auto vh_cur = boost::get<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          if ((*vh_cur) != (*vh_ref)) {
            std::cout << "Error: point location number " << pl << std::endl;
            std::cout << "Expected: a vertex: "<< (*vh_ref)->point()
                      << std::endl;
            std::cout << "Actual: a different vertex: "<< (*vh_cur)->point()
                      << std::endl;
            result += -1;
          }
          continue;
        }
        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a vertex: "<< (*vh_ref)->point() << std::endl;
        result += -1;
	auto fh_cur = boost::get<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
	auto hh_cur = boost::get<Halfedge_const_handle>(&(objs[pl][qi]));
        if (hh_cur) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknown object." << std::endl;
      }
      continue;
    }

    std::cout << "Error: Unknown!" << std::endl;
    result += -1;
  }
  return result;
}

#include <CGAL/enable_warnings.h>

#endif
