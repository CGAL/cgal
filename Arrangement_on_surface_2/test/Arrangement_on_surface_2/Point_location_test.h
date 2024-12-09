#ifndef CGAL_POINT_LOCATION_TEST_H
#define CGAL_POINT_LOCATION_TEST_H

#define NAIVE_PL_ENABLED 1
#define SIMPLE_PL_ENABLED 1
#define WALK_PL_ENABLED 1
#define LM_PL_ENABLED 1
#define LM_RANDOM_PL_ENABLED 1
#define LM_GRID_PL_ENABLED 1
#define LM_HALTON_PL_ENABLED 1
#define LM_MIDDLE_EDGES_PL_ENABLED 1
#define LM_SPECIFIED_POINTS_PL_ENABLED 1
#define TRIANGULATION_PL_ENABLED 1
#define TRAPEZOID_RIC_PL_ENABLED 1
#define TRAPEZOID_RIC_NO_GUARANTEE_PL_ENABLED 1

#include <vector>

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

  typedef typename std::variant<Vertex_const_handle,
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

  typedef std::variant<Naive_pl*,
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

  struct Locator {
    Pl_variant m_variant;

    //! The name of the locator.
    std::string m_name;

    bool m_active;
  };

  std::vector<Locator> m_locators;

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
    typedef InputIterator               Input_iterator;
    CGAL::Timer timer;
    timer.reset(); timer.start();
    for (Input_iterator piter = begin; piter != end; ++piter)
      *oi++ = pl.locate(*piter);
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
  void allocate_pl()
  {
    Strategy* locator = new Strategy();
    assert(locator);
    m_locators[id].m_variant = locator;
  }

  //! Construct point location.
  template <typename Strategy, Pl_strategy id>
  bool construct_pl()
  {
    Strategy* locator = new Strategy(*m_arr);
    if (! locator) return false;
    m_locators[id].m_variant = locator;
    return true;
  }

  //! Construct landmark point-location with a generator.
  template <typename Strategy, Pl_strategy id, typename Generator>
  bool construct_pl(Generator* generator)
  {
    Strategy* locator = new Strategy(*m_arr, generator);
    if (! locator) return false;
    m_locators[id].m_variant = locator;
    return true;
  }

  //! Delete point location.
  template <typename Strategy, Pl_strategy id>
  void deallocate_pl()
  {
    Strategy* strategy = std::get<Strategy*>(m_locators[id].m_variant);
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
    Strategy* strategy = std::get<Strategy*>(m_locators[id].m_variant);
    strategy->attach(*m_arr);
  }

  //! Attach landmark point location with a generator.
  template <typename Strategy, Pl_strategy id, typename Generator>
  void attach_pl(Generator* generator)
  {
    if (! m_locators[id].m_active) return;
    Strategy* strategy = std::get<Strategy*>(m_locators[id].m_variant);
    strategy->attach(*m_arr, generator);
  }

  //! Query using point location.
  template <typename Strategy, Pl_strategy id, typename T>
  void query_pl(T& objs)
  {
    if (! m_locators[id].m_active) return;
    const std::string& name = m_locators[id].m_name;
    Strategy* strategy = std::get<Strategy*>(m_locators[id].m_variant);
    query(*strategy, name.c_str(), m_query_points.begin(), m_query_points.end(),
          std::back_inserter(objs));
  }

  //! Measure the time consumption of an operation
  template <Pl_strategy id, typename Timer, typename UnaryOperation>
  void measure(Timer& timer, UnaryOperation op)
  {
    if (! m_locators[id].m_active) return;
    const std::string& name = m_locators[id].m_name;
    timer.reset();
    timer.start();
    op();
    timer.stop();
    std::cout << name.c_str() << " took " << timer.time() << std::endl;
  }
};

// Naive
#if (NAIVE_PL_ENABLED)
#define INIT_PL_NATIVE()        init_pl<Naive_pl, NAIVE_PL>("Naive")
#define ALLOCATE_PL_NATIVE()    allocate_pl<Naive_pl, NAIVE_PL>()
#define CONSTRUCT_PL_NATIVE()   construct_pl<Naive_pl, NAIVE_PL>()
#define DEALLOCATE_PL_NATIVE()  deallocate_pl<Naive_pl, NAIVE_PL>()
#define ATTACH_PL_NATIVE()      attach_pl<Naive_pl, NAIVE_PL>()
#define QUERY_PL_NATIVE(obj)    query_pl<Naive_pl, NAIVE_PL>(obj)
#else
#define INIT_PL_NATIVE()
#define ALLOCATE_PL_NATIVE()
#define CONSTRUCT_PL_NATIVE()
#define DEALLOCATE_PL_NATIVE()
#define ATTACH_PL_NATIVE()
#define QUERY_PL_NATIVE(obj)
#endif

// Simple
#if (SIMPLE_PL_ENABLED)
#define INIT_PL_SIMPLE()        init_pl<Simple_pl, SIMPLE_PL>("Simple")
#define ALLOCATE_PL_SIMPLE()    allocate_pl<Simple_pl, SIMPLE_PL>()
#define CONSTRUCT_PL_SIMPLE()   construct_pl<Simple_pl, SIMPLE_PL>()
#define DEALLOCATE_PL_SIMPLE()  deallocate_pl<Simple_pl, SIMPLE_PL>()
#define ATTACH_PL_SIMPLE()      attach_pl<Simple_pl, SIMPLE_PL>()
#define QUERY_PL_SIMPLE(obj)    query_pl<Simple_pl, SIMPLE_PL>(obj)
#else
#define INIT_PL_SIMPLE()
#define ALLOCATE_PL_SIMPLE()
#define CONSTRUCT_PL_SIMPLE()
#define DEALLOCATE_PL_SIMPLE()
#define ATTACH_PL_SIMPLE()
#define QUERY_PL_SIMPLE(obj)
#endif

// Walk
#if (WALK_PL_ENABLED)
#define INIT_PL_WALK()          init_pl<Walk_pl, WALK_PL>("Walk")
#define ALLOCATE_PL_WALK()      allocate_pl<Walk_pl, WALK_PL>()
#define CONSTRUCT_PL_WALK()     construct_pl<Walk_pl, WALK_PL>()
#define DEALLOCATE_PL_WALK()    deallocate_pl<Walk_pl, WALK_PL>()
#define ATTACH_PL_WALK()        attach_pl<Walk_pl, WALK_PL>()
#define QUERY_PL_WALK(obj)      query_pl<Walk_pl, WALK_PL>(obj)
#else
#define INIT_PL_WALK()
#define ALLOCATE_PL_WALK()
#define CONSTRUCT_PL_WALK()
#define DEALLOCATE_PL_WALK()
#define ATTACH_PL_WALK()
#define QUERY_PL_WALK(obj)
#endif

// Landmarks (vertices)
#if (LM_PL_ENABLED)
#define INIT_PL_LM()            init_pl<Lm_pl, LM_PL>("Landmarks (vertices)");
#define ALLOCATE_PL_LM()        allocate_pl<Lm_pl, LM_PL>()
#define CONSTRUCT_PL_LM()       construct_pl<Lm_pl, LM_PL>()
#define DEALLOCATE_PL_LM()      deallocate_pl<Lm_pl, LM_PL>()
#define ATTACH_PL_LM()          attach_pl<Lm_pl, LM_PL>()
#define QUERY_PL_LM(obj)        query_pl<Lm_pl, LM_PL>(obj)
#else
#define INIT_PL_LM()
#define ALLOCATE_PL_LM()
#define CONSTRUCT_PL_LM()
#define DEALLOCATE_PL_LM()
#define ATTACH_PL_LM()
#define QUERY_PL_LM(obj)
#endif

// Landmarks random
#if (LM_RANDOM_PL_ENABLED)
#define INIT_PL_LM_RANDOM()       \
  init_pl<Lm_random_pl, LM_RANDOM_PL>("Landmarks random");
#define ALLOCATE_PL_LM_RANDOM()   \
  allocate_pl<Lm_random_pl, LM_RANDOM_PL>()

#define CONSTRUCT_PL_LM_RANDOM()                \
  m_random_g = new Random_lm_generator(*m_arr); \
  construct_pl<Lm_random_pl, LM_RANDOM_PL>(m_random_g)

#define DEALLOCATE_PL_LM_RANDOM()               \
  if (m_random_g) {                             \
    delete m_random_g;                          \
    m_random_g = NULL;                          \
  }                                             \
  deallocate_pl<Lm_random_pl, LM_RANDOM_PL>()

#define ATTACH_PL_LM_RANDOM()                   \
  m_random_g = new Random_lm_generator(*m_arr); \
  attach_pl<Lm_random_pl, LM_RANDOM_PL>(m_random_g)

#define QUERY_PL_LM_RANDOM(obj)   query_pl<Lm_random_pl, LM_RANDOM_PL>(obj)
#else
#define INIT_PL_LM_RANDOM()
#define ALLOCATE_PL_LM_RANDOM()
#define CONSTRUCT_PL_LM_RANDOM()
#define DEALLOCATE_PL_LM_RANDOM()
#define ATTACH_PL_LM_RANDOM()
#define QUERY_PL_LM_RANDOM(obj)
#endif

// Landmarks grid
#if (LM_GRID_PL_ENABLED)
#define INIT_PL_LM_GRID()       init_pl<Lm_grid_pl, LM_GRID_PL>("Landmarks grid")
#define ALLOCATE_PL_LM_GRID()   allocate_pl<Lm_grid_pl, LM_GRID_PL>()

#define CONSTRUCT_PL_LM_GRID()              \
  m_grid_g = new Grid_lm_generator(*m_arr); \
  construct_pl<Lm_grid_pl, LM_GRID_PL>(m_grid_g)

#define DEALLOCATE_PL_LM_GRID() \
  if (m_grid_g) {               \
    delete m_grid_g;            \
    m_grid_g = NULL;            \
  }                             \
  deallocate_pl<Lm_grid_pl, LM_GRID_PL>()

#define ATTACH_PL_LM_GRID()                 \
  m_grid_g = new Grid_lm_generator(*m_arr); \
  attach_pl<Lm_grid_pl, LM_GRID_PL>(m_grid_g)

#define QUERY_PL_LM_GRID(obj)   query_pl<Lm_grid_pl, LM_GRID_PL>(obj)
#else
#define INIT_PL_LM_GRID()
#define ALLOCATE_PL_LM_GRID()
#define CONSTRUCT_PL_LM_GRID()
#define DEALLOCATE_PL_LM_GRID()
#define ATTACH_PL_LM_GRID()
#define QUERY_PL_LM_GRID(obj)
#endif

// Landmarks Halton
#if (LM_HALTON_PL_ENABLED)
#define INIT_PL_LM_HALTON()       \
  init_pl<Lm_halton_pl, LM_HALTON_PL>("Landmarks Halton")
#define ALLOCATE_PL_LM_HALTON()   allocate_pl<Lm_halton_pl, LM_HALTON_PL>()

#define CONSTRUCT_PL_LM_HALTON()                \
  m_halton_g = new Halton_lm_generator(*m_arr); \
  construct_pl<Lm_halton_pl, LM_HALTON_PL>(m_halton_g)

#define DEALLOCATE_PL_LM_HALTON() \
    if (m_halton_g) {             \
    delete m_halton_g;            \
    m_halton_g = NULL;            \
  }                               \
  deallocate_pl<Lm_halton_pl, LM_HALTON_PL>()

#define ATTACH_PL_LM_HALTON()                   \
  m_halton_g = new Halton_lm_generator(*m_arr); \
  attach_pl<Lm_halton_pl, LM_HALTON_PL>(m_halton_g)

#define QUERY_PL_LM_HALTON(obj)   query_pl<Lm_halton_pl, LM_HALTON_PL>(obj)
#else
#define INIT_PL_LM_HALTON()
#define ALLOCATE_PL_LM_HALTON()
#define CONSTRUCT_PL_LM_HALTON()
#define DEALLOCATE_PL_LM_HALTON()
#define ATTACH_PL_LM_HALTON()
#define QUERY_PL_LM_HALTON(obj)
#endif

// Landmarks middle edges
#if (LM_MIDDLE_EDGES_PL_ENABLED)
#define INIT_PL_LM_MIDDLE_EDGES()       \
  init_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>("Landmarks middle edges")
#define ALLOCATE_PL_LM_MIDDLE_EDGES()   \
  allocate_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>()

#define CONSTRUCT_PL_LM_MIDDLE_EDGES()                   \
  m_middle_edges_g = new Middle_edges_generator(*m_arr); \
  construct_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(m_middle_edges_g)

#define DEALLOCATE_PL_LM_MIDDLE_EDGES() \
  if (m_middle_edges_g) {               \
    delete m_middle_edges_g;            \
    m_middle_edges_g = NULL;            \
  }                                     \
  deallocate_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>()

#define ATTACH_PL_LM_MIDDLE_EDGES()                      \
  m_middle_edges_g = new Middle_edges_generator(*m_arr); \
  attach_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(m_middle_edges_g)

#define QUERY_PL_LM_MIDDLE_EDGES(obj)   \
  query_pl<Lm_middle_edges_pl, LM_MIDDLE_EDGES_PL>(obj)
#else
#define INIT_PL_LM_MIDDLE_EDGES()
#define ALLOCATE_PL_LM_MIDDLE_EDGES()
#define CONSTRUCT_PL_LM_MIDDLE_EDGES()
#define DEALLOCATE_PL_LM_MIDDLE_EDGES()
#define ATTACH_PL_LM_MIDDLE_EDGES()
#define QUERY_PL_LM_MIDDLE_EDGES(obj)
#endif

// Landmarks specified points
#if (LM_SPECIFIED_POINTS_PL_ENABLED)
#define INIT_PL_LM_SPECIFIED_POINTS()       \
  init_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>("Landmarks specified points");
#define ALLOCATE_PL_LM_SPECIFIED_POINTS()   \
  allocate_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>()

#define CONSTRUCT_PL_LM_SPECIFIED_POINTS()                       \
  m_specified_points_g = new Specified_points_generator(*m_arr); \
  construct_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>(m_specified_points_g)

#define DEALLOCATE_PL_LM_SPECIFIED_POINTS() \
  if (m_specified_points_g) {               \
    delete m_specified_points_g;            \
    m_specified_points_g = NULL;            \
  }                                         \
  deallocate_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>()

#define ATTACH_PL_LM_SPECIFIED_POINTS()                          \
  m_specified_points_g = new Specified_points_generator(*m_arr); \
  attach_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>(m_specified_points_g)

#define QUERY_PL_LM_SPECIFIED_POINTS(obj)   \
  query_pl<Lm_specified_points_pl, LM_SPECIFIED_POINTS_PL>(obj)
#else
#define INIT_PL_LM_SPECIFIED_POINTS()
#define ALLOCATE_PL_LM_SPECIFIED_POINTS()
#define CONSTRUCT_PL_LM_SPECIFIED_POINTS()
#define DEALLOCATE_PL_LM_SPECIFIED_POINTS()
#define ATTACH_PL_LM_SPECIFIED_POINTS()
#define QUERY_PL_LM_SPECIFIED_POINTS(obj)
#endif

// Triangulation_pl
#if (TRIANGULATION_PL_ENABLED)
#define INIT_PL_TRIANGULATION_PL()       \
  init_pl<Triangulation_pl, TRIANGULATION_PL>("Triangulation")
#define ALLOCATE_PL_TRIANGULATION_PL()   \
  allocate_pl<Triangulation_pl, TRIANGULATION_PL>()
#define CONSTRUCT_PL_TRIANGULATION_PL()  \
  construct_pl<Triangulation_pl, TRIANGULATION_PL>()
#define DEALLOCATE_PL_TRIANGULATION_PL() \
  deallocate_pl<Triangulation_pl, TRIANGULATION_PL>()
#define ATTACH_PL_TRIANGULATION_PL()     \
  attach_pl<Triangulation_pl, TRIANGULATION_PL>()
#define QUERY_PL_TRIANGULATION_PL(obj)   \
  query_pl<Triangulation_pl, TRIANGULATION_PL>(obj);
#else
#define INIT_PL_TRIANGULATION_PL()
#define ALLOCATE_PL_TRIANGULATION_PL()
#define CONSTRUCT_PL_TRIANGULATION_PL()
#define DEALLOCATE_PL_TRIANGULATION_PL()
#define ATTACH_PL_TRIANGULATION_PL()
#define QUERY_PL_TRIANGULATION_PL(obj)
#endif

// Trapezoidal RIC
#if (TRAPEZOID_RIC_PL_ENABLED)
#define INIT_PL_TRAPEZOID_RIC_PL()       \
  init_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>("Trapezoidal RIC")
#define ALLOCATE_PL_TRAPEZOID_RIC_PL()   \
  allocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>()
#define CONSTRUCT_PL_TRAPEZOID_RIC_PL()  \
  construct_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>()
#define DEALLOCATE_PL_TRAPEZOID_RIC_PL() \
  deallocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>()
#define ATTACH_PL_TRAPEZOID_RIC_PL()     \
  attach_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>()
#define QUERY_PL_TRAPEZOID_RIC_PL(obj)   \
  query_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_PL>(obj)
#else
#define INIT_PL_TRAPEZOID_RIC_PL()
#define ALLOCATE_PL_TRAPEZOID_RIC_PL()
#define CONSTRUCT_PL_TRAPEZOID_RIC_PL()
#define DEALLOCATE_PL_TRAPEZOID_RIC_PL()
#define ATTACH_PL_TRAPEZOID_RIC_PL()
#define QUERY_PL_TRAPEZOID_RIC_PL(obj)
#endif

#if (TRAPEZOID_RIC_NO_GUARANTEE_PL_ENABLED)
#define INIT_PL_TRAPEZOID_RIC_NO_GUARANTEE()       \
  init_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>("Trapezoidal RIC without guarantees")
#define ALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE()   \
  allocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>()
#define CONSTRUCT_PL_TRAPEZOID_RIC_NO_GUARANTEE()  \
  construct_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>()
#define DEALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE() \
  deallocate_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>()

#define ATTACH_PL_TRAPEZOID_RIC_NO_GUARANTEE()                           \
  Pl_variant& var = m_locators[TRAPEZOID_RIC_NO_GUARANTEE_PL].m_variant; \
  Trapezoid_ric_pl* strategy = std::get<Trapezoid_ric_pl*>(var);       \
  strategy->with_guarantees(false);                                      \
  attach_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>()

#define QUERY_PL_TRAPEZOID_RIC_NO_GUARANTEE(obj)   \
  query_pl<Trapezoid_ric_pl, TRAPEZOID_RIC_NO_GUARANTEE_PL>(obj)
#else
#define INIT_PL_TRAPEZOID_RIC_NO_GUARANTEE()
#define ALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE()
#define CONSTRUCT_PL_TRAPEZOID_RIC_NO_GUARANTEE()
#define DEALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE()
#define ATTACH_PL_TRAPEZOID_RIC_NO_GUARANTEE()
#define QUERY_PL_TRAPEZOID_RIC_NO_GUARANTEE(obj)
#endif

#define MEASURE_NAIVE_PL(timer, op) measure<NAIVE_PL>(timer, [&](){ op; });
#define MEASURE_SIMPLE_PL(timer, op) measure<SIMPLE_PL>(timer, [&](){ op; });
#define MEASURE_WALK_PL(timer, op) measure<WALK_PL>(timer, [&](){ op; });
#define MEASURE_LM_PL(timer, op) measure<LM_PL>(timer, [&](){ op; });
#define MEASURE_LM_RANDOM_PL(timer, op) \
  measure<LM_RANDOM_PL>(timer, [&](){ op; });
#define MEASURE_LM_GRID_PL(timer, op) \
  measure<LM_GRID_PL>(timer, [&](){ op; });
#define MEASURE_LM_HALTON_PL(timer, op) \
  measure<LM_HALTON_PL>(timer, [&](){ op; });
#define MEASURE_LM_MIDDLE_EDGES_PL(timer, op) \
  measure<LM_MIDDLE_EDGES_PL>(timer, [&](){ op; });
#define MEASURE_LM_SPECIFIED_POINTS_PL(timer, op) \
  measure<LM_SPECIFIED_POINTS_PL>(timer, [&](){ op; });
#define MEASURE_TRIANGULATION_PL(timer, op) \
  measure<TRIANGULATION_PL>(timer, [&](){ op; });
#define MEASURE_TRAPEZOID_RIC_PL(timer, op) \
  measure<TRAPEZOID_RIC_PL>(timer, [&](){ op; });
#define MEASURE_TRAPEZOID_RIC_NO_GUARANTEE_PL(timer, op) \
  measure<TRAPEZOID_RIC_NO_GUARANTEE_PL>(timer, [&](){ op; });

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
  m_locators.resize(NUM_PL_STRATEGIES);

  INIT_PL_NATIVE();
  INIT_PL_SIMPLE();
  INIT_PL_WALK();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  INIT_PL_LM();
  INIT_PL_LM_RANDOM();
  INIT_PL_LM_GRID();
  INIT_PL_LM_HALTON();
  INIT_PL_LM_MIDDLE_EDGES();
  INIT_PL_LM_SPECIFIED_POINTS();

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  INIT_PL_TRIANGULATION_PL();
#endif

  INIT_PL_TRAPEZOID_RIC_PL();
  INIT_PL_TRAPEZOID_RIC_NO_GUARANTEE();
}

//! Set the file names.
template <typename GeomTraits, typename TopolTraits>
void Point_location_test<GeomTraits, TopolTraits>::
set_filenames(const char* points_filename, const char* xcurves_filename,
              const char* curves_filename, const char* queries_filename)
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
  DEALLOCATE_PL_NATIVE();
  DEALLOCATE_PL_SIMPLE();
  DEALLOCATE_PL_WALK();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  DEALLOCATE_PL_LM();
  DEALLOCATE_PL_LM_RANDOM();
  DEALLOCATE_PL_LM_GRID();
  DEALLOCATE_PL_LM_HALTON();
  DEALLOCATE_PL_LM_MIDDLE_EDGES();
  DEALLOCATE_PL_LM_SPECIFIED_POINTS();

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  DEALLOCATE_PL_TRIANGULATION_PL();
#endif

  DEALLOCATE_PL_TRAPEZOID_RIC_PL();
  DEALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE();
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
{ if (m_arr) m_arr->clear(); }

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::allocate_pl_strategies()
{
  // Allocate all point location strategies.
  ALLOCATE_PL_NATIVE();
  ALLOCATE_PL_SIMPLE();
  ALLOCATE_PL_WALK();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  ALLOCATE_PL_LM();
  ALLOCATE_PL_LM_RANDOM();
  ALLOCATE_PL_LM_GRID();
  ALLOCATE_PL_LM_HALTON();
  ALLOCATE_PL_LM_MIDDLE_EDGES();
  ALLOCATE_PL_LM_SPECIFIED_POINTS();

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  ALLOCATE_PL_TRIANGULATION_PL();
#endif

  ALLOCATE_PL_TRAPEZOID_RIC_PL();
  ALLOCATE_PL_TRAPEZOID_RIC_NO_GUARANTEE();

  return true;
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::construct_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  CONSTRUCT_PL_NATIVE();
  CONSTRUCT_PL_SIMPLE();
  CONSTRUCT_PL_WALK();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  MEASURE_LM_PL(timer, CONSTRUCT_PL_LM());
  MEASURE_LM_RANDOM_PL(timer, CONSTRUCT_PL_LM_RANDOM());
  MEASURE_LM_GRID_PL(timer, CONSTRUCT_PL_LM_GRID());
  MEASURE_LM_HALTON_PL(timer, CONSTRUCT_PL_LM_HALTON());
  MEASURE_LM_MIDDLE_EDGES_PL(timer, CONSTRUCT_PL_LM_MIDDLE_EDGES());
  MEASURE_LM_SPECIFIED_POINTS_PL(timer, CONSTRUCT_PL_LM_SPECIFIED_POINTS());

#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  MEASURE_TRIANGULATION_PL(timer, CONSTRUCT_PL_TRIANGULATION_PL());
#endif

  MEASURE_TRAPEZOID_RIC_PL(timer, CONSTRUCT_PL_TRAPEZOID_RIC_PL());
  MEASURE_TRAPEZOID_RIC_NO_GUARANTEE_PL(timer,
                                        CONSTRUCT_PL_TRAPEZOID_RIC_NO_GUARANTEE());

  std::cout << std::endl;

  return true;
}

template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::attach_pl_strategies()
{
  // Initialize all point location strategies.
  CGAL::Timer timer;

  ATTACH_PL_NATIVE();
  ATTACH_PL_SIMPLE();
  ATTACH_PL_WALK();

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  MEASURE_LM_PL(timer, ATTACH_PL_LM());
  MEASURE_LM_RANDOM_PL(timer, ATTACH_PL_LM_RANDOM());
  MEASURE_LM_GRID_PL(timer, ATTACH_PL_LM_GRID());
  MEASURE_LM_HALTON_PL(timer, ATTACH_PL_LM_HALTON());
  MEASURE_LM_MIDDLE_EDGES_PL(timer, ATTACH_PL_LM_MIDDLE_EDGES());
  MEASURE_LM_SPECIFIED_POINTS_PL(timer, ATTACH_PL_LM_SPECIFIED_POINTS());
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  MEASURE_TRIANGULATION_PL(timer, ATTACH_PL_TRIANGULATION_PL());
#endif

  MEASURE_TRAPEZOID_RIC_PL(timer, ATTACH_PL_TRAPEZOID_RIC_PL());
  MEASURE_TRAPEZOID_RIC_NO_GUARANTEE_PL(timer,
                                        ATTACH_PL_TRAPEZOID_RIC_NO_GUARANTEE());

  std::cout << std::endl;

  return true;
}

// Perform the test
template <typename GeomTraits, typename TopolTraits>
bool Point_location_test<GeomTraits, TopolTraits>::perform()
{
  Variants_vector objs[NUM_PL_STRATEGIES];

  // Locate the points in the list using all point location strategies.

  // std::cout << "Time in seconds" << std::endl;
  std::cout << std::endl;

  size_t pl_index = 0;

  QUERY_PL_NATIVE(objs[pl_index++]);
  QUERY_PL_SIMPLE(objs[pl_index++]);
  QUERY_PL_WALK(objs[pl_index++]);

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS)

  QUERY_PL_LM(objs[pl_index++]);
  QUERY_PL_LM_RANDOM(objs[pl_index++]);
  QUERY_PL_LM_GRID(objs[pl_index++]);
  QUERY_PL_LM_HALTON(objs[pl_index++]);
  QUERY_PL_LM_MIDDLE_EDGES(objs[pl_index++]);
  QUERY_PL_LM_SPECIFIED_POINTS(objs[pl_index++]);
#endif

#if (TEST_GEOM_TRAITS == SEGMENT_GEOM_TRAITS)
  QUERY_PL_TRIANGULATION_PL(objs[pl_index++]);
#endif

  QUERY_PL_TRAPEZOID_RIC_PL(objs[pl_index++]);
  QUERY_PL_TRAPEZOID_RIC_NO_GUARANTEE(objs[pl_index++]);

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
            std::cout << "Expected: a face." << std::endl;
            std::cout << "Actual: a different face" << std::endl;
            result += -1;
          }
          continue;
        }

        std::cout << "Error: point location number " << pl << std::endl;
        std::cout << "Expected: a face." << std::endl;
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
    Face_const_handle* fh_ref = std::get_if<Face_const_handle>(&(objs[0][qi]));
    if (fh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        Face_const_handle* fh_cur =
          std::get_if<Face_const_handle>(&(objs[pl][qi]));
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
        Halfedge_const_handle* hh_cur =
          std::get_if<Halfedge_const_handle>(&(objs[pl][qi]));
        if (hh_cur) {
          std::cout << "Actual: a halfedge." << std::endl;
          continue;
        }
        Vertex_const_handle* vh_cur =
          std::get_if<Vertex_const_handle>(&(objs[pl][qi]));
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
    Halfedge_const_handle* hh_ref =
      std::get_if<Halfedge_const_handle>(&(objs[0][qi]));
    if (hh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        Halfedge_const_handle* hh_cur =
          std::get_if<Halfedge_const_handle>(&(objs[pl][qi]));
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
        Face_const_handle* fh_cur =
          std::get_if<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        Vertex_const_handle* vh_cur =
          std::get_if<Vertex_const_handle>(&(objs[pl][qi]));
        if (vh_cur) {
          std::cout << "Actual: a vertex." << std::endl;
          continue;
        }
        std::cout << "Actual: an unknowen object." << std::endl;
      }
      continue;
    }

    // Assign object to a vertex
    Vertex_const_handle* vh_ref =
      std::get_if<Vertex_const_handle>(&(objs[0][qi]));
    if (vh_ref) {
      for (size_t pl = 1; pl < pls_num; ++pl) {
        Vertex_const_handle* vh_cur =
          std::get_if<Vertex_const_handle>(&(objs[pl][qi]));
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
        Face_const_handle* fh_cur =
          std::get_if<Face_const_handle>(&(objs[pl][qi]));
        if (fh_cur) {
          std::cout << "Actual: a face." << std::endl;
          continue;
        }
        Halfedge_const_handle* hh_cur =
          std::get_if<Halfedge_const_handle>(&(objs[pl][qi]));
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
