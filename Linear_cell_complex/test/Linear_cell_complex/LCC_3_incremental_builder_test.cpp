#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_incremental_builder_3.h>
#include <CGAL/Combinatorial_map_save_load.h>

#include "Linear_cell_complex_3_test.h"

///////////////////////////////////////////////////////////////////////////////
/*     3
 *    /|\
 *   0-|-2
 *    \|/
 *     1
 */
template<typename IncrementalBuilder>
void make_tetrahedron_with_builder(IncrementalBuilder& ib,
                                   std::size_t i0,
                                   std::size_t i1,
                                   std::size_t i2,
                                   std::size_t i3)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2});
  ib.add_facet({i1,i0,i3});
  ib.add_facet({i2,i1,i3});
  ib.add_facet({i0,i2,i3});
  ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      4
 *     /|\
 *    0-|-3
 *    | | |
 *    1---2
 */
template<typename IncrementalBuilder>
void make_pyramid_with_builder(IncrementalBuilder& ib,
                               std::size_t i0,
                               std::size_t i1,
                               std::size_t i2,
                               std::size_t i3,
                               std::size_t i4)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3});
  ib.add_facet({i1,i0,i4});
  ib.add_facet({i2,i1,i4});
  ib.add_facet({i3,i2,i4});
  ib.add_facet({i0,i3,i4});
  ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      3
 *     /|\
 *    4---5
 *    | | |
 *    | 0 |
 *    |/ \|
 *    1---2
 */
template<typename IncrementalBuilder>
void make_prism_with_builder(IncrementalBuilder& ib,
                             std::size_t i0,
                             std::size_t i1,
                             std::size_t i2,
                             std::size_t i3,
                             std::size_t i4,
                             std::size_t i5)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2});
  ib.add_facet({i1,i0,i3,i4});
  ib.add_facet({i2,i1,i4,i5});
  ib.add_facet({i0,i2,i5,i3});
  ib.add_facet({i5,i4,i3});
  ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
/*      7----6
 *     /|   /|
 *    4----5 |
 *    | 3--|-2
 *    |/   |/
 *    0----1
 */
template<typename IncrementalBuilder>
void make_hexahedron_with_builder(IncrementalBuilder& ib,
                                  std::size_t i0,
                                  std::size_t i1,
                                  std::size_t i2,
                                  std::size_t i3,
                                  std::size_t i4,
                                  std::size_t i5,
                                  std::size_t i6,
                                  std::size_t i7)
{
  ib.begin_surface();
  ib.add_facet({i0,i1,i2,i3});
  ib.add_facet({i1,i0,i4,i5});
  ib.add_facet({i2,i1,i5,i6});
  ib.add_facet({i3,i2,i6,i7});
  ib.add_facet({i0,i3,i7,i4});
  ib.add_facet({i7,i6,i5,i4});
  ib.end_surface();
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
bool test_ib(const char* filename)
{
  typedef typename LCC::Point Point;
  LCC lcc;
  CGAL::Linear_cell_complex_incremental_builder_3<LCC> ib(lcc);

  ib.add_vertex(Point(0,0,0)); // 0
  ib.add_vertex(Point(1,0,0)); // 1
  ib.add_vertex(Point(1,2,0)); // 2
  ib.add_vertex(Point(0,2,0)); // 3

  ib.add_vertex(Point(0,0,1.5)); // 4
  ib.add_vertex(Point(1,0,1.5)); // 5
  ib.add_vertex(Point(1,2,1.5)); // 6
  ib.add_vertex(Point(0,2,1.5)); // 7

  ib.add_vertex(Point(0.5,1,2.5)); // 8

  ib.add_vertex(Point(2,0,0)); // 9
  ib.add_vertex(Point(2,0,1)); // 10

  trace_test_begin();
  make_hexahedron_with_builder(ib, 0,1,2,3,4,5,6,7);
  make_pyramid_with_builder(ib, 4,5,6,7,8);
  make_prism_with_builder(ib, 2,1,9,6,5,10);
  make_tetrahedron_with_builder(ib, 6,5,10,8);

  if ( !check_number_of_cells_3(lcc, 11, 22, 16, 4, 1) )
  { return false; }

  LCC lcc2;
  std::ifstream input(std::string("data/")+filename);
  if (!input)
  {
    std::cout<<"Problem to load LCC data/"<<filename<<std::endl;
    return false;
  }
  input>>lcc2;
  input.close();

  if (!lcc.is_isomorphic_to(lcc2))
  { return false; }
  trace_test_end();

  return true;
}

int main()
{
  std::cout<<"LCC_3_incremental_builder_test (v1)."<<std::flush;

  // ****************** TEST FOR CMAP ******************
  trace_display_msg("test_LCC_3<LCC3>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC3;
  if ( !test_ib<LCC3>("lcc3_ib_test.cmap") )
  {
    std::cout<<" Error during Test_LCC_3<LCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  // ****************** TEST FOR GMAP ******************
  trace_display_msg("test_LCC_3<GLCC3>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<3> GLCC3;
  if ( !test_ib<GLCC3>("lcc3_ib_test.gmap") )
  {
    std::cout<<" Error during Test_LCC_3<GLCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
