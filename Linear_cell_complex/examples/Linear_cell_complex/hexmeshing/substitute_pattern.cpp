#include "hexmeshing.h"
#include "utils.h"

int pattern_substitution()
{
  LCC lcc;
  auto d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,5), Point(0,0,5),
                        Point(5,0,5), Point(5,5,5));

  auto d2=
    lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
                        Point(10,5,0), Point(5,5,0),
                        Point(5,5,5), Point(5,0,5),
                        Point(10,0,5), Point(10,5,5));

  auto d3=
    lcc.make_hexahedron(Point(0,0,-5), Point(5,0,-5),
                        Point(5,5,-5), Point(0,5,-5),
                        Point(0,5,0), Point(0,0,0),
                        Point(5,0,0), Point(5,5,0));

  lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));
  lcc.sew<3>(lcc.beta(d1, 0, 2), lcc.beta(d3, 1, 2));

  lcc.display_characteristics(std::cout);

  auto toprocess_mark = lcc.get_new_mark();
  debug_edge_mark = lcc.get_new_mark();

  lcc.mark_cell<0>(d1, toprocess_mark);
  auto md1 = lcc.beta(d1, 1);
  lcc.mark_cell<0>(md1, toprocess_mark);
  auto md2 = lcc.beta(d2, 1);
  lcc.mark_cell<0>(md2, toprocess_mark);
  auto md3 = lcc.beta(d2, 1, 1);
  lcc.mark_cell<0>(md3, toprocess_mark);
  auto md4 = d3;
  lcc.mark_cell<0>(md4, toprocess_mark);

  // lcc.mark_cell<0>(md3, toprocess_mark);

  std::vector<Dart_handle> marked_cells;

  // We assume we already filled out marked_cells earlier

  marked_cells.push_back(d1);
  marked_cells.push_back(md1);
  marked_cells.push_back(md2);
  marked_cells.push_back(md3);
  marked_cells.push_back(md4);

  CGAL::HexRefinement::tworefinement_one_step_from_marked(lcc, marked_cells, toprocess_mark);

  render(lcc, toprocess_mark, debug_edge_mark);
  return EXIT_SUCCESS;
}

int grid_test() {
  using namespace CGAL::HexRefinement::TwoRefinement;

  LCC lcc;
  auto grid = generate_grid(lcc, Point(0,0,0), 5, 10);
  
  lcc.display_characteristics(std::cout);

  auto m1 = lcc.get_new_mark();
  auto m2 = lcc.get_new_mark();

  mark_nodes_in_cell_pair(lcc, grid, Plane::ZX, m1, m1);

  // auto a =  lcc.beta(lcc.first_dart(), 0, 2, 0);
  // mark_face(lcc, a, m2);
  // lcc.mark_cell<0>(a, m1);

  render(lcc, m1, m2);

  return EXIT_SUCCESS;
}


int main(){
  return grid_test();
}
