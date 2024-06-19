#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;
typedef LCC::Dart_descriptor Dart_descriptor;
typedef LCC::Point           Point;


void mark_1template(LCC& pattern, typename LCC::size_type mark_to_preserve) {
  std::cout << "mark: " << mark_to_preserve << std::endl;
}

int main()
{
  LCC lcc;
  Pattern_substituer<LCC> ps;
  Dart_descriptor d1=
    lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                        Point(5,5,0), Point(0,5,0),
                        Point(0,5,5), Point(0,0,5),
                        Point(5,0,5), Point(5,5,5));

  ps.load_vpatterns(CGAL::data_file_path("hexmeshing/vpattern"));
  
  if (ps.query_replace_one_volume(lcc, d1) == std::numeric_limits<std::size_t>::max()){
    std::cout << "No replacement was found" << std::endl;
    return EXIT_FAILURE;
  };

  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
