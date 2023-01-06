#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <vector>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>   LCC_3_gmap;

template<typename LCC,
         typename Map=typename LCC::Combinatorial_data_structure>
struct Alpha1
{
  static typename LCC::Dart_descriptor run(LCC&, typename LCC::Dart_descriptor d)
  { return d; }
};
template<typename LCC>
struct Alpha1<LCC, CGAL::Generalized_map_tag>
{
  static typename LCC::Dart_descriptor run(LCC& lcc, typename LCC::Dart_descriptor d)
  { return lcc.template alpha<1>(d); }
};

template<typename LCC>
void run_test()
{
  typedef typename LCC::Point Point;
  typedef typename LCC::Dart_descriptor Dart_descriptor;

  LCC lcc;

  Dart_descriptor d1 = lcc.
    make_hexahedron(Point(0,0,0),Point(1,0,0),
                    Point(1,2,0),Point(0,2,0),
                    Point(0,3,4),Point(0,0,4),
                    Point(6,0,4),Point(6,3,4));
  Dart_descriptor d2 = lcc.
    make_hexahedron(Point(0,-5,0),Point(2,-5,0),
                    Point(2,-2,0),Point(0,-2,0),
                    Point(1,-1,5),Point(1,-2,5),
                    Point(5,-2,5),Point(5,-2,5));
  Dart_descriptor d3 = lcc.
    make_hexahedron(Point(1,0,5),Point(0,0,6),
                    Point(0,2,5),Point(1,2,6),
                    Point(1,3,8),Point(0,0,8),
                    Point(5,0,9),Point(7,3,9));

  lcc.template sew<3>(d1,lcc.other_orientation
                      (lcc.template opposite<2>
                       (lcc.next(lcc.next(lcc.template opposite<2>(d2))))));
  lcc.template sew<3>(lcc.template opposite<2>(lcc.next(d1)),
                      lcc.other_orientation(lcc.template opposite<2>(lcc.previous(d3))));

  lcc.insert_cell_1_in_cell_2(lcc.next(d1),
                              Alpha1<LCC>::run(lcc, lcc.previous(d1)));
  d2=lcc.template opposite<2>(lcc.next(lcc.next
                                        (lcc.template opposite<2>(d1))));
  lcc.insert_cell_1_in_cell_2(d2, Alpha1<LCC>::run
                              (lcc, lcc.next(lcc.next(d2))));

  std::vector<Dart_descriptor> path;
  path.push_back(lcc.next(d1));
  path.push_back(lcc.next(lcc.template opposite<2>(lcc.previous(d1))));
  path.push_back(lcc.previous(d2));
  path.push_back(lcc.next(lcc.template opposite<2>(d2)));
  lcc.insert_cell_2_in_cell_3(path.begin(),path.end());

  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;
}

int main()
{
  run_test<LCC_3_cmap>();
  run_test<LCC_3_gmap>();

  return EXIT_SUCCESS;
}
