
#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Arr_segment_traits_2.h>
#include<CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>

#include "utils.h"

typedef CGAL::Exact_rational                          Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Base_traits_2;
typedef Base_traits_2::Curve_2                        Base_curve_2;
typedef Base_traits_2::X_monotone_curve_2             Base_x_monotone_curve_2;
typedef CGAL::Arr_curve_data_traits_2<Base_traits_2,
                                      unsigned int,
                                      std::plus<unsigned int> >  
                                                      Traits_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_iterator                Vertex_iterator;
typedef Arrangement_2::Edge_iterator                  Edge_iterator;
typedef std::vector<Curve_2>                          Curve_container;



int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cout << "Missing input file" << std::endl;
    return 0;
  }

  std::ifstream in_file(argv[1]);
  if (!in_file.is_open()) {
    std::cout << "Can't open file: " << argv[1] << std::endl;
    return -1;
  }

  if (argc < 3) {
    std::cout << "Missing output file" << std::endl;
    return 0;
  }
  std::ifstream test (argv[2]);
  if (test.is_open()) {
    std::cout << argv[2] << " already exist, overwrite? (y/n)" << std::endl;
    char c = std::cin.get();
    if (c != 'y') return 0;
  }
  std::ofstream out_file(argv[2]);
  if (!out_file.is_open()) {
    std::cout << "Can't open file: " << argv[2] << std::endl;
    return -1;
  }

  Traits_2 traits;
  Arrangement_2 arr(&traits);
  Curve_container curves;

  unsigned int num_of_curves;
  in_file >> num_of_curves;
  out_file << num_of_curves << std::endl;
  curves.resize(num_of_curves);
  for (size_t i = 0; i < num_of_curves; ++i) {
    Base_curve_2 base_cv;
    in_file >> base_cv;
    out_file << base_cv << std::endl;
    curves[i] = Curve_2(base_cv, 1);
  }
  out_file << std::endl;
  std::cout << "Inserting the curves to the map ... " << std::endl;

  CGAL::insert(arr, curves.begin(), curves.end());

  std::cout<< "Finished insertion...\n";
  std::cout<<CGAL::is_valid(arr)<<"\n";

  std::cout<<"|V| = " << arr.number_of_vertices() << std::endl;
  std::cout<<"|E| = " << arr.number_of_edges() << std::endl;
  std::cout<<"|F| = " << arr.number_of_faces() << std::endl;

  out_file << arr.number_of_vertices() << std::endl;
  out_file << arr.number_of_edges() << std::endl;
  out_file << arr.number_of_faces() << std::endl;
  out_file << std::endl;

  std::vector<Point_2> pts(arr.number_of_vertices());
  size_t i = 0;
  for (Vertex_iterator vit = arr.vertices_begin(); vit != arr.vertices_end();
       ++vit, ++i)
    pts[i] = vit->point();
  std::sort(pts.begin(), pts.end(), Point_compare<Traits_2>(traits));
  std::copy(pts.begin(), pts.end(),
            std::ostream_iterator<Point_2>(out_file, "\n"));
  out_file << std::endl;

  std::vector<X_monotone_curve_2> curves_res(arr.number_of_edges());
  i = 0;
  for (Edge_iterator eit = arr.edges_begin(); eit != arr.edges_end();
       ++eit, ++i)
    curves_res[i] = eit->curve();
  std::sort(curves_res.begin(), curves_res.end(),
            Curve_compare<Traits_2>(traits));

  for (i = 0; i < curves_res.size(); ++i)
    out_file << curves_res[i] <<' '<<curves_res[i].data()<< std::endl;

  std::cout << argv[2]
            << " was generated successfully"
            << ", dont forget to add it to test_construction.cmd"
            << std::endl;
  return 0;
}
