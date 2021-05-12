#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/iterator.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>


// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational         Number_type;

typedef CGAL::Simple_cartesian<Number_type> Kernel;

typedef CGAL::Polygon_2<Kernel>                        Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>             Polygon_with_holes_2;
typedef CGAL::Gps_segment_traits_2<Kernel>             Traits_2;

typedef std::vector<Polygon_2>               Polygons_vec;
typedef std::vector<Polygon_with_holes_2>    Polygons_with_holes_vec;
typedef std::back_insert_iterator<Polygons_with_holes_vec>
                                             Output_itr;

template <class Vec>
bool are_polygons_valid(Vec& vec)
{
  Traits_2 tr;
  unsigned int i=0;
  for(; i < vec.size(); ++i)
  {
    if(!is_valid_unknown_polygon(vec[i], tr))
      return false;
  }
  return true;
}

void write_result_to_file(Polygons_with_holes_vec& result,
                          std::ofstream& out)
{
  std::ostream_iterator<Polygon_with_holes_2> ostream_itr(out, "\n");
  out << result.size() << std::endl;
  std::copy(result.begin(), result.end(), ostream_itr);
  out<<std::endl;
  result.clear();
}

int main(int argc, char *argv[])
{
  if(argc<3)
  {
    std::cerr << "Missing input file!" << std::endl
              << "When performed as part of the test suite no input file is expected."
              << std::endl;
    return 0;
  }

  std::ifstream inp (argv[1]);
  if(!inp.is_open())
  {
    std::cerr << "Failed to open file!" << std::endl;
    std::exit(-1);
  }

  std::ifstream test (argv[2]);
  if(test.is_open())
  {
      std::cout<<argv[2] << " already exist, overwrite? (y/n)" << std::endl;
      char c = std::cin.get();
      if(c != 'y')
        return 0;

  }
  std::ofstream out (argv[2]);
  if(!out.is_open())
  {
    std::cerr << "Failed to create output file" << std::endl;
    std::exit(-1);
  }

  Polygons_vec            polygons;
  Polygons_with_holes_vec polygons_with_holes;
  Polygons_with_holes_vec result;

  Output_itr oi(result);


  int n;
  // read polygons
  inp >> n;
  out << n << std::endl;
  polygons.resize(n);
  int i=0;
  for(; i<n; ++i)
  {
    inp >> polygons[i];
    out << polygons[i]<< std::endl;
  }

  if(!are_polygons_valid(polygons))
  {
    std::cout << "invalid input polygons !!!" << std::endl;;
  }

  //read polygons with holes
  inp >> n;
  out << n << std::endl;
  polygons_with_holes.resize(n);
  for(i=0 ; i<n; ++i)
  {
    inp >> polygons_with_holes[i];
    out << polygons_with_holes[i] << std::endl;
  }
  if(!are_polygons_valid(polygons_with_holes))
  {
    std::cout<<"invalid input polygon with hole !!!" << std::endl;;
  }

  CGAL::join(polygons.begin(), polygons.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);

  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);

  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);


  CGAL::intersection(polygons.begin(), polygons.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);

  CGAL::intersection(polygons_with_holes.begin(),
                     polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);

  CGAL::intersection(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;;
  }
  write_result_to_file(result, out);


  CGAL::symmetric_difference(polygons.begin(), polygons.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;
  }
  write_result_to_file(result, out);

  CGAL::symmetric_difference(polygons_with_holes.begin(),
                             polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout<<"invalid polygon was created !!!" << std::endl;
  }
  write_result_to_file(result, out);

  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(!are_polygons_valid(result))
  {
    std::cout << "invalid polygon was created !!!" << std::endl;
  }
  write_result_to_file(result, out);

  return 0;
}
