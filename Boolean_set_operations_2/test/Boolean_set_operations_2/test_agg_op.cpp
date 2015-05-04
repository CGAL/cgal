#include <string>

#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <vector>
#include <iostream>
#include <cstdlib>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Polygon_set_2.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational        Number_type;

// instead of
//typedef CGAL::Simple_cartesian<Number_type>            Kernel;
// workaround for VC++ 
struct Kernel : public CGAL::Simple_cartesian<Number_type> {};

typedef CGAL::Polygon_2<Kernel>                       Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>            Polygon_with_holes_2;

typedef CGAL::Gps_segment_traits_2<Kernel>            Traits;
typedef CGAL::Polygon_set_2<Kernel>                   Ps;

typedef CGAL::Arr_segment_traits_2<Kernel>            Arr_traits;
typedef CGAL::Gps_traits_2<Arr_traits>                General_traits;
typedef CGAL::General_polygon_set_2<General_traits>   Gps;

template <class Container>
bool are_equal(const Container& l1, Container& l2)
{
  bool eq = ((l1.size() == l2.size()) &&
             (std::equal(l1.begin(), l1.end(), l2.begin())));
  if(!eq)
  {
    std::ostream_iterator<Polygon_with_holes_2> osi(std::cout, "\n");
    std::cout<<"infile: " << std::endl;
    std::copy(l1.begin(), l1.end(), osi);
    std::cout<<std::endl<<"calculated: " <<std::endl;
    std::copy(l2.begin(), l2.end(), osi);
  }
  l2.clear();

  return (eq);
}

template <class Vector>
void read_polygons_vector(std::istream& inp,
                          Vector& vec)
{
  unsigned int n = 0;
  inp >> n;
  vec.resize(n);
  unsigned int i;
  for(i=0; i<n; ++i)
  {
    inp >> vec[i];
  }

}
void read_file(std::istream&                       inp,
               std::vector<Polygon_2>&             polygons,
               std::vector<Polygon_with_holes_2>&  polygons_with_holes,
               std::vector<Polygon_with_holes_2>&  join1_res,
               std::vector<Polygon_with_holes_2>&  join2_res,
               std::vector<Polygon_with_holes_2>&  join3_res,
               std::vector<Polygon_with_holes_2>&  intersection1_res,
               std::vector<Polygon_with_holes_2>&  intersection2_res,
               std::vector<Polygon_with_holes_2>&  intersection3_res,
               std::vector<Polygon_with_holes_2>&  symm_diff1_res,
               std::vector<Polygon_with_holes_2>&  symm_diff2_res,
               std::vector<Polygon_with_holes_2>&  symm_diff3_res)
{
  read_polygons_vector(inp, polygons);
  read_polygons_vector(inp, polygons_with_holes);
  read_polygons_vector(inp, join1_res);
  read_polygons_vector(inp, join2_res);
  read_polygons_vector(inp, join3_res);
  read_polygons_vector(inp, intersection1_res);
  read_polygons_vector(inp, intersection2_res);
  read_polygons_vector(inp, intersection3_res);
  read_polygons_vector(inp, symm_diff1_res);
  read_polygons_vector(inp, symm_diff2_res);
  read_polygons_vector(inp, symm_diff3_res);
}

bool test_one_file(std::istream& inp)
{
  std::vector<Polygon_2>              polygons;
  std::vector<Polygon_with_holes_2>   polygons_with_holes;
  std::vector<Polygon_with_holes_2>   join1_res;
  std::vector<Polygon_with_holes_2>   join2_res;
  std::vector<Polygon_with_holes_2>   join3_res;
  std::vector<Polygon_with_holes_2>   intersection1_res;
  std::vector<Polygon_with_holes_2>   intersection2_res;
  std::vector<Polygon_with_holes_2>   intersection3_res;
  std::vector<Polygon_with_holes_2>   symm_diff1_res;
  std::vector<Polygon_with_holes_2>   symm_diff2_res;
  std::vector<Polygon_with_holes_2>   symm_diff3_res;

  read_file(inp,
            polygons,
            polygons_with_holes,
            join1_res,
            join2_res,
            join3_res,
            intersection1_res,
            intersection2_res,
            intersection3_res,
            symm_diff1_res,
            symm_diff2_res,
            symm_diff3_res);

  std::vector<Polygon_with_holes_2>    temp;
  std::back_insert_iterator<std::vector<Polygon_with_holes_2> > oi(temp);

  CGAL::join(polygons.begin(), polygons.end(), oi);
  if(! are_equal(join1_res, temp))
  {
    std::cout<<"join failed...\n";
    return false;
  }
  CGAL::join(polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(join2_res, temp))
  {
    std::cout<<"join failed...\n";
    return false;
  }
  CGAL::join(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(join3_res, temp))
  {
    std::cout<<"join failed...\n";
    return false;
  }

  CGAL::intersection(polygons.begin(), polygons.end(), oi);
  if(! are_equal(intersection1_res, temp))
  {
    std::cout<<"intersection failed...\n";
    return false;
  }
  CGAL::intersection(polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(intersection2_res, temp))
  {
    std::cout<<"intersection failed...\n";
    return false;
  }
  CGAL::intersection(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(intersection3_res, temp))
  {
    std::cout<<"intersection failed...\n";
    return false;
  }

  CGAL::symmetric_difference(polygons.begin(), polygons.end(), oi);
  if(! are_equal(symm_diff1_res, temp))
  {
    std::cout<<"symmetric_difference failed...\n";
    return false;
  }
  CGAL::symmetric_difference(polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(symm_diff2_res, temp))
  {
    std::cout<<"symmetric_difference failed...\n";
    return false;
  }
  CGAL::symmetric_difference(polygons.begin(), polygons.end(),
             polygons_with_holes.begin(), polygons_with_holes.end(), oi);
  if(! are_equal(symm_diff3_res, temp))
  {
    std::cout<<"symmetric_difference failed...\n";
    return false;
  }

  return true;
}

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cerr<<"Missing input file\n";
    std::exit (-1);
  }

  int success = 0;
  for(int i=1; i<argc; ++i)
  {
    std::string str(argv[i]);
    if(str.empty())
      continue;

    std::string::iterator itr = str.end();
    --itr;
    while(itr != str.begin())
    {
      std::string::iterator tmp = itr;
      --tmp;
      if(*itr == 't')
        break;
      
      str.erase(itr);
      itr = tmp;
      
    }
    if(str.size() <= 1)
      continue;
    std::ifstream inp(str.c_str());
    if(!inp.is_open())
    {
      std::cerr<<"Failed to open " <<str<<"\n";
      return (-1);
    }
    if (! test_one_file(inp))
    {
      inp.close();
      std::cout<<str<<": ERROR\n";
      success = -1;
    }
    else
    {
      std::cout<<str<<": succeeded\n";
    }
    inp.close();
  }
  
  return success;
}
