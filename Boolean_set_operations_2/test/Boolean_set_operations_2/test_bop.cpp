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
//typedef CGAL::Simple_cartesian<Number_type>           Kernel;
// workaround for VC++ 
struct Kernel : public CGAL::Simple_cartesian<Number_type> {};

typedef CGAL::Polygon_2<Kernel>                       Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>            Polygon_with_holes_2;

typedef CGAL::Gps_segment_traits_2<Kernel>            Traits;
typedef CGAL::Polygon_set_2<Kernel>                   Ps;

typedef CGAL::Arr_segment_traits_2<Kernel>            Arr_traits;
typedef CGAL::Gps_traits_2<Arr_traits>                General_traits;
typedef CGAL::General_polygon_set_2<General_traits>   Gps;


void read_file(std::istream& inp,
               bool& intersect,
               Polygon_with_holes_2& join_res,
               std::vector<Polygon_with_holes_2>&  intersection_res,
               std::vector<Polygon_with_holes_2>&  diff1_res,
               std::vector<Polygon_with_holes_2>&  diff2_res,
               std::vector<Polygon_with_holes_2>&  symm_diff_res,
               std::vector<Polygon_with_holes_2>&  comp1_res,
               std::vector<Polygon_with_holes_2>&  comp2_res,
               CGAL::Oriented_side& or_side)
{
  unsigned int x;
  inp >> x;
  intersect = (x != 0);
  if (intersect)
    inp >> join_res;
  
  unsigned int n_pgns;
  inp >> n_pgns;
  intersection_res.resize(n_pgns);
  unsigned int i;
  for (i = 0; i < n_pgns; ++i)
    inp >> intersection_res[i];

  inp >> n_pgns;
  diff1_res.resize(n_pgns);
  for (i = 0; i < n_pgns; ++i)
    inp >> diff1_res[i];

  inp >> n_pgns;
  diff2_res.resize(n_pgns);
  for (i = 0; i < n_pgns; ++i)
    inp >> diff2_res[i];

  inp >> n_pgns;
  symm_diff_res.resize(n_pgns);
  for (i=0; i<n_pgns; ++i)
    inp >> symm_diff_res[i];

  inp >> n_pgns;
  comp1_res.resize(n_pgns);
  for (i = 0; i < n_pgns; ++i)
    inp >> comp1_res[i];

  inp >> n_pgns;
  comp2_res.resize(n_pgns);
  for (i = 0; i < n_pgns; ++i)
    inp >> comp2_res[i];

  int or_side_int;
  inp >> or_side_int;
  or_side = CGAL::sign(or_side_int);
}

template <class Container>
bool are_equal(const Container& l1, Container& l2)
{
  bool eq = ((l1.size() == l2.size()) &&
             (std::equal(l1.begin(), l1.end(), l2.begin())));
  if (!eq) {
    std::ostream_iterator<Polygon_with_holes_2> osi(std::cout, "\n");
    std::cout << "infile: " << std::endl;
    std::copy(l1.begin(), l1.end(), osi);
    std::cout<<std::endl<<"calculated: " <<std::endl;
    std::copy(l2.begin(), l2.end(), osi);
  }
  l2.clear();

  return eq;
}

template <class Traits_>
void my_complement(const typename Traits_::Polygon_2& p,
                   std::vector<typename Traits_::Polygon_with_holes_2>& vec,
                   Traits_& )
{
  vec.resize(1);
  CGAL::complement(p, vec[0]);
}

template <class Traits_>
void my_complement(const typename Traits_::Polygon_with_holes_2& p,
                   std::vector<typename Traits_::Polygon_with_holes_2>& vec,
                   Traits_& )
{
  CGAL::complement(p, std::back_inserter(vec));
}

template <class Polygon1, class Polygon2>
bool test(std::istream& inp, const Polygon1& p1, const Polygon2& p2)
{
  std::vector<Polygon_with_holes_2>  intersection_res_from_file;
  std::vector<Polygon_with_holes_2>  diff1_res_from_file;
  std::vector<Polygon_with_holes_2>  diff2_res_from_file;
  std::vector<Polygon_with_holes_2>  symm_diff_res_from_file;
  std::vector<Polygon_with_holes_2>  comp1_res_from_file;
  std::vector<Polygon_with_holes_2>  comp2_res_from_file;

  Polygon_with_holes_2 join_res_from_file;

  bool intersect;

  CGAL::Oriented_side or_side_from_file;
  
  read_file(inp,
            intersect,
            join_res_from_file,
            intersection_res_from_file,
            diff1_res_from_file,
            diff2_res_from_file,
            symm_diff_res_from_file,
            comp1_res_from_file,
            comp2_res_from_file,
            or_side_from_file);

  std::vector<Polygon_with_holes_2>  temp_result;
  std::back_insert_iterator<std::vector<Polygon_with_holes_2> > oi(temp_result);
  
  CGAL::intersection(p1, p2, oi);
  if (! are_equal(intersection_res_from_file, temp_result))
  {
    std::cout << "intersection 1 failed..." << std::endl;
    return false;
  }

  CGAL::intersection(p2, p1, oi);
  if (! are_equal(intersection_res_from_file, temp_result))
  {
    std::cout << "intersection 2 failed..." << std::endl;
    return false;
  }

  Polygon_with_holes_2 join_res;
  bool do_x;
  
  do_x = CGAL::join(p1, p2, join_res);
  if (do_x != intersect)
  {
    std::cout << "join 11 failed..." << std::endl;
    return false;
  }
  if (do_x)
  {
    if (join_res != join_res_from_file)
    {
      std::cout << "join 12 failed..." << std::endl;
      return false;
    }  
  }

  do_x = CGAL::join(p2, p1, join_res);
  if (do_x != intersect)
  {
    std::cout << "join 21 failed..." << std::endl;
    return false;
  }
  if (do_x)
  {
    if (join_res != join_res_from_file)
    {
      std::cout << "join 22 failed..." << std::endl;
      return false;
    }  
  }

  CGAL::difference(p1 ,p2, oi);
  if (! are_equal(diff1_res_from_file, temp_result))
  {
    std::cout << "diff 1 failed..." << std::endl;
    return false;
  }

  CGAL::difference(p2 ,p1, oi);
  
  if (! are_equal(diff2_res_from_file, temp_result))
  {
    std::cout << "diff 2 failed" << std::endl;
    return false;
  }

  CGAL::symmetric_difference(p1 ,p2, oi);
  if (! are_equal(symm_diff_res_from_file, temp_result))
  {
    std::cout << "symmetric_difference 1 failed" << std::endl;
    return false;
  }

  CGAL::symmetric_difference(p2 ,p1, oi);
  if (! are_equal(symm_diff_res_from_file, temp_result))
  {
    std::cout << "symmetric_difference 2 failed" << std::endl;
    return false;
  }

  Traits tr;
  my_complement(p1, temp_result, tr);

  if (! are_equal(comp1_res_from_file, temp_result))
  {
    std::cout << "complement 1 failed" << std::endl;
    return false;
  }

  my_complement(p2, temp_result, tr);
  if (! are_equal(comp2_res_from_file, temp_result))
  {
    std::cout << "complement 2 failed" << std::endl;
    return false;
  }

  CGAL::Oriented_side or_side = CGAL::oriented_side(p1, p2);
  if (or_side != or_side_from_file) {
    std::cout << "oriented_side 1 failed" << std::endl;
    return false;
  }

  or_side = CGAL::oriented_side(p2, p1);
  if (or_side != or_side_from_file) {
    std::cout << "oriented_side 2 failed" << std::endl;
    return false;
  }
    
  /////////////////////////////////////////

  Ps ps;
  Ps ps1(p1);
  Ps ps2(p2);
 
  ps.intersection(ps1, ps2);
  ps.polygons_with_holes(oi);
  if (!ps.is_valid() || ! are_equal(intersection_res_from_file, temp_result))
  {
    std::cout << "intersection failed" << std::endl;
    return false;
  }

  ps.clear();
  ps.join(ps1, ps2);
  ps.polygons_with_holes(oi);
  if (!ps.is_valid())
  {
    std::cout << "join failed" << std::endl;
    return false;
  }
  if (intersect)
  {
    if (ps.number_of_polygons_with_holes() != 1)
    {
      std::cout << "join failed" << std::endl;
      return false;
    }

    if (temp_result[0] != join_res_from_file)
    {
      std::cout << "join failed" << std::endl;
      return false;
    }
  }
  else
  {
    if (temp_result.size() > 2)
    {
       std::cout << "join failed" << std::endl;
       return false;
    }
    if (temp_result.size() == 1)
    {
      if (! (temp_result[0] == p1 || temp_result[0]==p2))
      {
        std::cout << "join failed" << std::endl;
        return false;
      }
    }
    else
      if (temp_result.size() == 2)
      {
        if (! (temp_result[0]==p1 && temp_result[1]==p2) || 
             (temp_result[0]==p2 && temp_result[1]==p1))
        {
          std::cout << "join failed" << std::endl;
          return false;
        }
      }
  }
  temp_result.clear();

  ps.difference(ps1, ps2);
  ps.polygons_with_holes(oi);
  if (!ps.is_valid() ||! are_equal(diff1_res_from_file, temp_result))
  {
    std::cout << "diff1 failed" << std::endl;
    return false;
  }
  ps.clear();

  ps.difference(ps2, ps1);
  ps.polygons_with_holes(oi);
  if (!ps.is_valid() ||! are_equal(diff2_res_from_file, temp_result))
  {
    std::cout << "diff2 failed" << std::endl;
    return false;
  }
  ps.clear();

  ps.symmetric_difference(ps1, ps2);
  ps.polygons_with_holes(oi);
  if (!ps.is_valid() ||! are_equal(symm_diff_res_from_file, temp_result))
  {
    std::cout << "symmetric_difference failed" << std::endl;
    return false;
  }
  ps.clear();

  return true;
}

bool test_one_file(std::ifstream & inp)
{
  Polygon_2 p1;
  Polygon_with_holes_2 pwh1;

  unsigned int type1;
  inp >> type1;
  if (type1 == 0)
    inp >> p1;
  else
    inp >> pwh1;

  Polygon_2 p2;
  Polygon_with_holes_2 pwh2;

  unsigned int type2;
  inp >> type2;
  if (type2 == 0)
    inp >> p2;
  else
    inp >> pwh2;

  if (type1 == 0 && type2 == 0) return ::test(inp, p1 ,p2);
  if (type1 == 0 && type2 == 1) return ::test(inp, p1 ,pwh2);
  if (type1 == 1 && type2 == 0) return ::test(inp, pwh1, p2);
  if (type1 == 1 && type2 == 1) return ::test(inp, pwh1, pwh2);
  return false;
}


int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << "Missing input file" << std::endl;
    std::exit (-1);
  }

  int success = 0;
  for (int i = 1; i < argc; ++i) {
    std::string str(argv[i]);
    if (str.empty()) continue;

    std::string::iterator itr = str.end();
    --itr;
    while (itr != str.begin()) {
      std::string::iterator tmp = itr;
      --tmp;
      if (*itr == 't') break;
      
      str.erase(itr);
      itr = tmp;
    }
    if (str.size() <= 1) continue;
    std::ifstream inp(str.c_str());
    if (!inp.is_open()) {
      std::cerr << "Failed to open " << str << std::endl;
      return (-1);
    }
    if (! test_one_file(inp)) {
      std::cout << str << ": ERROR" << std::endl;
      ++success;
    }
    else
    {
      std::cout << str << ": succeeded" << std::endl;
    }
    inp.close();
  }
  
  return success;
}
