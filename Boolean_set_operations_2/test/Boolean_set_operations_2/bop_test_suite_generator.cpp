#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/iterator.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>


// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational        Number_type;

typedef CGAL::Simple_cartesian<Number_type>            Kernel;

typedef CGAL::Polygon_2<Kernel>                        Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>             Polygon_with_holes_2;
typedef CGAL::Gps_segment_traits_2<Kernel>             Traits_2;

bool are_polygons_valid(const std::vector<Polygon_with_holes_2>& vec)
{
  Traits_2 tr;
  unsigned int i=0;
  for(; i < vec.size(); ++i)
  {
    if(!is_valid_polygon_with_holes(vec[i],tr))
      return false;
  }
  return true;
}

template <class T_P1, class T_P2>
std::ostream& write_result_to_file(std::ostream& out, const T_P1& p1,
                                   const T_P2& p2)
{
  Traits_2 tr;
  std::vector<Polygon_with_holes_2> res;
  std::back_insert_iterator<std::vector<Polygon_with_holes_2> > oi(res);

  Polygon_with_holes_2 res_pgn;
  bool intersect = CGAL::join(p1, p2, res_pgn);
  if(intersect)
  {
    out << 1 << std::endl;
    out << res_pgn <<std::endl;
    if(!is_valid_polygon_with_holes(res_pgn,tr))
      std::cout << "warning: invalid polygon was generated" << std::endl;
  }
  else
  {
    out << 0 << std::endl;
  }

  out<<std::endl;
  CGAL::intersection        (p1, p2, oi);
  out<< res.size()<<std::endl;
  std::copy(res.begin(),
            res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(out, "\n"));
  if(!are_polygons_valid(res))
    std::cout << "warning: invalid polygon was generated" << std::endl;
  res.clear();

  out<<std::endl;
  CGAL::difference          (p1, p2, oi);
  out<< res.size()<<std::endl;
  std::copy(res.begin(),
            res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(out, "\n"));
  if(!are_polygons_valid(res))
    std::cout << "warning: invalid polygon was generated" << std::endl;
  res.clear();

  out<<std::endl;
  CGAL::difference          (p2, p1, oi);
  out<< res.size()<<std::endl;
  std::copy(res.begin(),
            res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(out, "\n"));
  if(!are_polygons_valid(res))
    std::cout << "warning: invalid polygon was generated" << std::endl;
  res.clear();

  out<<std::endl;
  CGAL::symmetric_difference(p2, p1, oi);
  out<< res.size()<<std::endl;
  std::copy(res.begin(),
            res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(out, "\n"));
  if(!are_polygons_valid(res))
    std::cout << "warning: invalid polygon was generated" << std::endl;
  res.clear();

  return out;
}

std::ostream& write_complement_to_file(std::ostream& out, const Polygon_2& pgn)
{
  out<<std::endl;
  Polygon_with_holes_2 res_pgn;
  CGAL::complement(pgn, res_pgn);
  //result of a complement of Polygon_2 is one Polygon_with_holes_2
  out << 1 << std::endl;
  out << res_pgn << std::endl;

  Traits_2 tr;
  if(!is_valid_polygon_with_holes(res_pgn,tr))
    std::cout << "warning: invalid polygon was generated" << std::endl;

  return out;
}

std::ostream& write_complement_to_file(std::ostream& out,
                                       const Polygon_with_holes_2& pgn)
{
  out<<std::endl;
  std::vector<Polygon_with_holes_2> res;
  std::back_insert_iterator<std::vector<Polygon_with_holes_2> > oi(res);
  CGAL::complement(pgn, oi);

  out<< res.size()<<std::endl;
  std::copy(res.begin(),
            res.end(),
            std::ostream_iterator<Polygon_with_holes_2>(out, "\n"));
  if(!are_polygons_valid(res))
    std::cout << "warning: invalid polygon was generated" << std::endl;
  res.clear();

  return out;
}

int main(int argc, char *argv[])
{
  if(argc<4)
  {
    std::cerr << "Missing input files!"
              << "When performed as part of the test suite no input file is expected."
              << std::endl;
    return 0;
  }

  std::ifstream inp1 (argv[1]);
  std::ifstream inp2 (argv[2]);
  if(!inp1.is_open() || !inp2.is_open())
  {
    std::cerr << "Failed to open file" << std::endl;
    std::exit(-1);
  }

  std::ifstream test (argv[3]);
  if(test.is_open())
  {
      std::cout<<argv[3] << " already exist, overwrite? (y/n)"<<std::endl;
      char c = std::cin.get();
      if(c != 'y')
        return 0;

  }
  std::ofstream out (argv[3]);
  if(!out.is_open())
  {
    std::cerr<<"Failed to create output file"<<std::endl;
    std::exit(-1);
  }

  Traits_2 tr;
  int type1;
  Polygon_2 p1;
  Polygon_with_holes_2 pwh1;

  inp1 >> type1;
  out << type1 << std::endl;
  if(type1 == 0)
  {
    inp1 >> p1;
    out  << p1;
    if(! is_valid_polygon(p1,tr))
    {
      std::cout << "warning: first input polygon is invalid!!!" << std::endl;
    }
  }
  else
  {
    inp1 >> pwh1;
    out  << pwh1;
    if(! is_valid_polygon_with_holes(pwh1,tr))
    {
      std::cout << "warning: first input polygon is invalid!!!" << std::endl;
    }
  }

  int type2;
  Polygon_2 p2;
  Polygon_with_holes_2 pwh2;

  inp2 >> type2;
  out <<std::endl<< type2 << std::endl;
  if(type2 == 0)
  {
    inp2 >> p2;
    out  << p2;
    if(!is_valid_polygon(p2,tr))
    {
      std::cout << "warning: second input polygon is invalid!!!" << std::endl;
    }
  }
  else
  {
    inp2 >> pwh2;
    out << pwh2;
    if(! is_valid_polygon_with_holes(pwh2,tr))
    {
      std::cout << "warning: second input polygon is invalid!!!" << std::endl;
    }
  }
  out <<std::endl;


  if(type1 == 0 && type2 == 0)
  {
    write_result_to_file(out, p1, p2);
    write_complement_to_file(out, p1);
    write_complement_to_file(out, p2);
    std::cout<<argv[3] <<
      " was generated successfully, dont forget to add it to test_bop.cmd"
             << std::endl;

    return (0);
  }

  if(type1 == 0 && type2 == 1)
  {
    write_result_to_file(out, p1, pwh2);
    write_complement_to_file(out, p1);
    write_complement_to_file(out, pwh2);
    std::cout<<argv[3] <<
      " was generated successfully, dont forget to add it to test_bop.cmd"
             << std::endl;

    return (0);
  }

  if(type1 == 1 && type2 == 0)
  {
    write_result_to_file(out, pwh1, p2);
    write_complement_to_file(out, pwh1);
    write_complement_to_file(out, p2);
    std::cout<<argv[3] <<
      " was generated successfully, dont forget to add it to test_bop.cmd"
             << std::endl;

    return (0);
  }

  if(type1 == 1 && type2 == 1)
  {
    write_result_to_file(out, pwh1, pwh2);
    write_complement_to_file(out, pwh1);
    write_complement_to_file(out, pwh2);
    std::cout<<argv[3] <<
      " was generated successfully, dont forget to add it to test_bop.cmd"
             << std::endl;

    return (0);
  }

  std::cout << "unknown polygon type" << std::endl;

  return (0);
}
