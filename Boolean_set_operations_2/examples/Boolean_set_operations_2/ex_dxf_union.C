
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Dxf_bsop_reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>             NT;
typedef CGAL::Simple_cartesian<NT>                  Kernel;
 
typedef CGAL::Gps_circle_segment_traits_2<Kernel>     Traits_2;
typedef Traits_2::Polygon_2                           Circ_polygon;
typedef Traits_2::Polygon_with_holes_2                Circ_polygon_with_holes;
typedef std::vector<Circ_polygon>                     Circ_pgn_vec;
typedef std::vector<Circ_polygon_with_holes>          Circ_pgn_with_holes_vec;

typedef CGAL::General_polygon_set_2<Traits_2>         Gps;

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cout<<"Missing DXF file"<<std::endl;
    return 0;
  }

  std::ifstream input_file (argv[1]);
  if(!input_file.is_open())
  {
    std::cout<<"Failed to open the file"<<std::endl;
    return 0;
  }

  bool simplify = true;
  if(argc >= 3)
  {
    int i = atoi(argv[2]);
    if(i == 0)
      simplify = false;
  }

  Gps gps;

  Circ_pgn_vec pgns;
  Circ_pgn_with_holes_vec pgns_with_holes;
  CGAL::Dxf_bsop_reader<Kernel>  reader;
  reader(input_file,
         std::back_inserter(pgns),
         std::back_inserter(pgns_with_holes),
         simplify);

  CGAL::Timer t;
  
  std::cout<<"Performing union\n";
  
  t.start();
  gps.join(pgns.begin(), pgns.end(),
           pgns_with_holes.begin(), pgns_with_holes.end());
  t.stop();
  
  std::cout<<"Union time : "<< t.time()<<" seconds\n";

  std::cout<<"|V| = " << gps.arrangement().number_of_vertices()<<"\n";
  std::cout<<"|E| = " << gps.arrangement().number_of_edges()<<"\n";
  std::cout<<"|F| = " << gps.arrangement().number_of_faces()<<"\n";

  input_file.close();  
  return 0;
}
