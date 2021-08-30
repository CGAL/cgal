#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include "dump_to_eps.h"

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;

typedef boost::shared_ptr<Polygon_with_holes> Polygon_with_holes_ptr ;

typedef std::vector<Polygon_with_holes_ptr> Polygon_with_holes_ptr_vector ;

int main( int argc, char* argv[] )
{
  Polygon_with_holes input ;

  if ( argc > 1 )
  {
    std::string name = argv[1] ;

    std::cout << "Input file: " << name << std::endl ;

    std::ifstream is(name.c_str()) ;
    if ( is )
    {
      is >> input ;

      assert(input.outer_boundary().is_counterclockwise_oriented());
      for(Polygon_with_holes::Hole_const_iterator it = input.holes_begin();
          it != input.holes_end();
          ++it){
        assert(it->is_counterclockwise_oriented());
      }

      double lOffset = 0.25 ;

      if ( argc > 2 )
        lOffset = std::atof(argv[2]);

      std::cout << "Offsetting at: " << lOffset << std::endl ;

      Polygon_with_holes_ptr_vector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(lOffset,input);

      std::string eps_name ;
      if ( argc > 3  )
           eps_name = argv[3];
      else eps_name = name + ".offset.eps" ;

      std::ofstream eps(eps_name.c_str()) ;
      if ( eps )
      {
        std::cerr << "Result: " << eps_name << std::endl ;
        dump_to_eps(input,offset_polygons,eps);
      }
      else
      {
        std::cerr << "Could not open result file: " << eps_name << std::endl ;
      }
    }
    else
    {
      std::cerr << "Could not open input file: " << name << std::endl ;
    }
  }
  else
  {
    std::cerr << "Computes the interior offset of a polygon with holes and draws the result in an EPS file." << std::endl
              << std::endl
              << "Usage: show_offset_polygon <intput_file> [offset_distance] [output_eps_file]" << std::endl
              << std::endl
              << "       intput_file  Text file describing the input polygon with holes." << std::endl
              << "         (See inputfile_format.txt for details)" << std::endl
              << "       offset_distance [default=0.25]." << std::endl
              << "       output_file     [default='innput_file.offset.eps']" << std::endl ;
  }

  return 0;
}
