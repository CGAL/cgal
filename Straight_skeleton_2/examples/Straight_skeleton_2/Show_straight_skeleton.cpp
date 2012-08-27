#include<vector>
#include<iterator>
#include<iostream>
#include<iomanip>
#include<string>
#include <fstream>

#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include "dump_to_eps.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                    Point ;
typedef CGAL::Polygon_2<K>            Polygon_2 ;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes ;
typedef CGAL::Straight_skeleton_2<K>  Straight_skeleton ;

typedef boost::shared_ptr<Straight_skeleton> Straight_skeleton_ptr ;

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
      
      Straight_skeleton_ptr ss = CGAL::create_interior_straight_skeleton_2(input);
      if ( ss )
      {
        std::string eps_name ;
        if ( argc > 2  )
             eps_name = argv[2];
        else eps_name = name + ".skeleton.eps" ;
        
        std::ofstream eps(eps_name.c_str()) ;
        if ( eps )  
        {
          std::cerr << "Result: " << eps_name << std::endl ;
          dump_to_eps(input,*ss,eps);
        }
        else
        {
          std::cerr << "Could not open result file: " << eps_name << std::endl ;
        }  
      }
      else
      {
        std::cerr << "ERROR creating interior straight skeleton" << std::endl ;
      }
    }  
    else
    {
      std::cerr << "Could not open input file: " << name << std::endl ;
    }  
  }
  else
  {
    std::cerr << "Computes the straight skeleton in the interior of a polygon with holes and draws it in an EPS file." << std::endl 
              << std::endl 
              << "Usage: show_straight_skeleton <intput_file> [output_eps_file]" << std::endl 
              << std::endl 
              << "       intput_file  Text file describing the input polygon with holes." << std::endl
              << "         (See input_file_format.txt for details)" << std::endl
              << "       output_file     [default='innput_file.skeleton.eps']" << std::endl ; 
  }

  return 0;
}
