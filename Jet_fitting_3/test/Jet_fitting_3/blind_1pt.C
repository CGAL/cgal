#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <vector>
#include "../../include/CGAL/Monge_via_jet_fitting.h" 
#include "include/LinAlg_lapack.h" 

typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     DPoint;
typedef Data_Kernel::Vector_3     DVector;
typedef CGAL::Monge_rep<Data_Kernel>   My_Monge_rep;

typedef double                   LFT;
typedef CGAL::Cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_info<Local_Kernel> My_Monge_info;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel, Local_Kernel, Lapack> My_Monge_via_jet_fitting;

       
int main()
{
  //open the input file
  std::ifstream inFile( "data/in_points_d4.txt", std::ios::in);
  if ( !inFile ) 
    {
      std::cerr << "cannot open file for input\n";
      exit(-1);
    }
  //initalize the in_points container
  double x, y, z;
  std::vector<DPoint> in_points;
  char ch[40];
  while (inFile >> ch) {
    x = atof(ch);
    inFile >> ch;
    y = atof(ch); 
    inFile >> ch;
    z = atof(ch);
    DPoint p(x,y,z);
    in_points.push_back(p);
  }
  inFile.close();

  // fct parameters
  int d_fitting = 4;
  int d_monge = 4;
  My_Monge_rep monge_rep;
  My_Monge_info monge_info;
  //run the main fct
  My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				 d_fitting, d_monge, 
				 monge_rep, monge_info);

  monge_rep.comply_wrt_given_normal( -monge_rep.n() );
  //OUTPUT on std::cout
  monge_rep.dump_verbose( std::cout );
  monge_rep.dump_4ogl( std::cout, 1 );
  monge_info.dump_verbose( std::cout );
  std::cout << "success\n";
}
