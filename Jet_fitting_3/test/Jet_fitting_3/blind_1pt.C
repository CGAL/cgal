#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <vector>
#include <CGAL/Monge_via_jet_fitting.h>

typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     DPoint;
typedef Data_Kernel::Vector_3     DVector;

typedef double                   LFT;
typedef CGAL::Cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;
typedef My_Monge_via_jet_fitting::Monge_form_condition_numbers My_Monge_form_condition_numbers;
     
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
  My_Monge_form monge_form;
  My_Monge_form_condition_numbers monge_form_condition_numbers;
  //run the main fct
  My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				 d_fitting, d_monge, 
				 monge_form, monge_form_condition_numbers);

  monge_form.comply_wrt_given_normal( -monge_form.n() );
  //OUTPUT on std::cout
  std::cout << monge_form
	    << monge_form_condition_numbers;

  monge_form.dump_4ogl( std::cout, 1 );
  double precision = 0.01;
  assert(monge_form.coefficients()[0] >= -0.2 - precision);
  assert(monge_form.coefficients()[0] <= -0.2 + precision);
  assert(monge_form.coefficients()[1] >= -0.4 - precision);
  assert(monge_form.coefficients()[1] <= -0.4 + precision);
  std::cout << "success\n";
}
