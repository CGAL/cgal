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

typedef double                   LFT;
typedef CGAL::Cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form  My_Monge_form;
       
int main(int argc, char *argv[])
{
  //check command line  
  if (argc<4)
    {
      std::cout << " Usage : " 
		<< argv[0]
		<< " <inputPoints.txt> <d_fitting> <d_monge>" 
		<< std::endl;
      exit(-1);
    }
  //open the input file
  std::ifstream inFile( argv[1]);
  if ( !inFile ) 
    {
      std::cerr << "cannot open file for input\n";
      exit(-1);
    }
  //initalize the in_points container
  double x, y, z;
  std::vector<DPoint> in_points;
  while (inFile) { 
    inFile >> x >> y >> z;
    DPoint p(x,y,z);
    in_points.push_back(p);
  }
  in_points.pop_back();//the last point is inserted twice.. to fix!
  inFile.close();
    for (size_t j=0;j<in_points.size();j++)
  std::cout << in_points[j]  << std::endl  ; 
  // fct parameters
  size_t d_fitting = std::atoi(argv[2]);
  size_t d_monge = std::atoi(argv[3]);
  My_Monge_form monge_form;
  //run the main fct
 //  My_Monge_via_jet_fitting monge_fit(in_points.begin(), in_points.end(),
// 				 d_fitting, d_monge, 
// 				 monge_form);

   My_Monge_via_jet_fitting monge_fit; 
   //  monge_form = monge_fit()(in_points.begin(), in_points.end(), d_fitting, d_monge);
   monge_form = monge_fit.run(in_points.begin(), in_points.end(), d_fitting, d_monge);

  //OUTPUT on std::cout
  CGAL::set_pretty_mode(std::cout);
  std::cout << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl
	    << monge_form;

  std::cout << "check new access fct"  << std::endl;
  if ( monge_form.coefficients().size() >= 2) 
    std::cout << "d1 : " << monge_form.maximal_principal_direction() << std::endl 
	      << "d2 : " << monge_form.minimal_principal_direction() << std::endl
	      << "k1 : " << monge_form.principal_curvatures(0) << std::endl 
	      << "k2 : " << monge_form.principal_curvatures(1) << std::endl;	      
  if ( monge_form.coefficients().size() >= 6) 
     std::cout << "b0 : " << monge_form.third_order_coefficients(0) << std::endl 
	       << "b1 : " << monge_form.third_order_coefficients(1) << std::endl 
 	       << "b2 : " << monge_form.third_order_coefficients(2) << std::endl 
 	       << "b3 : " << monge_form.third_order_coefficients(3) << std::endl;
  if ( monge_form.coefficients().size() >= 11) 
     std::cout << "c0 : " << monge_form.fourth_order_coefficients(0) << std::endl 
	       << "c1 : " << monge_form.fourth_order_coefficients(1) << std::endl 
 	       << "c2 : " << monge_form.fourth_order_coefficients(2) << std::endl 
 	       << "c3 : " << monge_form.fourth_order_coefficients(3) << std::endl 
 	       << "c4 : " << monge_form.fourth_order_coefficients(4) << std::endl;
  std::cout  << "condition_number : " << monge_fit.condition_number() << std::endl 
	     << "pca_eigen_vals and associated pca_eigen_vecs :"  << std::endl;
  for (int i=0; i<3; i++)
    std::cout << monge_fit.pca_basis(i).first << std::endl 
	      << monge_fit.pca_basis(i).second  << std::endl;
  
  return 1;
}
