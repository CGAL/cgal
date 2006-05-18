#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <vector>
#include "../../include/CGAL/Monge_via_jet_fitting.h" 
#include "LinAlg_lapack.h" 
 
typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     DPoint;
typedef CGAL::Monge_rep<Data_Kernel>   My_Monge_rep;

typedef double                   LFT;
typedef CGAL::Cartesian<LFT>     Local_Kernel;
typedef CGAL::Monge_info<Local_Kernel> My_Monge_info;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel, Local_Kernel, Lapack> My_Monge_via_jet_fitting;

       
int main(int argc, char *argv[])
{
  //check command line  
  if (argc<5)
    {
      std::cout << "Usage : blind_1pt <inputPoints.txt> <output.txt> <d_fitting>, <d_monge>" 
		<< std::endl;
      exit(-1);
    }
  //open the input file
  char name_in[20];
  sprintf(name_in, "%s", argv[1]);
  std::cout << name_in << '\n';

  std::ifstream inFile( name_in, std::ios::in);
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
  int d_fitting = std::atoi(argv[3]);
  int d_monge = std::atoi(argv[4]);
  My_Monge_rep monge_rep;
  My_Monge_info monge_info;
  //run the main fct
  My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				 d_fitting, d_monge, 
				 monge_rep, monge_info);

  //open a file for output
  char name_out[20];
  sprintf(name_out, "%s", argv[2]);
  std::cout << name_out << '\n';

  std::ofstream outFile( name_out, std::ios::out);
  if ( !outFile ) 
    {
      std::cerr << "cannot open file for output\n";
      exit(-1);
    }

  //OUTPUT on outFile
  CGAL::set_pretty_mode(outFile);
  outFile   << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl
	    << "origin : " << monge_rep.origin_pt() << std::endl
	    << "d1 : " << monge_rep.d1() << std::endl 
	    << "d2 : " << monge_rep.d2() << std::endl
	    << "n : " << monge_rep.n() << std::endl
	    << "cond_nb : " << monge_info.cond_nb() << std::endl 
	    << "pca_eigen_vals " << monge_info.pca_eigen_vals()[0] 
	    << " " << monge_info.pca_eigen_vals()[2] 
	    << " " << monge_info.pca_eigen_vals()[3]  << std::endl 
	    << "pca_eigen_vecs : " << std::endl 
	    << monge_info.pca_eigen_vecs()[0] << std::endl 
	    << monge_info.pca_eigen_vecs()[1] << std::endl 
	    << monge_info.pca_eigen_vecs()[2] << std::endl
	    << std::endl ;
  if ( d_monge >= 2) 
    outFile   << "k1 : " << monge_rep.coefficients()[0] << std::endl 
	      << "k2 : " << monge_rep.coefficients()[1] << std::endl;
  if ( d_monge >= 3) 
    outFile   << "b0 : " << monge_rep.coefficients()[2] << std::endl 
	      << "b1 : " << monge_rep.coefficients()[3] << std::endl
	      << "b2 : " << monge_rep.coefficients()[4] << std::endl 
	      << "b3 : " << monge_rep.coefficients()[5] << std::endl; 
  if ( d_monge >= 4) 
    outFile   << "c0 : " << monge_rep.coefficients()[6] << std::endl 
	      << "c1 : " << monge_rep.coefficients()[7] << std::endl
	      << "c2 : " << monge_rep.coefficients()[8] << std::endl 
	      << "c3 : " << monge_rep.coefficients()[9] << std::endl 
	      << "c4 : " << monge_rep.coefficients()[10] << std::endl; 
  
   //OUTPUT on std::cout
  CGAL::set_pretty_mode(std::cout);
  std::cout << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl
	    << "origin : " << monge_rep.origin_pt() << std::endl
	    << "d1 : " << monge_rep.d1() << std::endl 
	    << "d2 : " << monge_rep.d2() << std::endl
	    << "n : " << monge_rep.n() << std::endl
	    << "cond_nb : " << monge_info.cond_nb() << std::endl 
	    << "pca_eigen_vals " << monge_info.pca_eigen_vals()[0] 
	    << " " << monge_info.pca_eigen_vals()[2] 
	    << " " << monge_info.pca_eigen_vals()[3]  << std::endl 
	    << "pca_eigen_vecs : " << std::endl 
	    << monge_info.pca_eigen_vecs()[0] << std::endl 
	    << monge_info.pca_eigen_vecs()[1] << std::endl 
	    << monge_info.pca_eigen_vecs()[2] << std::endl
	    << std::endl ;
  if ( d_monge >= 2) 
    std::cout << "k1 : " << monge_rep.coefficients()[0] << std::endl 
	      << "k2 : " << monge_rep.coefficients()[1] << std::endl;
  if ( d_monge >= 3) 
    std::cout << "b0 : " << monge_rep.coefficients()[2] << std::endl 
	      << "b1 : " << monge_rep.coefficients()[3] << std::endl
	      << "b2 : " << monge_rep.coefficients()[4] << std::endl 
	      << "b3 : " << monge_rep.coefficients()[5] << std::endl; 
  if ( d_monge >= 4) 
    std::cout << "c0 : " << monge_rep.coefficients()[6] << std::endl 
	      << "c1 : " << monge_rep.coefficients()[7] << std::endl
	      << "c2 : " << monge_rep.coefficients()[8] << std::endl 
	      << "c3 : " << monge_rep.coefficients()[9] << std::endl 
	      << "c4 : " << monge_rep.coefficients()[10] << std::endl; 
  return 1;
}
