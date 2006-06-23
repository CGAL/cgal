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
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;
typedef My_Monge_via_jet_fitting::Monge_form_condition_numbers My_Monge_form_condition_numbers;
       
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
  My_Monge_form monge_form;
  My_Monge_form_condition_numbers monge_form_condition_numbers;
  //run the main fct
  My_Monge_via_jet_fitting do_it(in_points.begin(), in_points.end(),
				 d_fitting, d_monge, 
				 monge_form, monge_form_condition_numbers);

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
	    << "number of points used : " << in_points.size() << std::endl;
  monge_form.dump_verbose(outFile);
  monge_form_condition_numbers.dump_verbose(outFile);
  
  //OUTPUT on std::cout
  CGAL::set_pretty_mode(std::cout);
  std::cout << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl;
  monge_form.dump_verbose(std::cout);
  monge_form_condition_numbers.dump_verbose(std::cout);

  return 1;
}
