#include <CGAL/Cartesian.h>

#ifndef CGAL_USE_LAPACK
int main()
{
  std::cerr << "Skip since LAPACK is not installed" << std::endl;
  std::cerr << std::endl;
  return 0;
}
#else

#include <fstream>
#include <vector>

#include <CGAL/Monge_via_jet_fitting.h>
typedef double                   DFT;
typedef CGAL::Cartesian<DFT>     Data_Kernel;
typedef Data_Kernel::Point_3     DPoint;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form     My_Monge_form;

int main(int argc, char *argv[])
{
  //check command line
  if (argc<4)
    {
      std::cout << " Usage : " << argv[0]
                << " <inputPoints.txt> <d_fitting> <d_monge>" << std::endl;
      return 0;
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
  while (inFile >> x) {
    inFile >> y >> z;
    DPoint p(x,y,z);
    in_points.push_back(p);
  }
  inFile.close();
  // fct parameters
  size_t d_fitting = std::atoi(argv[2]);
  size_t d_monge = std::atoi(argv[3]);

  My_Monge_form monge_form;
  My_Monge_via_jet_fitting monge_fit;
  monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);

  //OUTPUT on std::cout
  CGAL::set_pretty_mode(std::cout);
  std::cout << "vertex : " << in_points[0] << std::endl
	    << "number of points used : " << in_points.size() << std::endl
	    << monge_form;
  std::cout  << "condition_number : " << monge_fit.condition_number() << std::endl
	     << "pca_eigen_vals and associated pca_eigen_vecs :"  << std::endl;
  for (int i=0; i<3; i++)
    std::cout << monge_fit.pca_basis(i).first << std::endl
	      << monge_fit.pca_basis(i).second  << std::endl;
  return 0;
}
#endif //CGAL_USE_LAPACK
