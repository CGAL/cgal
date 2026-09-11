#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/shared_ptr.hpp>
#include <CGAL/Reconstruction_from_parallel_slices_3/check_and_fix_input.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

std::ostream& write_polygon(std::vector<Point_3>& points, std::ostream& out){
  out << std::setprecision(17) <<  points.size() << "\n";
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point_3>(out,"\n"));
  return out;
}

bool read_polygon(std::vector<Point_3>& points, std::istream& input){
  std::size_t nbpt;
  Kernel::Point_3 p;
  if (input >> nbpt && input)
    input >> p;
  else
    return false;
  points.resize(nbpt);
  points[0]=p;

  for (std::size_t i=1;i<nbpt;++i)
    input >> points[i];
  return true;
}

int main(int argc,char** argv){
  if (argc!=3)
  {
    std::cerr << "Please provide a file of polyline (*.cgal) and the constant coordinate index (0,1,2)\n";
    return 1;
  }

  std::ifstream input(argv[1]);
  int cc=atoi(argv[2]);

  typedef std::vector< Point_3 > Contour;

  namespace RFPS = CGAL::Reconstruction_from_parallel_slices;

  RFPS::Contour_checker_and_fixer<Kernel,true> checker(cc);
  RFPS::Error_code errcode;
  double last_elev;
  int k=0;
  std::vector<int> slice;

  std::ofstream ok("ok.cgal");
  std::ofstream bad("bad.cgal");

  do{
    boost::shared_ptr<Contour> ctr_ptr( new Contour() );
    if ( read_polygon(*ctr_ptr,input) )
    {
      double elev = ctr_ptr->front()[cc];
      if ( checker.empty() )
      {
        last_elev=elev;
        checker.begin_slice();
      }
      else{
        if (last_elev!=elev)
        {
          // a slice is completed
          if ( !checker.end_slice() )
          {
            std::cerr << "Popping the slice made of polygons {";
            std::copy(slice.begin(), slice.end(), std::ostream_iterator<int>(std::cerr," ") );
            std::cerr << " } that have intersecting polygons (or is empty)\n";
            for (std::size_t i=0; i<checker.size_of_contours(); ++i)
              write_polygon(*checker.contour(i), bad);
            checker.pop_slice_back();
          }
          slice.clear();
          checker.begin_slice();
          last_elev=elev;
        }
      }

      errcode=checker.add_contour_to_slice(ctr_ptr);
      switch (errcode)
      {
        case RFPS::INCONSISTENT_CST_COORD:
          std::cerr << "Polygon " << k << " has inconsistent constant coordinate\n";
          write_polygon(*ctr_ptr, bad);
        break;
        case RFPS::DEGENERATE_POLYGON:
          std::cerr << "Polygon " << k << " is degenerate\n";
          write_polygon(*ctr_ptr, bad);
        break;
        case RFPS::POLYGON_NOT_SIMPLE:
          std::cerr << "Polygon " << k << " is not simple\n";
          write_polygon(*ctr_ptr, bad);
          checker.pop_contour_back();
        break;
        case RFPS::VALID_OR_TRIVIALLY_FIXED_POLYGON:
          slice.push_back(k);
          write_polygon(*ctr_ptr, ok);
        break;
        case RFPS::LOGIC_ERROR:
          std::cerr << "Logic error!!!!" << std::endl;
      };
    }
    else
    {
      if ( !checker.empty() && !checker.end_slice() )
      {
        std::cerr << "Popping the slice made of polygons {";
        std::copy(slice.begin(), slice.end(), std::ostream_iterator<int>(std::cerr," ") );
        std::cerr << "} that have intersecting polygons (or is empty)\n";
        checker.pop_slice_back();
      }
      slice.clear();
      break;
    }
    ++k;
  }
  while (true);

  ok.close();
  bad.close();

  if (!checker.empty())
  {
    std::ofstream slices_ok("slices_ok.cgal");
    checker.output_slices(slices_ok);
  }
}
