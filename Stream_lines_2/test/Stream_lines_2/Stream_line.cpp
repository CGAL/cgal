#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Euler_integrator_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>
#include <CGAL/Triangular_field_2.h>

#include <fstream>
#include <sstream>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> CK;
typedef CGAL::Filtered_kernel<CK> K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Euler_integrator_2<Field> Euler_integrator;
typedef CGAL::Stream_lines_2<Field, Euler_integrator> Stl_euler;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl_rungk;
typedef Stl_euler::Stream_line_iterator_2 Stl_iterator;
typedef Stl_euler::Point_iterator_2 Point_iterator;
typedef Stl_euler::Point_2 Point_2;
typedef Stl_euler::Vector_2 Vector_2;
typedef Stl_rungk::Stream_line_iterator_2 Stl_iterator;
typedef Stl_rungk::Point_iterator_2 Point_iterator;
typedef Stl_rungk::Point_2 Point_2;
typedef Stl_rungk::Vector_2 Vector_2;

typedef CGAL::Triangular_field_2<K>                                        TRField;
typedef CGAL::Runge_kutta_integrator_2<TRField>                            TRRunge_kutta_integrator;
typedef CGAL::Stream_lines_2<TRField, TRRunge_kutta_integrator>            TRStl;
typedef TRStl::Stream_line_iterator_2                                      TRStl_iterator;


int main()
{
  for (int i=1;i<=23;i++)
  {
    std::ostringstream os;
    os << "data/" << i << ".vec.cin";
    std::string file_name = os.str();

    std::ifstream infile(file_name.c_str(), std::ios::in);
    if (!infile) {
      std::cout << "Unable to open file " << os.str()
                << "  ... skipping it." << std::endl;
      continue;
    }

    double iXSize, iYSize;
    unsigned int x_samples, y_samples;
    iXSize = iYSize = 512;
    infile >> x_samples >> y_samples;
    Field regular_grid_2(x_samples, y_samples, iXSize, iYSize);
    /*fill the grid with the appropriate values*/
    for (unsigned int i=0;i<x_samples;i++)
      for (unsigned int j=0;j<y_samples;j++)
    {
      double xval, yval;
      infile >> xval >> yval;
      regular_grid_2.set_field(i, j, Vector_2(xval, yval));
    }
    if (!infile) {
      std::cout << "An error occurred while parsing file "
                << os.str() << std::endl;
    }
    infile.close();
    std::cout << "processing... (" << file_name << ")\n";
    /* the placement of streamlines using Runge-Kutta integrator */
    std::cout << "using Runge-Kutta integrator\n";
    Runge_kutta_integrator runge_kutta_integrator;
    double dSep = 3.5;
    double dRat = 1.6;
    Stl_rungk Stream_lines_rungk(regular_grid_2, runge_kutta_integrator,dSep,dRat,1);
    std::cout << "placement generated\n";
    std::cout << "updating parameters...\n";
    Stream_lines_rungk.set_separating_distance(12.0);
    Stream_lines_rungk.set_saturation_ratio(1.7);
    Stream_lines_rungk.update();
    std::cout << "placement generated\n";
    std::ofstream fw("result.stl",std::ios::out);
    Stream_lines_rungk.print_stream_lines(fw);
    fw.close();
    /* the placement of streamlines using Euler integrator */
    std::cout << "using Euler integrator\n";
    Euler_integrator euler_integrator;
    dSep = 3.5;
    dRat = 1.6;
    Stl_euler Stream_lines(regular_grid_2, euler_integrator,dSep,dRat,1);
    std::cout << "placement generated\n";
  }
  
  std::cout << "placement of streamlines on irregular field\n";
  TRRunge_kutta_integrator trrunge_kutta_integrator(1);
  std::ifstream inp("data/datap.tri.cin");
  std::ifstream inv("data/datav.tri.cin");
  std::istream_iterator<Point_2> beginp(inp);
  std::istream_iterator<Vector_2> beginv(inv);
  std::istream_iterator<Point_2> endp;
  TRField triangular_field(beginp, endp, beginv);
  /* the placement of streamlines */
  std::cout << "processing...\n";
  double dSep = 38.0;
  double dRat = 1.6;
  TRStl TRStream_lines(triangular_field, trrunge_kutta_integrator,dSep,dRat);
  std::cout << "placement generated\n";
  std::cout << "updating parameters...\n";
  TRStream_lines.set_separating_distance(45.0);
  TRStream_lines.set_saturation_ratio(1.7);
  TRStream_lines.update();
  std::cout << "placement generated\n";
  
  std::cout << "success\n";
}
