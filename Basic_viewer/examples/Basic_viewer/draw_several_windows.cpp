#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <QMainWindow>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
using PS3=CGAL::Point_set_3<Point>;

int main(void)
{
  /// (1) Some CGAL code that create data structure and fill Graphics_scene.
  std::vector<Pwn> points;

  if(!CGAL::IO::read_points(CGAL::data_file_path("points_3/kitten.xyz"), std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                             .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file " << CGAL::data_file_path("points_3/kitten.xyz") << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron output_mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
    (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

  if (!CGAL::poisson_surface_reconstruction_delaunay
      (points.begin(), points.end(),
       CGAL::First_of_pair_property_map<Pwn>(),
       CGAL::Second_of_pair_property_map<Pwn>(),
       output_mesh, average_spacing))
  { return EXIT_FAILURE; }

  PS3 point_set;
  for(Pwn& it: points)
  { point_set.insert(it.first); }

  CGAL::Graphics_scene<float> graphic_buffer1, graphic_buffer2;
  CGAL::add_in_graphics_scene(point_set, graphic_buffer1);
  CGAL::add_in_graphics_scene(output_mesh, graphic_buffer2);

  /// (2) Qt code that create windows, add them in a layout, and create app.

#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (cgal_test_suite) { return EXIT_SUCCESS; }

  int argc=1;
  const char* argv[2]={"Draw several windows example","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  QMainWindow* mainWindow=new QMainWindow;
  QWidget *centralWidget = new QWidget(mainWindow);
  QHBoxLayout* layout = new QHBoxLayout(mainWindow);

  CGAL::Basic_viewer_qt bv1(mainWindow, graphic_buffer1);
  CGAL::Basic_viewer_qt bv2(mainWindow, graphic_buffer2);
  bv1.set_draw_vertices(true);

  layout->addWidget(&bv1);
  layout->addWidget(&bv2);

  centralWidget->setLayout(layout);
  mainWindow->setCentralWidget(centralWidget);

  mainWindow->show();
  app.exec();

  return EXIT_SUCCESS;
}
