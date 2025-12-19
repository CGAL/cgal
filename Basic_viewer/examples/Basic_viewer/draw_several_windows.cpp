#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Basic_viewer.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <QMainWindow>
#endif

#include <vector>
#include <iostream>

using Kernel=CGAL::Exact_predicates_inexact_constructions_kernel;
using Point=Kernel::Point_3;
using Vector=Kernel::Vector_3;
using Pwn=std::pair<Point, Vector>;
using Polyhedron=CGAL::Polyhedron_3<Kernel>;
using PS3=CGAL::Point_set_3<Point>;

int main(void)
{
  /// (1) Some CGAL code that create data structures and fill two Graphics_scene.
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

  CGAL::Graphics_scene scene1, scene2;
  CGAL::add_to_graphics_scene(point_set, scene1);
  CGAL::add_to_graphics_scene(output_mesh, scene2);

  /// (2) Qt code that create windows, add them in a layout, and create app.
#ifdef CGAL_USE_BASIC_VIEWER

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

  CGAL::Qt::Basic_viewer bv1(mainWindow, scene1);
  CGAL::Qt::Basic_viewer bv2(mainWindow, scene2);
  bv1.draw_vertices(true);

  layout->addWidget(&bv1);
  layout->addWidget(&bv2);

  centralWidget->setLayout(layout);
  mainWindow->setCentralWidget(centralWidget);

  mainWindow->show();
  app.exec();
#endif

  return EXIT_SUCCESS;
}
