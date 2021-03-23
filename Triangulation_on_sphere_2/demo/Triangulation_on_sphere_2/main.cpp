#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <CGAL/utility.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <QApplication>
#include <QMainWindow>
#include <QFileDialog>

#include <boost/iterator/transform_iterator.hpp>

#include <list>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          Kernel;
typedef Kernel::FT                                                   FT;
typedef Kernel::Point_3                                              Point_3;
typedef Kernel::Segment_3                                            Segment_3;

typedef CGAL::Projection_on_sphere_traits_3<Kernel>                  Projection_traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Projection_traits>  Projected_DToS2;

typedef std::list<std::vector<Point_3> >                             Subsampled_arcs;

#include "Viewer.h"
#include "ui_Mainwindow.h"

template <class Output_iterator>
void read_points(const char* file_path, Output_iterator out)
{
  std::ifstream input(file_path);
  if(!input)
  {
    std::cerr << "Error while reading " << file_path << std::endl;
    std::exit(EXIT_FAILURE);
  }

  double x,y,z;
  while(input >> x >> y >> z)
    *out++ = Point_3(x, y, z);
}

class MainWindow
  : public CGAL::Qt::DemosMainWindow,
    public Ui::MainWindow
{
  Q_OBJECT
public:
  MainWindow()
  {
    setupUi(this);
  }

public slots:
  void open(QString filename)
  {
    std::vector<Point_3> lst_pt;
    read_points(filename.toUtf8().data(), std::back_inserter(lst_pt));

    const Point_3 center(0,0,0);
    const FT radius = 1.0;
    Projection_traits traits(center, radius);
    Projected_DToS2 dtos(lst_pt.begin(), lst_pt.end(), traits);

    std::cout << dtos.number_of_vertices() << " vertices" << std::endl;

    // Instantiate the viewer
    viewer->open(lst_pt.begin(), lst_pt.end(), dtos);
  }

  void on_action_Quit_triggered()
  {
    close();
  }

  void on_action_Open_triggered()
  {
    QString filename = QFileDialog::getOpenFileName(this);
    if(!filename.isNull())
      open(filename);
  }
};

int main(int argc, char** argv)
{
  // Read command lines arguments
  QApplication application(argc,argv);
  MainWindow mainWindow;
  mainWindow.show();

  // Run main loop
  return application.exec();
}

#include "main.moc"
