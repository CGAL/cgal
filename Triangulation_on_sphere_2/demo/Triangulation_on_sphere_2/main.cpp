//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/utility.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <QMainWindow>
#include <QFileDialog>
#include <list>
#include <fstream>

#include <qapplication.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

#include "simpleViewer.h"
#include "ui_Mainwindow.h"


typedef CGAL::Projection_sphere_traits_3<Kernel>					Projection_traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Projection_traits>                 Projected_DToS2;


struct Cell_info{
  Kernel::Point_3* pt;
  Cell_info():pt(NULL){}
  ~Cell_info(){if (pt!=NULL) delete pt;}
};

typedef Kernel::Point_3                             Point_3;
Projection_traits traits(Kernel::Point_3(0,0,0));
Projected_DToS2 dtos(traits);

template <class Output_iterator>
void read_points(const char* file_path,Output_iterator out){
  int nb;
  double x,y,z;
  std::ifstream input(file_path);
  if (!input){
    std::cerr << "Error while reading " << file_path << std::endl;
    exit(EXIT_FAILURE);
  }
  
  input >> nb;
  
  for (int i=0;i<nb;++i){
    #warning tmp : must handle the sphere on which are points (parameter of the traits?)
    input >> x >> y >> z;
    *out++=Point_3(x/100.,y/100.,z/100.);
  }
}

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::MainWindow
{
  Q_OBJECT
public:
  MainWindow() {
    setupUi(this);
  }
public slots:
  void open(QString filename) {

    std::vector<Point_3> lst_pt;
    read_points(filename.toUtf8().data(),
                std::back_inserter(lst_pt));
  
    	  
	Projection_traits::Construct_projected_point_3 cst =
	traits.construct_projected_point_3_object();
	dtos.insert(
		boost::make_transform_iterator(lst_pt.begin(), cst),
		boost::make_transform_iterator(lst_pt.end(), cst)
	 );
	  
	 	  

    Point_3 center(0,0,0);
    double scale=1;
  
    MainWindow mainWindow;
    mainWindow.show();
  
    // Instantiate the viewer.
    viewer->open(lst_pt.begin(),lst_pt.end(),dtos,center,scale);
	  
	    }

  void on_action_Quit_triggered() {
    close();
  }

  void on_action_Open_triggered() {
    QString filename = QFileDialog::getOpenFileName(this);
    if(!filename.isNull())
      open(filename);
  }
};

int main(int argc, char** argv)
{
  // Read command lines arguments.
  QApplication application(argc,argv);

  QStringList args = QApplication::arguments();
  args.removeAt(0);

  MainWindow mainWindow;
  if(!args.empty())
  {
    mainWindow.open(args[0]);
  }  
  mainWindow.show();

  // Run main loop.
  return application.exec();
}

#include "Mainwindow.moc"
