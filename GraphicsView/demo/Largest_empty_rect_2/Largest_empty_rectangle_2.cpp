#include <fstream>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsRectItem>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/PointsGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_Largest_empty_rectangle_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_iso_rectangle_2;


class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Largest_empty_rectangle_2
{
  Q_OBJECT
  
private:

  Iso_rectangle_2 square;
  Largest_empty_iso_rectangle_2 ler;
  CGAL::Qt::Converter<K> convert;
  std::vector<Point_2> points; 
  QGraphicsScene scene;  

  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * pgi;
  QGraphicsRectItem * rgi;
  QGraphicsLineItem* frame[4];
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  template <typename G>
  void
  on_actionGenerate_triggered()
  {
    Iso_rectangle_2 isor = square;
    Point_2 center = CGAL::midpoint(isor[0], isor[2]);
    Vector_2 offset = center - CGAL::ORIGIN;
    double w = isor.xmax() - isor.xmin();
    double h = isor.ymax() - isor.ymin();
    double radius = (w<h) ? w/2 : h/2;

    G pg(radius);
    bool ok = false;

    const int number_of_points = 
      QInputDialog::getInt(this, 
                               tr("Number of random points"),
                               tr("Enter number of random points"),
                               100,
                               0,
                               (std::numeric_limits<int>::max)(),
                               1,
                               &ok);

    if(!ok) {
      return;
    }

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    points.reserve(points.size() + number_of_points);
    for(int i = 0; i < number_of_points; ++i){
      points.push_back(*pg + offset);
      ler.insert(points.back());
      ++pg;
    }
    
    // default cursor
    QApplication::restoreOverrideCursor();
    Q_EMIT( changed());
  }

public:
  MainWindow();

public Q_SLOTS:

  void on_actionClear_triggered();

  void processInput(CGAL::Object);

  void on_actionRecenter_triggered();
  void on_actionGeneratePointsOnCircle_triggered();
  void on_actionGeneratePointsInSquare_triggered();
  void on_actionGeneratePointsInDisc_triggered();
  void clear();

  void update_largest_empty_rectangle();

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), square(Point_2(-1, -1), Point_2(1,1)), ler(square)
{
  setupUi(this);

  // Add a GraphicItem for the point set
  pgi = new CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> >(&points);

  rgi = new QGraphicsRectItem(convert(square));

  Point_2 bl(-1,-1), br(1,-1), tl(-1,1), tr(1,1);
  
  frame[0] = new QGraphicsLineItem(convert(Segment_2(bl, br)));
  frame[1] = new QGraphicsLineItem(convert(Segment_2(br, tr)));
  frame[2] = new QGraphicsLineItem(convert(Segment_2(tr, tl)));
  frame[3] = new QGraphicsLineItem(convert(Segment_2(tl, bl)));

  QObject::connect(this, SIGNAL(changed()),
		   pgi, SLOT(modelChanged()));

  QObject::connect(this, SIGNAL(changed()),
		   this, SLOT(update_largest_empty_rectangle()));

  pgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  rgi->setBrush(QBrush(Qt::cyan));
  scene.addItem(pgi);
  scene.addItem(rgi);
  scene.addItem(frame[0]);
  scene.addItem(frame[1]);
  scene.addItem(frame[2]);
  scene.addItem(frame[3]);

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

 
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 1, false); // inputs a list with one point
  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
   
  scene.installEventFilter(pi);
  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-1, -1, 1, 1);
  this->graphicsView->setScene(&scene);

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Largest_empty_rectangle_2.html");
  this->addAboutCGAL();


}


/* 
 *  Qt Automatic Connections
 *  http://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */



void
MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> input;
  if(CGAL::assign(input, o)){
    if(input.size() == 1) {
      Point_2 p = input.front();
      if(! square.has_on_unbounded_side(p)){
        points.push_back(p);
        ler.insert(p);
      }
    }
    Q_EMIT( changed());
  }

}


void
MainWindow::update_largest_empty_rectangle()
{
  rgi->setRect(convert(ler.get_largest_empty_iso_rectangle()));
}


void
MainWindow::on_actionClear_triggered()
{
  clear();
  Q_EMIT( changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(convert(square));
  this->graphicsView->fitInView(convert(square), Qt::KeepAspectRatio);  
}

void
MainWindow::on_actionGeneratePointsOnCircle_triggered()
{
  typedef CGAL::Random_points_on_circle_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInSquare_triggered()
{
  typedef CGAL::Random_points_in_square_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInDisc_triggered()
{
  typedef CGAL::Random_points_in_disc_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}



void
MainWindow::clear()
{
  points.clear();
  ler.clear();
  rgi->setRect(convert(square));
}


#include "Largest_empty_rectangle_2.moc"
#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Largest_empty_rectangle_2 demo");

  // Import resources from libCGAL (Qt5).
  // See http://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Largest_empty_rectangle_2);

  MainWindow mainWindow;
  mainWindow.show();
  mainWindow.on_actionRecenter_triggered();
  return app.exec();
}
