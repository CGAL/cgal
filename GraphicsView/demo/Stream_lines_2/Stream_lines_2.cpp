#include <fstream>
#include <boost/config.hpp>
#include <boost/version.hpp>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>


// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/StreamLinesGraphicsItem.h>
#include <CGAL/Qt/RegularGridVectorFieldGraphicsItem.h>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif
// the two base classes
#include "ui_Stream_lines_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_grid_2<K> Regular_grid;
typedef CGAL::Runge_kutta_integrator_2<Regular_grid> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator> Stream_lines;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Stream_line_iterator_2 Stream_line_iterator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Point_iterator_2 Point_iterator;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Point_2 Point_2;
typedef CGAL::Stream_lines_2<Regular_grid, Runge_kutta_integrator>::Vector_2 Vector;

typedef K::Iso_rectangle_2 Iso_rectangle_2;



class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Stream_lines_2
{
  Q_OBJECT

private:
  Stream_lines * stream_lines;
  Runge_kutta_integrator * runge_kutta_integrator;
  Regular_grid * regular_grid;
  double density;
  double ratio;
  double integrating;
  int sampling;
  QGraphicsScene scene;

  CGAL::Qt::StreamLinesGraphicsItem<Stream_lines,K> * sli;
  CGAL::Qt::RegularGridVectorFieldGraphicsItem<Regular_grid,K> * rgi;

public:
  MainWindow();

public Q_SLOTS:

  void on_actionLoadPoints_triggered();

  void on_actionClear_triggered();

  void on_actionSavePoints_triggered();

  void on_actionRecenter_triggered();

  virtual void open(QString fileName);

private:
  void generate();

Q_SIGNALS:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), density(12.0), ratio(1.6), integrating(1.0), sampling(1)
{
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);


  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()),
                   this, SLOT(close()));

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&scene);

  // Turn the vertical axis upside down
  this->graphicsView->transform().scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Stream_lines_2.html");
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  connect(this, SIGNAL(openRecentFile(QString)),
          this, SLOT(open(QString)));
}



/*
 *  Qt Automatic Connections
 *  https://doc.qt.io/qt-5/designer-using-a-ui-file.html#automatic-connections
 *
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */


void
MainWindow::on_actionClear_triggered()
{
  Q_EMIT( changed());
}


void
MainWindow::generate()
{
  stream_lines = new Stream_lines(*regular_grid, *runge_kutta_integrator, density, ratio, sampling);

  sli = new CGAL::Qt::StreamLinesGraphicsItem<Stream_lines, K>(stream_lines);
  rgi = new CGAL::Qt::RegularGridVectorFieldGraphicsItem<Regular_grid, K>(regular_grid);

  QObject::connect(this, SIGNAL(changed()),
                   sli, SLOT(modelChanged()));


  rgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  rgi->setEdgesPen(QPen(Qt::gray, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  sli->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(sli);
  scene.addItem(rgi);

  on_actionRecenter_triggered();
  Q_EMIT( changed());
}



void
MainWindow::on_actionLoadPoints_triggered()
{
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#endif
  QString fileName = QFileDialog::getOpenFileName(this,
                                                  tr("Open grid file"),
                                                  "."
                                                #if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
                                                ,tr("WKT files (*.wkt *.WKT)")
                                                #endif
                                                  );
  if(! fileName.isEmpty()){
    open(fileName);
  }
}


void
MainWindow::open(QString fileName)
{
  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::ifstream ifs(qPrintable(fileName));

  runge_kutta_integrator = new Runge_kutta_integrator(integrating);
  double iXSize, iYSize;
  iXSize = iYSize = 512;
  if(fileName.endsWith(".wkt", Qt::CaseInsensitive))
  {
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
    std::vector<std::vector<Point_2> > mp;
    int size= -1;
    do
    {
      std::vector<Point_2> ps;
      CGAL::read_multi_point_WKT(ifs, ps);
      if(size == -1)
        size = static_cast<int>(ps.size());
      else if(ps.size() > 0 && size != static_cast<int>(ps.size()))
        ps.resize(size);
      else if(ps.size() == 0)
        continue;
      mp.push_back(ps);
    }while(ifs.good() && !ifs.eof());
    regular_grid = new Regular_grid(size, static_cast<int>(mp.size()), iXSize, iYSize);
    /*fill the grid with the appropriate values*/
    for (unsigned int i=0;i<static_cast<unsigned int>(size);++i)
      for (unsigned int j=0;j<mp.size();++j)
      {
        regular_grid->set_field(i, j, Vector(mp[j][i].x(), mp[j][i].y()));
      }
#else
    QApplication::restoreOverrideCursor();
    return;
#endif
  }
  else{
    unsigned int x_samples, y_samples;
    ifs >> x_samples;
    ifs >> y_samples;
    regular_grid = new Regular_grid(x_samples, y_samples, iXSize, iYSize);
    /*fill the grid with the appropriate values*/
    for (unsigned int i=0;i<x_samples;i++)
      for (unsigned int j=0;j<y_samples;j++)
      {
        double xval, yval;
        ifs >> xval;
        ifs >> yval;
        regular_grid->set_field(i, j, Vector(xval, yval));
      }
  }
  ifs.close();
  // default cursor
  QApplication::restoreOverrideCursor();
  this->addToRecentFiles(fileName);
  generate();
  Q_EMIT( changed());

}

void
MainWindow::on_actionSavePoints_triggered()
{
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
  QString fileName = QFileDialog::getSaveFileName(this,
                                                  tr("Save points"),
                                                  ".",
                                                  tr("WKT files (*.wkt *.WKT)"));
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));

    std::vector<std::vector<Point_2> >mp;
    mp.resize(regular_grid->get_dimension().second);
    for (int i=0;i<regular_grid->get_dimension().first;++i)
    {
      mp[i].reserve(regular_grid->get_dimension().second);
      for (int j=0;j<regular_grid->get_dimension().second;++j)
      {
        mp[i].push_back(Point_2(regular_grid->get_field(j,i).x(),
                                regular_grid->get_field(j,i).y()));
      }
      CGAL::write_multi_point_WKT(ofs, mp[i]);
    }
    ofs.close();
  }
#endif
}


void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(rgi->boundingRect());
  this->graphicsView->fitInView(rgi->boundingRect(), Qt::KeepAspectRatio);
}


#include <CGAL/Qt/resources.h>

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Stream_lines_2 demo");

  // Import resources from libCGAL (Qt5).
  // See https://doc.qt.io/qt-5/qdir.html#Q_INIT_RESOURCE
  CGAL_QT_INIT_RESOURCES;
  Q_INIT_RESOURCE(Stream_lines_2);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}

#include "Stream_lines_2.moc"
