// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>

int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/basic.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include "Mesh.h"


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>

#include <CGAL/IO/pixmaps/polygon.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/circle.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>

#include "Qt_layer_show_points.h"
#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_triangulation_constraints.h"
#include "Qt_layer_show_circles.h"


#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qstring.h>
//#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qpushbutton.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>

typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       Kernel;
struct K : public Kernel {}; 

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,
  CGAL::Exact_predicates_tag> Tr;

typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon;

typedef CGAL::Mesh<Tr> Mesh;

CGAL::Qt_widget& operator<< (CGAL::Qt_widget& w, const Mesh& mesh)
{
  w.lock();
  mesh.draw_triangulation(w);
  w.unlock();
  return w;
}

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow()
    {
      widget = new CGAL::Qt_widget(this,"Main widget");
      setCentralWidget(widget);
      //      widget->show();
      resize(700,500);
      widget->set_window(0.,1000.,0.,1000.);

      // LAYERS
      show_points = new CGAL::Qt_layer_show_points<Mesh>(mesh);
      show_triangulation = new
	CGAL::Qt_layer_show_triangulation<Mesh>(mesh);
      show_constraints = new
	CGAL::Qt_layer_show_triangulation_constraints<Mesh>(mesh);
      show_circles = new CGAL::Qt_layer_show_circles<Mesh>(mesh);
      widget->attach(show_triangulation);
      widget->attach(show_constraints);
      widget->attach(show_circles);
      widget->attach(show_points);
   
      // TOOLBARS
      get_point = new CGAL::Qt_widget_get_point<K>();
      widget->attach(get_point);
      get_point->deactivate();

      get_polygon = new CGAL::Qt_widget_get_polygon<Polygon>();
      widget->attach(get_polygon);
      get_polygon->deactivate();

      connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
	      this, SLOT(get_cgal_object(CGAL::Object)));
//       connect(widget, SIGNAL(custom_redraw()),
// 	      this, SLOT(redraw_win()));


      // TOOLBARS
      // actions: mesh only
      QToolBar *toolBarActions = new QToolBar("Actions",this);
      QPushButton *pbMesh = 
	new QPushButton("Mesh", toolBarActions);
      connect(pbMesh, SIGNAL(clicked()), this, SLOT(refineMesh()));

      // Inputs: polygons or points
      QToolBar *toolbarInputs = new QToolBar("Inputs",this);
      QButtonGroup *bgChooseInputs = 
	new QButtonGroup("Choose inputs type", 0,
			 "InputType");
      bgChooseInputs->setExclusive(true);
      QToolButton *pbPolygon = 
	new QToolButton(QPixmap( (const char**)polygon_xpm ),
			"Polygon", "Insert polygonal constraints", 
			this, SLOT(fake_slot()), 
			toolbarInputs, "polygon");
      QToolButton *pbPoint = 
	new QToolButton(QPixmap( (const char**)point_xpm ),
			"Point", "Insert points", 
			this, SLOT(fake_slot()), 
			toolbarInputs, "point");
      pbPoint->setToggleButton(true);
      pbPolygon->setToggleButton(true);
      bgChooseInputs->insert(pbPoint);
      bgChooseInputs->insert(pbPolygon);

      connect(pbPoint, SIGNAL(stateChanged(int)),
	      get_point, SLOT(stateChanged(int)));
      connect(pbPolygon, SIGNAL(stateChanged(int)),
	      get_polygon, SLOT(stateChanged(int)));

      pbPolygon->setOn(true);

      // layers: points, edges, constrained edges
      QToolBar *toolbarLayers = new QToolBar("Layers",this);

      QToolButton *pbShowPoints 
	= new QToolButton(QPixmap( (const char**)points_xpm ),
			  "Show points", "Display mesh vertices",
			  this, SLOT(fake_slot()), 
			  toolbarLayers, "show points");
      pbShowPoints->setToggleButton(true);
      pbShowPoints->setOn(true);
      connect(pbShowPoints, SIGNAL(stateChanged(int)),
	      show_points, SLOT(stateChanged(int)));

      QToolButton *pbShowTriangulation 
	= new QToolButton(QPixmap( (const char**)triangulation_xpm ),
			  "Show triangulation", "Display mesh edges", 
			  this, SLOT(fake_slot()), 
			  toolbarLayers, 
			  "show triangulation");
      pbShowTriangulation->setToggleButton(true);
      pbShowTriangulation->setOn(true);
      connect(pbShowTriangulation, SIGNAL(stateChanged(int)),
	      show_triangulation, SLOT(stateChanged(int)));

      QToolButton *pbShowConstraints 
	= new QToolButton(QPixmap( (const char**)line_xpm ),
			  "Show constraints", "Display mesh constraints edges",
			  this, SLOT(fake_slot()), 
			  toolbarLayers,
			  "show constraints");
      pbShowConstraints->setToggleButton(true);
      pbShowConstraints->setOn(true);
      connect(pbShowConstraints, SIGNAL(stateChanged(int)),
	      show_constraints, SLOT(stateChanged(int)));

      QToolButton *pbShowCircles
	= new QToolButton(QPixmap( (const char**)circle_xpm ),
			  "Show circles", "Display circumcircles of faces",
			  this, SLOT(fake_slot()), 
			  toolbarLayers,
			  "show circles");
      pbShowCircles->setToggleButton(true);
      pbShowCircles->setOn(true);
      connect(pbShowCircles, SIGNAL(stateChanged(int)),
	      show_circles, SLOT(stateChanged(int)));
      pbShowCircles->toggle(); // show_circles is desactivated by default

      // button group trick to connect to widget->redraw() slot
      QButtonGroup *bgLayers = 
	new QButtonGroup("Layers", 0, "layers");
      bgLayers->insert(pbShowPoints);
      bgLayers->insert(pbShowTriangulation);
      bgLayers->insert(pbShowConstraints);
      bgLayers->insert(pbShowCircles);
      connect(bgLayers, SIGNAL(clicked(int)),
	      widget, SLOT(redraw()));

      setUsesBigPixmaps(true);

      // the standard toolbar
      CGAL::Qt_widget_standard_toolbar *std_toolbar =
	new CGAL::Qt_widget_standard_toolbar(widget, this);
      this->addToolBar(std_toolbar->toolbar(), Top, FALSE);
      setUsesBigPixmaps(true);
      
      // MENUS
      QPopupMenu *pmMesh = new QPopupMenu(this);
      menuBar()->insertItem("&Mesh", pmMesh);
      pmMesh->insertItem("&Refine mesh", this, SLOT(refineMesh()),
			 CTRL+Key_R );
      pmMesh->insertItem("&Clear mesh", this, SLOT(clearMesh()),
			 CTRL+Key_C );
      pmMesh->insertItem("&Open constrained triangulation...", this,
			 SLOT(openTriangulation()),
			 CTRL+Key_O );
      pmMesh->insertItem("&Save constrained edges...", this,
			 SLOT(saveTriangulation()),
			 CTRL+Key_S );
      pmMesh->insertItem("&Quit", qApp, SLOT(closeAllWindows()),
			 CTRL+Key_Q );

      // STATUSBAR
      statusBar();
    };

public slots:
  void get_cgal_object(CGAL::Object obj)
    {
      Point p;
      Polygon poly;
      
      if(CGAL::assign(p,obj))
	mesh.insert(p);
      else
	if (CGAL::assign(poly,obj))
	  for(Polygon::Edge_const_iterator it=poly.edges_begin();
	      it!=poly.edges_end();
	      it++)
	    mesh.insert((*it).source(),(*it).target());
      widget->redraw();
    }

  void redraw_win()
    {
      widget->lock();
      widget->clear();

      for(Mesh::Edge_iterator it=mesh.edges_begin();
	  it!=mesh.edges_end();
	  it++)
	if(mesh.is_constrained(*it))
	  *widget << CGAL::RED << mesh.segment(*it);
	else
	  *widget << CGAL::BLUE << mesh.segment(*it);
      *widget << CGAL::GREEN;
      for(Mesh::Vertex_iterator it=mesh.vertices_begin();
	  it!=mesh.vertices_end();
	  it++)
	*widget << it->point();
      widget->unlock();
    }
  
  void refineMesh()
    {
      mesh.refine_mesh();
      widget->redraw();
    }
  
  void clearMesh()
    {
      mesh.reset();
      widget->redraw();
    }
  
  void openTriangulation()
    {    
      QString s( QFileDialog::getOpenFileName( QString::null,
        "Constrained edges (*.edg)", this ) );
      if ( s.isEmpty() )
        return;
      ifstream f(s);
      mesh.read(f);
      widget->redraw();
    }

  void saveTriangulation()
    {
      QString s( QFileDialog::getSaveFileName( "filename.edg", 
        "Constrained edges (*.edg)", this ) );
      if ( s.isEmpty() )
        return;
      ofstream of(s);
      mesh.write(of);
    }

  void saveUrgently()
    {
      ofstream of("dump.edg");
      mesh.write(of);
    }

  inline
  void fake_slot()
    {
    }

private:
  Mesh mesh;

  CGAL::Qt_widget* widget;
  CGAL::Qt_widget_get_point<K>* get_point;
  CGAL::Qt_widget_get_polygon<Polygon>* get_polygon;

  CGAL::Qt_layer_show_points<Mesh>* show_points;
  CGAL::Qt_layer_show_triangulation<Mesh>* show_triangulation;
  CGAL::Qt_layer_show_triangulation_constraints<Mesh>* show_constraints;
  CGAL::Qt_layer_show_circles<Mesh>* show_circles;
};


#include <CGAL/assertions.h>
#include <exception>

CGAL::Failure_function my_previous_failure_function;

class Cgal_exception : public std::exception {
public:
  Cgal_exception(const char *t,
		 const char *e,
		 const char* f,
		 int l,
		 const char* m)
    : type(t), expr(e), file(f), line(l), msg(m) {};

  const char *type;
  const char *expr;
  const char* file;
  int line;
  const char* msg;
};

void cgal_with_exceptions_failure_handler(
			const char *type,
			const char *expr,
			const char* file,
			int line,
			const char* msg)
{
  throw Cgal_exception(type,expr,file,line,msg);
}

int main(int argc, char** argv)
{
  QApplication app( argc, argv );
  MyWindow* W = new MyWindow();
  app.setMainWidget(W);
  W->show();

  my_previous_failure_function = 
    CGAL::set_error_handler(cgal_with_exceptions_failure_handler);

  try {
    return app.exec();
  }
  catch(Cgal_exception e) {
    W->saveUrgently();
    my_previous_failure_function(e.type, e.expr, e.file, e. line, e.msg);
  }
}

// moc_source_file: mesh_demo.C
#include "mesh_demo.moc"

#endif
