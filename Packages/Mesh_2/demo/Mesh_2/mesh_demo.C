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
#include <CGAL/Mesh_2.h>
#include <CGAL/Mesh_local_size_traits_2.h>
#include <CGAL/Mesh_size_traits_2.h>
#include <CGAL/Mesh_face_base_2.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>

#include "icons.h"

#include "Qt_layer_show_points.h"
#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_triangulation_constraints.h"
#include "Qt_layer_show_circles.h"
#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>


#include <qapplication.h>
#include <qmainwindow.h>
#include <qlayout.h>
#include <qvbox.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qstring.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qpushbutton.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qtimer.h>
#include <qcursor.h>

typedef CGAL::Simple_cartesian<double>  K1;
typedef CGAL::Filtered_kernel<K1>       Kernel;
struct K : public Kernel {};
typedef K::FT                           FT;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Mesh_local_size_traits_2<K> Meshtraits;
typedef CGAL::Constrained_Delaunay_triangulation_2<Meshtraits, Tds,
  CGAL::Exact_predicates_tag> Tr;

typedef K::Point_2 Point;
typedef K::Circle_2 Circle;
typedef CGAL::Polygon_2<K> Polygon;

typedef CGAL::Mesh_2<Tr> Mesh;

template <class M>
class Show_marked_faces : public CGAL::Qt_widget_layer
{
  M &mesh;
  CGAL::Color color;
public:
  Show_marked_faces(Mesh &m, CGAL::Color c=CGAL::GREEN) : mesh(m), color(c) {};

  typedef typename M::Finite_faces_iterator Face_iterator;

  void draw()
  {
    QColor old_fill_color = widget->fillColor();
    int old_line_width = widget->lineWidth();
    *widget << CGAL::FillColor(CGAL::GREEN) << CGAL::LineWidth(0);
    for(Face_iterator fit=mesh.finite_faces_begin();
	fit!=mesh.finite_faces_end();
	++fit)
      if(fit->is_marked())
	*widget << mesh.triangle(fit);
    widget->setFillColor(old_fill_color);
    widget->setLineWidth(old_line_width);
  }
};

class Follow_mouse : public CGAL::Qt_widget_layer
{
  QCursor oldcursor;
public:
  void mouseMoveEvent(QMouseEvent *e)
  {
    FT x=static_cast<FT>(widget->x_real(e->x()));
    FT y=static_cast<FT>(widget->y_real(e->y()));

    widget->new_object(make_object(Point(x, y)));
  }

  void activating()
  {	
    oldcursor = widget->cursor();
    widget->setCursor(crossCursor);
  };
  
  void deactivating()
  {
    widget->setCursor(oldcursor);
  };
};

namespace {

  class Seeds_wrapper_point : public Point {
  public:
    Seeds_wrapper_point(const Point& p): Point(p) {};

    inline Point point() const { return *this; };
  };

  class Seeds : public std::list<Seeds_wrapper_point> {
  public:
    typedef std::list<Seeds_wrapper_point> List;
    typedef List::const_iterator const_iterator;
    typedef const_iterator Vertex_iterator;

    inline const_iterator vertices_begin() { return begin(); }
    inline const_iterator vertices_end() { return end(); }
  };

};

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow() : is_mesh_initialized(false), traits(), mesh(traits)
    {
      QFrame* mainframe = new QFrame(this); 
      QHBoxLayout *hbox = new QHBoxLayout(mainframe);
      hbox->setAutoAdd(true);
      widget = new CGAL::Qt_widget(mainframe, "Main widget");
      widget->setSizePolicy(QSizePolicy( QSizePolicy::Expanding,
					 QSizePolicy::Expanding ));
      setCentralWidget(mainframe);
      resize(700,500);
      mainframe->show();

      QVBox* info = new QVBox(mainframe);
//       QLabel* label = new QLabel("essai", info);
//       new QLineEdit("prout", info);
//       label->setText("Prout");
//       label->show();
      info->show();
      info->adjustSize();

      // STATUSBAR
      statusBar();

      aspect_ratio_label = new QLabel(statusBar(), "aspect_ratio");
      statusBar()->addWidget(aspect_ratio_label, 0, true);

      // LAYERS
      show_points = new CGAL::Qt_layer_show_points<Mesh>(mesh);
      show_seeds = new CGAL::Qt_layer_show_points<Seeds>(seeds,
							 CGAL::BLUE,
							 5,
							 CGAL::CROSS);
      show_triangulation = 
	new CGAL::Qt_layer_show_triangulation<Mesh>(mesh);
      show_marked = 
	new Show_marked_faces<Mesh>(mesh);
      show_constraints = 
	new CGAL::Qt_layer_show_triangulation_constraints<Mesh>(mesh);
      show_circles = 
	new CGAL::Qt_layer_show_circles<Mesh>(mesh, *aspect_ratio_label);
      show_mouse = new CGAL::Qt_layer_mouse_coordinates(*this);

      // layers order, first attached are "under" last attached
      widget->attach(show_marked);
      widget->attach(show_triangulation);
      widget->attach(show_constraints);
      widget->attach(show_circles);
      widget->attach(show_points);
      widget->attach(show_seeds);
      widget->attach(show_mouse);

      show_circles->deactivate();

      get_point = new CGAL::Qt_widget_get_point<K>();
      widget->attach(get_point);
      get_point->deactivate();

      get_polygon = new CGAL::Qt_widget_get_polygon<Polygon>();
      widget->attach(get_polygon);
      get_polygon->deactivate();

      get_seed = new CGAL::Qt_widget_get_point<K>();
      widget->attach(get_seed);
      get_seed->deactivate();

      follow_mouse = new Follow_mouse();
      widget->attach(follow_mouse);
      follow_mouse->deactivate();

      connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
	      this, SLOT(get_cgal_object(CGAL::Object)));


      // TOOLBARS
      // actions: bouding box and mesh
      QToolBar *toolBarActions = new QToolBar("Actions",this);
      QPushButton *pbBounding = 
	new QPushButton("Insert bounding box", toolBarActions);
      connect(pbBounding, SIGNAL(clicked()), this,
	      SLOT(insert_bounding_box()));

      QPushButton *pbMesh = 
	new QPushButton("Mesh", toolBarActions);
      connect(pbMesh, SIGNAL(clicked()), this, SLOT(refineMesh()));

      QPushButton *pbConform = 
	new QPushButton("Conform", toolBarActions);
      connect(pbConform, SIGNAL(clicked()), this, SLOT(conformMesh()));

      QPushButton *pbMeshStep = 
	new QPushButton("Mesh step", toolBarActions);
      connect(pbMeshStep, SIGNAL(clicked()), this,
	      SLOT(refineMeshStep()));

      // TIMER
      timer = new QTimer(this);
      connect(timer, SIGNAL(timeout()),
	      this, SLOT(refineMeshStep()));

      pbMeshTimer = new QPushButton("Auto step", toolBarActions);
      pbMeshTimer->setToggleButton(true);
      connect(pbMeshTimer, SIGNAL(stateChanged(int)),
	      this, SLOT(updateTimer(int)));

      QLineEdit *leTimerInterval = 
	new QLineEdit("1000",toolBarActions);
      timer_interval=1000;
      connect(leTimerInterval, SIGNAL(textChanged(const QString&)),
	      this, SLOT(updateTimerInterval(const QString&)));

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
      QToolButton *pbSeed = 
	new QToolButton(QPixmap( (const char**)marked_xpm ),
			"Seed", "Insert a seed to define a region to mesh", 
			this, SLOT(fake_slot()), 
			toolbarInputs, "seed");

      pbPoint->setToggleButton(true);
      pbPolygon->setToggleButton(true);
      pbSeed->setToggleButton(true);
      bgChooseInputs->insert(pbPoint);
      bgChooseInputs->insert(pbPolygon);
      bgChooseInputs->insert(pbSeed);

      connect(pbPoint, SIGNAL(stateChanged(int)),
	      get_point, SLOT(stateChanged(int)));
      connect(pbPolygon, SIGNAL(stateChanged(int)),
	      get_polygon, SLOT(stateChanged(int)));
      connect(pbSeed, SIGNAL(stateChanged(int)),
	      get_seed, SLOT(stateChanged(int)));

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

      QToolButton *pbShowSeeds
	= new QToolButton(QPixmap( (const char**)seeds_xpm ),
			  "Show seeds", "Display seeds that define the "
			  "region not to mesh",
			  this, SLOT(fake_slot()), 
			  toolbarLayers, "show points");
      pbShowSeeds->setToggleButton(true);
      pbShowSeeds->setOn(true);
      connect(pbShowSeeds, SIGNAL(stateChanged(int)),
	      show_seeds, SLOT(stateChanged(int)));

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
	= new QToolButton(QPixmap( (const char**)contraints_xpm ),
			  "Show constraints", "Display mesh constraints edges",
			  this, SLOT(fake_slot()), 
			  toolbarLayers,
			  "show constraints");
      pbShowConstraints->setToggleButton(true);
      pbShowConstraints->setOn(true);
      connect(pbShowConstraints, SIGNAL(stateChanged(int)),
	      show_constraints, SLOT(stateChanged(int)));

      QToolButton *pbShowMarked 
	= new QToolButton(QPixmap( (const char**)marked_xpm ),
			  "Show marked faces", 
			  "Display faces that will be refined",
			  this, SLOT(fake_slot()), 
			  toolbarLayers,
			  "show marked");
      pbShowMarked->setToggleButton(true);
      pbShowMarked->setOn(true);
      connect(pbShowMarked, SIGNAL(stateChanged(int)),
	      show_marked, SLOT(stateChanged(int)));

      QToolButton *pbShowCircles
	= new QToolButton(QPixmap( (const char**)circle_xpm ),
			  "Show circles", "Display circumcircles of faces",
			  this, SLOT(fake_slot()), 
			  toolbarLayers,
			  "show circles");
      pbShowCircles->setToggleButton(true);
      connect(pbShowCircles, SIGNAL(stateChanged(int)),
	      show_circles, SLOT(stateChanged(int)));

      bgChooseInputs->insert(pbShowCircles);

      // button group trick to connect to widget->redraw() slot
      QButtonGroup *bgLayers = 
	new QButtonGroup("Layers", 0, "layers");
      bgLayers->insert(pbShowPoints);
      bgLayers->insert(pbShowMarked);
      bgLayers->insert(pbShowTriangulation);
      bgLayers->insert(pbShowConstraints);
      bgLayers->insert(pbShowSeeds);
      //      bgLayers->insert(pbShowCircles);
      connect(bgLayers, SIGNAL(clicked(int)),
	      widget, SLOT(redraw()));

      // the standard toolbar
      CGAL::Qt_widget_standard_toolbar *std_toolbar =
	new CGAL::Qt_widget_standard_toolbar(widget, this);
      this->addToolBar(std_toolbar->toolbar(), Top, FALSE);

      setUsesBigPixmaps(true);

      widget->clear_history();

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

      pmCriteria = new QPopupMenu(this);
      pmCriteria->setCheckable(true);
      menuBar()->insertItem("&Criteria", pmCriteria);
      pmCriteria->insertItem("Set &bound...", this, SLOT(setBound()),
			     CTRL+Key_B );
      pmCriteria->insertItem("Set &size bound...", this, SLOT(setSizeBound()),
			     CTRL+Key_S );
      menu_id = pmCriteria->insertItem("Size criteria only under &mouse",
				       this, SLOT(setLocal()),
				       CTRL+Key_L );

      widget->set_window(-500.,500.,-500.,500.);
      widget->setMouseTracking(TRUE);
    };

  // compute bounds of the mesh
  void bounds(FT &xmin, FT &ymin, 
	      FT &xmax, FT &ymax)
    {
      Mesh::Finite_vertices_iterator vi=mesh.finite_vertices_begin();
      xmin=xmax=vi->point().x();
      ymin=ymax=vi->point().y();
      vi++;
      while(vi != mesh.finite_vertices_end())
	{
	  if(vi->point().x() < xmin) xmin=vi->point().x();
	  if(vi->point().x() > xmax) xmax=vi->point().x();
	  if(vi->point().y() < ymin) ymin=vi->point().y();
	  if(vi->point().y() > ymax) ymax=vi->point().y();
	  vi++;
	}
    }

public slots:
  void get_cgal_object(CGAL::Object obj)
    {
      Point p;
      Polygon poly;
      
      if(CGAL::assign(p,obj))
	if(get_seed->is_active())
	  {
	    seeds.push_back(p);
	    mesh.mark_facets(seeds.begin(), seeds.end());
	  }
	else
	  if(follow_mouse->is_active())
	    {
	      typedef Mesh::Face_handle Face_handle;
	      std::list<Face_handle> l;

	      Face_handle fh = mesh.locate(p);
	      traits.set_point(p);
	      if( (fh!=NULL) )
		{
		  const Point&
		    a = fh->vertex(0)->point(),
		    b = fh->vertex(1)->point(),
		    c = fh->vertex(2)->point();
		  
		  if(traits.is_bad_object().operator()(a, b, c))
		    l.push_back(fh);
		}
	      mesh.set_geom_traits(traits, l.begin(), l.end());
	      while( mesh.refine_step() );
	      widget->redraw();
	    }
	  else
	    mesh.insert(p);
      else
	if (CGAL::assign(poly,obj))
	  for(Polygon::Edge_const_iterator it=poly.edges_begin();
	      it!=poly.edges_end();
	      it++)
	    mesh.insert((*it).source(),(*it).target());
      widget->redraw();
    }

//   void redraw_win()
//     {
//       widget->lock();
//       widget->clear();

//       for(Mesh::Edge_iterator it=mesh.edges_begin();
// 	  it!=mesh.edges_end();
// 	  it++)
// 	if(mesh.is_constrained(*it))
// 	  *widget << CGAL::RED << mesh.segment(*it);
// 	else
// 	  *widget << CGAL::BLUE << mesh.segment(*it);
//       *widget << CGAL::GREEN;
//       for(Mesh::Vertex_iterator it=mesh.vertices_begin();
// 	  it!=mesh.vertices_end();
// 	  it++)
// 	*widget << it->point();
//       widget->unlock();
//     }

    //insert a bounding box around the mesh
  void insert_bounding_box()
    {
      FT xmin, xmax, ymin, ymax;
      bounds(xmin, ymin, xmax, ymax);

      FT xcenter=(xmin+xmax)/2,
	ycenter=(ymin+ymax)/2;
      FT xspan = (xmax-xmin)/2,
	yspan = (ymax-ymin)/2;

      Point bb1(xcenter - 1.5*xspan, ycenter - 1.5*yspan);
      Point bb2(xcenter + 1.5*xspan, ycenter - 1.5*yspan);
      Point bb3(xcenter + 1.5*xspan, ycenter + 1.5*yspan);
      Point bb4(xcenter - 1.5*xspan, ycenter + 1.5*yspan);
      mesh.insert(bb1);
      mesh.insert(bb2);
      mesh.insert(bb3);
      mesh.insert(bb4);
      mesh.insert(bb1, bb2);
      mesh.insert(bb2, bb3);
      mesh.insert(bb3, bb4);
      mesh.insert(bb4, bb1);
      widget->redraw();
    }

  void refineMesh()
    {
      saveTriangulationUrgently("last_input.edg");
      mesh.refine(seeds.begin(), seeds.end());
      is_mesh_initialized=true;
      widget->redraw();
    }

  void conformMesh()
    {
      if(!is_mesh_initialized)
	{
	  saveTriangulationUrgently("last_input.edg");
	  mesh.init(seeds.begin(), seeds.end());
	  is_mesh_initialized=true;
	}
      mesh.conform();
      widget->redraw();
    }

  void refineMeshStep()
    {
      if(!is_mesh_initialized)
	{
	  mesh.init(seeds.begin(), seeds.end());
	  is_mesh_initialized=true;
	  saveTriangulationUrgently("last_input.edg");
	}
      if(!mesh.refine_step())
	pbMeshTimer->setOn(false);
      widget->redraw();
    }

  void updateTimer(int i)
    {
      if(i==0)
	timer->stop();
      else
	timer->start(timer_interval);
    }

  void updateTimerInterval(const QString& s)
    {
      bool ok;
      int interval = s.toInt(&ok);
      if(ok)
	{
	  timer_interval=interval;
	  if(timer->isActive())
	    timer->changeInterval(timer_interval);
	}
    }

  void clearMesh()
    {
      mesh.reset();
      seeds.clear();
      is_mesh_initialized=false;
      widget->redraw();
    }

  void openTriangulation()
    {    
      QString s( QFileDialog::getOpenFileName( QString::null,
        "Constrained edges (*.edg);;All files (*)", this ) );
      if ( s.isEmpty() )
        return;
      std::ifstream f(s);
      mesh.read(f);

      // compute bounds
      FT xmin, xmax, ymin, ymax;
      bounds(xmin, ymin, xmax, ymax);
      
      FT xspan = (xmax-xmin)/2,
	yspan = (ymax-ymin)/2;

      widget->set_window(xmin-1.1*xspan, xmax+1.1*xspan, 
			 ymin-1.1*yspan, ymax+1.1*yspan);

      seeds.clear();
      is_mesh_initialized=false;
      widget->redraw();
    }

  void saveTriangulation()
    {
      QString s( QFileDialog::getSaveFileName( "filename.edg", 
        "Constrained edges (*.edg)", this ) );
      if ( s.isEmpty() )
        return;
      std::ofstream of(s);
      mesh.write(of);
    }

  void saveTriangulationUrgently(QString s=QString("dump.edg"))
    {
      std::ofstream of(s);
      mesh.write(of);
    }

  inline
  void fake_slot()
    {
    }

  void setBound() 
    {
      double bound = QInputDialog::getDouble( "Set bound",
					      "Please enter a new bound",
					      mesh.geom_traits().bound());
      traits.set_bound(bound);
      mesh.set_geom_traits(traits);
    }

  void setSizeBound() 
    {
      double bound = QInputDialog::getDouble( "Set bound",
					      "Please enter a new bound",
					      mesh.geom_traits().size_bound());
      traits.set_size_bound(bound);
      mesh.set_geom_traits(traits);
    }

  void setLocal()
    {
      traits.set_local_size(!mesh.geom_traits().is_local_size());
      pmCriteria->setItemChecked(menu_id, traits.is_local_size());
      mesh.set_geom_traits(traits);
      if(traits.is_local_size())
	follow_mouse->activate();
      else
	follow_mouse->deactivate();
    }

private:
  bool is_mesh_initialized;
  Meshtraits traits;
  Mesh mesh;

  Seeds seeds;

  QPopupMenu *pmCriteria;
  int menu_id;

  CGAL::Qt_widget* widget;
  CGAL::Qt_widget_get_point<K>* get_point;
  CGAL::Qt_widget_get_point<K>* get_seed;
  CGAL::Qt_widget_get_polygon<Polygon>* get_polygon;
  Follow_mouse* follow_mouse;

  CGAL::Qt_layer_show_points<Mesh>* show_points;
  CGAL::Qt_layer_show_points<Seeds>* show_seeds;
  CGAL::Qt_layer_show_triangulation<Mesh>* show_triangulation;
  CGAL::Qt_layer_show_triangulation_constraints<Mesh>* show_constraints;
  CGAL::Qt_layer_show_circles<Mesh>* show_circles;
  CGAL::Qt_layer_mouse_coordinates* show_mouse;
  Show_marked_faces<Mesh>* show_marked;

  QLabel* aspect_ratio_label;
  QTimer* timer;
  QPushButton *pbMeshTimer;
  int timer_interval;
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

//   my_previous_failure_function = 
//     CGAL::set_error_handler(cgal_with_exceptions_failure_handler);

  try {
    return app.exec();
  }
  catch(Cgal_exception e) {
    W->saveTriangulationUrgently();
    my_previous_failure_function(e.type, e.expr, e.file, e. line, e.msg);
  }

  return 0;
}

// moc_source_file: mesh_demo.C
#include "mesh_demo.moc"

#endif
