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
#include <climits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#define CGAL_MESH_2_USE_TIMERS
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_local_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>

#include "icons.h"

#include "Qt_layer_show_points.h"
#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_triangulation_constraints.h"
#include "Qt_layer_show_circles.h"
#include "Show_clusters.h"
#include <CGAL/IO/Qt_widget_show_mouse_coordinates.h>


#include <qapplication.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qlayout.h>
#include <qvbox.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qstring.h>
#include <qregexp.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qpushbutton.h>
#include <qbuttongroup.h>

#include <qlabel.h>
#include <qlineedit.h>
#include <qvalidator.h>
#include <qtimer.h>
#include <qcursor.h>
#include <qslider.h>
#include <qspinbox.h>
#include <qcheckbox.h>

#ifdef TESTING
#  include <CGAL/MP_Float.h>
   typedef CGAL::Simple_cartesian<CGAL::MP_Float> Kernel;
#else // !TESTING
#  include <CGAL/Filtered_kernel.h>
   typedef CGAL::Simple_cartesian<double>  K1;
   typedef CGAL::Filtered_kernel<K1>       Kernel;
#endif // #ifdef TESTING

struct K : public Kernel {};
typedef K::FT                           FT;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,
  CGAL::Exact_predicates_tag> Tr;
typedef CGAL::Delaunay_mesh_local_size_criteria_2<Tr> Criteria;

typedef K::Point_2 Point_2;
typedef K::Circle_2 Circle;
typedef CGAL::Polygon_2<K> CGALPolygon;

typedef CGAL::Delaunay_mesher_2<Tr, Criteria> Mesher;
typedef Tr::Vertex Vertex;

template <class CDT>
void
read_constraints(CDT& t, std::istream &f)
{
  typedef typename CDT::Point Point;

  t.clear();

  int nedges = 0;
  f>>nedges;

  for(int n = 0; n<nedges; n++) {
    Point p1, p2;
    f >> p1 >> p2;
    t.insert_constraint(p1, p2);
  }
}

template <class CDT>
void
write_constraints(const CDT& t, std::ostream &f)
{
  typedef typename CDT::Finite_edges_iterator Finite_edges_iterator;

  int number_of_constrained_edges = 0;
  for(Finite_edges_iterator it = t.finite_edges_begin();
      it != t.finite_edges_end();
      ++it)
    if((*it).first->is_constrained((*it).second))
      ++number_of_constrained_edges;

  f << number_of_constrained_edges << std::endl;
  for(Finite_edges_iterator eit = t.finite_edges_begin();
      eit!=t.finite_edges_end();
      ++eit)
    if((*eit).first->is_constrained((*eit).second))
      {
        f << (*eit).first->vertex(t.cw((*eit).second))->point() << " "
          << (*eit).first->vertex(t.ccw((*eit).second))->point() <<std::endl;
      }
}

template <class Tr>
class Show_marked_faces : public CGAL::Qt_widget_layer
{
  Tr *&tr;
  CGAL::Color color;
public:
  Show_marked_faces(Tr *&t, CGAL::Color c=CGAL::GREEN) : tr(t),
    color(c) {};

  typedef typename Tr::Finite_faces_iterator Face_iterator;

  void draw()
  {
    QColor old_fill_color = widget->fillColor();
    int old_line_width = widget->lineWidth();
    *widget << CGAL::FillColor(CGAL::GREEN) << CGAL::LineWidth(0);
    for(Face_iterator fit=tr.finite_faces_begin();
        fit!=tr.finite_faces_end();
        ++fit)
      if(fit->is_marked())
        *widget << tr.triangle(fit);
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

    widget->new_object(CGAL::make_object(Point_2(x, y)));
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

struct Vertex_to_point {
  const Point_2& operator()(const Vertex& v) const
    {  return v.point(); }
};

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow() : criteria()
    {
      mesher = new Mesher(criteria);

      // --- DEMO WINDOW'S LAYOUT ---
      // a main frame, with a QHBoxLayout (border=0, space=0)
      QFrame* mainframe = new QFrame(this, "mainframe");
      QHBoxLayout *mainlayout = new QHBoxLayout(mainframe, 0, 0,
                                                "mainlayout");

      // a Qt_widget at the left of the main frame
      widget = new CGAL::Qt_widget(mainframe, "Main widget");
      widget->setSizePolicy(QSizePolicy( QSizePolicy::Expanding,
                                         QSizePolicy::Expanding ));
      mainlayout->addWidget(widget);

      // an other frame "infoframe" at the right of the main frame
      // with a QVBoxLayout (border=10, space=5)
      QFrame* infoframe = new QFrame(mainframe, "infoframe");
      infoframe->setFrameStyle( QFrame::Box | QFrame::Raised );
      infoframe->setLineWidth(2);
      QVBoxLayout *infoframelayout = new QVBoxLayout(infoframe, 10, 5,
                                                     "infoframelayout");
      mainlayout->addWidget(infoframe);

      // a 3x2 QGridLayout in the info frame (space=5)
      QGridLayout *numbers_layout = new QGridLayout(infoframelayout,
                                                    3, 2,
                                                    5, // space
                                                    "infolayout");
      //   number of points
      numbers_layout->addWidget(new QLabel("Number of points: ",
                                           infoframe),
                                0, 0,
                                AlignRight | AlignTop );

      nb_of_points = new QLabel("0", infoframe);
      numbers_layout->addWidget(nb_of_points, 0, 1,
                                AlignLeft | AlignTop );

      //   number of clusters
      numbers_layout->addWidget(new QLabel("Number of clusters: ",
                                           infoframe),
                                1, 0,
                                AlignRight | AlignTop );

      nb_of_clusters = new QLabel("", infoframe);
      numbers_layout->addWidget(nb_of_clusters, 1, 1,
                                AlignLeft | AlignTop );

      //   initialization status
      numbers_layout->addWidget(new QLabel("init status: ", infoframe),
                                2, 0,
                                AlignRight | AlignTop );

      init_status = new QLabel("no", infoframe);
      numbers_layout->addWidget(init_status, 2, 1,
                                AlignLeft | AlignTop );

      // a vertical spacer in the info frame
      infoframelayout->addItem(new
                               QSpacerItem( 0, 0,
                                            QSizePolicy::Minimum,
                                            QSizePolicy::Expanding ));

      // another grid: a 2x3 QGridLayout (space=5)
      QGridLayout *criteria_layout =
        new QGridLayout(infoframelayout, 2, 3,
                        5, // space
                        "criteria_layout");

      //    angle bound
      criteria_layout->addWidget(new QLabel("Angle bound: ",
                                            infoframe),
                                 0, 0,
                                 AlignRight | AlignTop);
      angle_bound = new QLineEdit("0.125", infoframe,
                                  "angle_bound");
      angle_bound->setValidator(new QDoubleValidator(angle_bound));
      criteria_layout->addWidget(angle_bound, 0, 1,
                                 AlignLeft | AlignTop );
      connect(angle_bound, SIGNAL(textChanged(const QString&)),
              this, SLOT(setBound(const QString&)));


      //    size bound
      criteria_layout->addWidget(new QLabel("Size bound: ",
                                            infoframe),
                                 1, 0,
                                 AlignRight | AlignTop);
      size_bound = new QLineEdit("0", infoframe,
                                 "size_bound");
      size_bound->setValidator(new QDoubleValidator(size_bound));
      criteria_layout->addWidget(size_bound, 1, 1,
                                 AlignLeft | AlignTop );
      connect(size_bound, SIGNAL(textChanged(const QString&)),
              this, SLOT(setSizeBound(const QString&)));

      //    under mouse
      under_mouse = new QCheckBox("Under mouse only", infoframe,
                                  "under_mouse");
      criteria_layout->addMultiCellWidget(under_mouse,
                                          3, 3, 0, 1);
      connect(under_mouse, SIGNAL(toggled(bool)),
              this, SLOT(setLocal(bool)));

      setCentralWidget(mainframe);
      resize(700,500);
      mainframe->show();

      // --- STATUSBAR ---
      statusBar();

      // --- LAYERS ---

      show_points =
        new Show_points_from_triangulation(tr,
                                           &Tr::finite_vertices_begin,
                                           &Tr::finite_vertices_end,
                                           CGAL::BLACK, 2);

      show_seeds = new Show_seeds(mesher,
                                  &Mesher::seeds_begin,
                                  &Mesher::seeds_end,
                                  CGAL::BLUE,
                                  5,
                                  CGAL::CROSS);
      show_triangulation =
        new CGAL::Qt_layer_show_triangulation<Tr>(tr,
                                                    CGAL::BLUE,1);
      show_marked =
        new Show_marked_faces<Tr>(tr, CGAL::GREEN);
      show_constraints =
        new CGAL::Qt_layer_show_triangulation_constraints<Tr>
        (tr, CGAL::RED, 1);
      show_circles =
        new CGAL::Qt_layer_show_circles<Tr>(tr, CGAL::GRAY, 1);
      show_mouse = new CGAL::Qt_widget_show_mouse_coordinates(*this);

      show_clusters = new Show_clusters<Mesher>(*mesher,
                                              CGAL::BLACK,3,CGAL::RECT,
                                              CGAL::BLACK,2);

      // layers order, first attached are "under" last attached
      widget->attach(show_marked);
      widget->attach(show_triangulation);
      widget->attach(show_constraints);
      widget->attach(show_circles);
      widget->attach(show_points);
      widget->attach(show_clusters);
      widget->attach(show_seeds);
      widget->attach(show_mouse);

      show_circles->deactivate();
      show_clusters->deactivate();

      get_point = new CGAL::Qt_widget_get_point<K>();
      widget->attach(get_point);
      get_point->deactivate();

      get_polygon = new CGAL::Qt_widget_get_polygon<CGALPolygon>();
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


      // --- TOOLBARS ---

      // Actions: bouding box and mesh
      QToolBar *toolBarActions = new QToolBar("Actions", this);
      QPushButton *pbBounding =
        new QPushButton("Insert bounding box", toolBarActions);
      connect(pbBounding, SIGNAL(clicked()), this,
              SLOT(insert_bounding_box()));

      QPushButton *pbMesh =
        new QPushButton("Mesh", toolBarActions);
      connect(pbMesh, SIGNAL(clicked()), this, SLOT(refineMesh()));

      QPushButton *pbConform =
        new QPushButton("Conform", toolBarActions);
      connect(pbConform, SIGNAL(clicked()), this,
              SLOT(conformMesh()));

      QPushButton *pbAdvanced =
        new QPushButton("Advanced", toolBarActions);
      pbAdvanced->setToggleButton(true);
      pbAdvanced->setOn(false);
      connect(pbAdvanced, SIGNAL(stateChanged(int)),
              this, SLOT(advanced(int)));

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

      // Layers: points, edges, constrained edges
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

      // Steps actions: step by step meshing operations
      toolBarAdvanced = new QToolBar("Advanced operations",this);
      toolBarAdvanced->hide();

      pbMeshStep =
        new QPushButton("Mesh 1 step", toolBarAdvanced);
      connect(pbMeshStep, SIGNAL(clicked()), this,
              SLOT(refineMeshStep()));

      QSpinBox *sbStepLenght =
        new QSpinBox(1, INT_MAX, 1, toolBarAdvanced);
      sbStepLenght->setValue(1);
      step_lenght = 1;
      connect(sbStepLenght, SIGNAL(valueChanged(int)),
              this, SLOT(updateStepLenght(int)));

      timer = new QTimer(this);
      connect(timer, SIGNAL(timeout()),
              this, SLOT(refineMeshStep()));

      pbMeshTimer = new QPushButton("Auto step", toolBarAdvanced);
      pbMeshTimer->setToggleButton(true);
      connect(pbMeshTimer, SIGNAL(stateChanged(int)),
              this, SLOT(updateTimer(int)));

      QSpinBox *sbTimerInterval =
        new QSpinBox(0, INT_MAX, 10, toolBarAdvanced);
      sbTimerInterval->setValue(1000);
      sbTimerInterval->setSuffix("ms");
      timer_interval=1000;
      connect(sbTimerInterval, SIGNAL(valueChanged(int)),
              this, SLOT(updateTimerInterval(int)));

      QPushButton* pbShowCluster =
        new QPushButton("Show clusters", toolBarAdvanced);
      pbShowCluster->setToggleButton(true);
      connect(pbShowCluster, SIGNAL(stateChanged(int)),
              show_clusters, SLOT(stateChanged(int)));
      connect(pbShowCluster, SIGNAL(stateChanged(int)),
              widget, SLOT(redraw()));

      setUsesBigPixmaps(true);

      // --- MENUS ---
      QPopupMenu *pmMesh = new QPopupMenu(this);
      menuBar()->insertItem("&File", pmMesh);
      pmMesh->insertItem("&Refine mesh", this, SLOT(refineMesh()),
                         CTRL+Key_R );
      pmMesh->insertItem("&Clear mesh", this, SLOT(clearMesh()),
                         CTRL+Key_C );
      pmMesh->insertItem("Clear seeds", this, SLOT(clearSeeds()));
      pmMesh->insertItem("&Open constrained triangulation...", this,
                         SLOT(openTriangulation()),
                         CTRL+Key_O );
      pmMesh->insertItem("&Save constrained edges...", this,
                         SLOT(saveTriangulation()),
                         CTRL+Key_S );
      pmMesh->insertItem("&Quit", qApp, SLOT(closeAllWindows()),
                         CTRL+Key_Q );

      connect(this, SIGNAL(insertedInput()),
              this, SLOT(after_inserted_input()));
      connect(this, SIGNAL(insertedInput()),
              show_clusters, SLOT(reinitClusters()));
      connect(this, SIGNAL(initializedMesh()),
              this, SLOT(after_initialized_mesh()));

      widget->set_window(-1.,1.,-1.,1.);
      widget->setMouseTracking(TRUE);
    };

  // compute bounds of the mesh
  void bounds(FT &xmin, FT &ymin,
              FT &xmax, FT &ymax)
    {
      Tr::Finite_vertices_iterator vi=tr.finite_vertices_begin();
      xmin=xmax=vi->point().x();
      ymin=ymax=vi->point().y();
      vi++;
      while(vi != tr.finite_vertices_end())
        {
          if(vi->point().x() < xmin) xmin=vi->point().x();
          if(vi->point().x() > xmax) xmax=vi->point().x();
          if(vi->point().y() < ymin) ymin=vi->point().y();
          if(vi->point().y() > ymax) ymax=vi->point().y();
          vi++;
        }
    }

signals:
  void insertedInput();
  void initializedMesh();

private slots:
  void after_initialized_mesh()
  {
    updatePointCounter();
    switch( mesh->get_initialized() )
      {
      case Mesh::NONE: init_status->setText("no"); break;
      case Mesh::DELAUNAY: init_status->setText("Delaunay"); break;
      case Mesh::GABRIEL: init_status->setText("Gabriel"); break;
      default: init_status->setText("unknown"); break;
      }
    show_clusters->reinitClusters();
  }

  void after_inserted_input()
  {
    mesh->set_initialized(Mesh::NONE);
    emit initializedMesh();
    nb_of_clusters_has_to_be_updated = true;
    mesh->mark_facets();
  }

public slots:

  void get_cgal_object(CGAL::Object obj)
    {
      Point_2 p;
      CGALPolygon poly;

      if(CGAL::assign(p,obj))
        if(follow_mouse->is_active())
          {
            typedef Mesh::Face_handle Face_handle;
            std::list<Face_handle> l;

            Face_handle fh = mesh->locate(p);
            criteria.set_point(p);
            if( (fh!=NULL) && (!mesh->is_infinite(fh)) && fh->is_marked() )
              {
                const Point_2&
                  a = fh->vertex(0)->point(),
                  b = fh->vertex(1)->point(),
                  c = fh->vertex(2)->point();

                Mesh::Quality q;
                if(criteria.is_bad_object().operator()(a, b, c, q))
                  l.push_back(fh);
              }
            mesh->set_criteria(criteria);
            mesh->set_bad_faces(l.begin(), l.end());
            while( mesh->step_by_step_refine_mesh() );
          }
        else
          if(get_seed->is_active())
            {
              Mesh::Seeds seeds;
              std::copy(mesh->seeds_begin(), mesh->seeds_end(),
                        std::back_inserter(seeds));
              seeds.push_back(p);
              mesh->set_seeds(seeds.begin(), seeds.end());
              mesh->mark_facets();
            }
          else // get_point is active
            {
              mesh->insert(p);
              emit( insertedInput() );
            }
      else
        if (CGAL::assign(poly,obj))
          {
            for(CGALPolygon::Edge_const_iterator it=poly.edges_begin();
                it!=poly.edges_end();
                it++)
              mesh->insert((*it).source(),(*it).target());
            emit( insertedInput() );
          }
        else // obj should be a polygon of or point!
          CGAL_assertion(false);
      updatePointCounter();
      widget->redraw();
    }

  //insert a bounding box around the mesh
  void insert_bounding_box()
    {
      FT xmin, xmax, ymin, ymax;
      bounds(xmin, ymin, xmax, ymax);

      FT xcenter=(xmin+xmax)/2,
        ycenter=(ymin+ymax)/2;
      FT xspan = (xmax-xmin)/2,
        yspan = (ymax-ymin)/2;

      Point_2 bb1(xcenter - FT(1.5)*xspan, ycenter - FT(1.5)*yspan);
      Point_2 bb2(xcenter + FT(1.5)*xspan, ycenter - FT(1.5)*yspan);
      Point_2 bb3(xcenter + FT(1.5)*xspan, ycenter + FT(1.5)*yspan);
      Point_2 bb4(xcenter - FT(1.5)*xspan, ycenter + FT(1.5)*yspan);
      mesh->insert(bb1);
      mesh->insert(bb2);
      mesh->insert(bb3);
      mesh->insert(bb4);
      mesh->insert(bb1, bb2);
      mesh->insert(bb2, bb3);
      mesh->insert(bb3, bb4);
      mesh->insert(bb4, bb1);
      emit( insertedInput() );
      widget->redraw();
    }

  void updatePointCounter()
    {
      nb_of_points->setNum(static_cast<int>(mesh->number_of_vertices()));
      if(nb_of_clusters_has_to_be_updated &&
         mesh->get_initialized() != Mesh::NONE)
        {
          nb_of_clusters_has_to_be_updated = false;
          nb_of_clusters->setNum(mesh->number_of_clusters_vertices());
        }
    }

  void refineMesh()
    {
      dumpTriangulation("last_input.edg");
      mesh->refine_mesh();
      emit initializedMesh();
      widget->redraw();
    }

  void conformMesh()
    {
      if(mesh->get_initialized() != Mesh::GABRIEL)
        dumpTriangulation("last_input.edg");
      mesh->make_conforming_Gabriel();
      initializedMesh();
      updatePointCounter();
      widget->redraw();
    }

  void refineMeshStep()
    {
      CGAL_assertion(mesh->is_valid(true)); // @todo: remove this
      int counter = step_lenght;
      if(mesh->get_initialized() != Mesh::GABRIEL)
        {
          mesh->init_Gabriel();
          initializedMesh();
          dumpTriangulation("last_input.edg");
        }
      while(counter>0)
        {
          --counter;
          if(!mesh->step_by_step_refine_mesh())
            {
              pbMeshTimer->setOn(false);
              counter = 0;
            }
        }
      CGAL_assertion(mesh->is_valid(true)); // @todo: remove this
      updatePointCounter();
      widget->redraw();
    }

  void updateTimer(int i)
    {
      if(i==0)
        timer->stop();
      else
        timer->start(timer_interval);
    }

  void updateStepLenght(int i)
    {
      step_lenght = i;
      QString s;
      s = "Mesh " + QString::number(i) + " step";
      if(i > 1)
        s+="s";
      pbMeshStep->setText(s);
    }

  void updateTimerInterval(int i)
    {
      timer_interval=i;
      if(timer->isActive())
        timer->changeInterval(timer_interval);
    }

  void clearMesh()
    {
      mesh->clear();
      emit( insertedInput() );
      updatePointCounter();
      widget->clear_history();
      widget->redraw();
    }

  void clearSeeds()
    {
      mesh->clear_seeds();
      mesh->mark_facets();
      widget->redraw();
    }

  void openTriangulation() { openTriangulation(QString()); }

  void openTriangulation(QString filename)
    {
      QString s;
      if( filename.isEmpty() )
        s = QFileDialog::getOpenFileName( QString::null,
                                          my_filters, this );
      else
        s = filename;

      if ( s.isEmpty() )
        return;
      std::ifstream f(s);
      if (!f) return;

      if(s.right(5) == ".poly")
        {
          Mesh::Seeds seeds;
          CGAL::read_poly(*mesh, f, std::back_inserter(seeds));
          mesh->set_seeds(seeds.begin(), seeds.end(), false);
        }
      else if(s.right(5) == ".data")
        {
          int nx, ny, niso, use_threshold;
          float threshold;

          std::ifstream ins(s);
          ins >> nx >> ny >> niso >> use_threshold >> threshold;
          for(int c = 0; c < niso; c++) {
            float f;
            ins >> f;
          }

          std::vector<Point_2> points(nx * ny);
          double xmin,xmax,ymin,ymax;
          ins >> xmin >> xmax >> ymin >> ymax;

          double dx = (xmax-xmin)/(nx-1);
          double dy = (ymax-ymin)/(ny-1);

          int k2=0;
          for (int i2=0; i2<nx; i2++) {
            for (int j=0; j<ny; j++) {
              points[k2] = Point_2(xmin + i2*dx,  ymin + j*dy);
              k2++;
            }
          }

          std::random_shuffle(points.begin(), points.end());
          mesh->clear();
          //      std::copy(points.begin(), points.end(), std::back_inserter(*mesh));

          s.replace(QRegExp("\\.data$"), "_fault.data");

          std::ifstream ins2(s);
          int num_lines;
          ins2 >> num_lines;
          std::vector<int> num_vertex_per_line(num_lines);
          for(int n = 0; n < num_lines; n++){
            ins2 >> num_vertex_per_line[n];
          }

          CGAL::Bbox_2 b;
          for(int i = 0; i < num_lines; i++){
            Point_2 p, q;
            ins2 >> p;
            if(i == 0){
              b = p.bbox();
            } else {
              b = b + p.bbox();
            }
            for(int j = 1; j < num_vertex_per_line[i]; j++){
              ins2 >> q;
              mesh->insert_constraint(p, q);
              p = q;
              b = b + p.bbox();
            }
          }

          for(unsigned int k = 0; k < points.size(); k++)
            if (CGAL::do_overlap(b,points[k].bbox()))
              mesh->insert(points[k]);

          xmax = b.xmax();
          xmin = b.xmin();
          ymax = b.ymax();
          ymin = b.ymin();

          dx = (xmax - xmin)/20.0;
          dy = (ymax - ymin)/20.0;
          xmin -= dx;
          ymin -= dy;
          xmax += dx;
          ymax += dy;
          Point_2 bl(xmin, ymin);
          Point_2 br(xmax, ymin);
          Point_2 tl(xmin, ymax);
          Point_2 tr(xmax, ymax);
          mesh->insert_constraint(bl, br);
          mesh->insert_constraint(br, tr);
          mesh->insert_constraint(tr, tl);
          mesh->insert_constraint(tl, bl);

          clearSeeds();
        }
      else
        {
          read_constraints(*mesh, f);
          clearSeeds();
        }

      // compute bounds
      FT xmin, xmax, ymin, ymax;
      bounds(xmin, ymin, xmax, ymax);

      FT xspan = (xmax-xmin)/2,
        yspan = (ymax-ymin)/2;

      widget->set_window(CGAL::to_double(xmin-FT(1.1)*xspan),
                         CGAL::to_double(xmax+FT(1.1)*xspan),
                         CGAL::to_double(ymin-FT(1.1)*yspan),
                         CGAL::to_double(ymax+FT(1.1)*yspan));
      widget->clear_history();

      emit( insertedInput() );
      updatePointCounter();
      widget->redraw();
    }

  void saveTriangulation()
    {
      QString s( QFileDialog::getSaveFileName( "filename.edg",
                                               my_filters, this ) );
      if ( s.isEmpty() )
        return;
      std::ofstream of(s);
      if(s.right(5) == ".poly")
        CGAL::write_poly(*mesh, of);
      else
        write_constraints(*mesh, of);
    }

  void dumpTriangulation(QString s=QString("dump.edg"))
    {
      std::ofstream of(s);
      write_constraints(*mesh, of);
    }

  inline
  void fake_slot()
    {
    }

  void setBound(const QString& bound)
    {
      criteria.set_bound(bound.toDouble());
      if( mesh != 0 )
        mesh->set_criteria(criteria);
    }

  void setSizeBound(const QString& size_bound)
    {
      criteria.set_size_bound(size_bound.toDouble());
      if ( mesh != 0 )
        mesh->set_criteria(criteria);
    }

  void setLocal(bool checked)
    {
      criteria.set_local_size(checked);
      mesh->set_criteria(criteria);
      if(criteria.is_local_size())
        follow_mouse->activate();
      else
        follow_mouse->deactivate();
    }

  void advanced(int state)
  {
    if( state == 0 )
      toolBarAdvanced->hide();
    else
      toolBarAdvanced->show();
  }

private:
  static const QString my_filters;
  Criteria criteria;
  Tr tr;
  Mesh* mesh;

  QPopupMenu *pmCriteria;
  int menu_id;

  CGAL::Qt_widget* widget;
  CGAL::Qt_widget_get_point<K>* get_point;
  CGAL::Qt_widget_get_point<K>* get_seed;
  CGAL::Qt_widget_get_polygon<CGALPolygon>* get_polygon;
  Follow_mouse* follow_mouse;

  typedef CGAL::Qt_layer_show_points<Tr, Tr::Finite_vertices_iterator,
  Vertex_to_point>
    Show_points_from_triangulation;

  typedef CGAL::Qt_layer_show_points<Tr, Tr::Seeds_const_iterator>
    Show_seeds;

  Show_points_from_triangulation* show_points;
  Show_seeds* show_seeds;
  CGAL::Qt_layer_show_triangulation<Tr>* show_triangulation;
  CGAL::Qt_layer_show_triangulation_constraints<Tr>* show_constraints;
  CGAL::Qt_layer_show_circles<Tr>* show_circles;
  CGAL::Qt_widget_show_mouse_coordinates* show_mouse;
  Show_marked_faces<Tr>* show_marked;

  bool nb_of_clusters_has_to_be_updated;
  QLabel *nb_of_clusters;
  Show_clusters<Mesher>* show_clusters;

  QLabel *nb_of_points;
  QLabel *init_status;
  QLineEdit* angle_bound;
  QLineEdit* size_bound;
  QCheckBox* under_mouse;
  QTimer* timer;
  QPushButton *pbMeshTimer;
  QPushButton *pbMeshStep;
  QToolBar *toolBarAdvanced;
  int timer_interval;
  int step_lenght;
};

const QString MyWindow::my_filters =
"Constrained edges (*.edg);;"
"Shewchuk Triangle .poly files (*.poly);;"
"All files (*)";


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

  if( argc == 2 )
    W->openTriangulation(QString(argv[1]));

  my_previous_failure_function =
    CGAL::set_error_handler(cgal_with_exceptions_failure_handler);

  try {
    return app.exec();
  }
  catch(Cgal_exception e) {
    std::cerr << "catch(Cgal_exception e)" << std::endl;
    try {
      W->dumpTriangulation();
    }
    catch(...) {
      std::cerr << "PANIC !!" << std::endl;
    }
    my_previous_failure_function(e.type, e.expr, e.file, e. line, e.msg);
  }

  return 0;
}

// moc_source_file: mesh_demo.C
#include "mesh_demo.moc"

// moc_source_file: Show_clusters.h
#include "Show_clusters.moc"

#endif // CGAL_USE_QT
