// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#include <CGAL/basic.h>


#include <climits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#define CGAL_MESH_2_USE_TIMERS
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_local_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

#include <CGAL/IO/File_poly.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_get_point.h>
#include <CGAL/IO/Qt_widget_get_polygon.h>

#include "icons.h"

#include "Show_points.h"
#include "Show_segments.h"
#include "Qt_layer_show_triangulation.h"
#include "Qt_layer_show_triangulation_constraints.h"
#include "Qt_layer_show_circles.h"
#include "Show_clusters.h"
#include <CGAL/IO/Qt_widget_show_mouse_coordinates.h>
#include "Qt_widget_styled_layer.h"

#ifdef CGAL_MESH_2_DEBUG_DRAW
#include "Debug_layer.h"
#endif

#include "Qt_widget_style_editor.h"

#include <qapplication.h>
#include <qmainwindow.h>
#include <qcolor.h>
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
#include <qtabwidget.h>

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
typedef K::Segment_2 Segment_2;
typedef K::Triangle_2 Triangle_2;
typedef K::Circle_2 Circle;
typedef CGAL::Polygon_2<K> CGALPolygon;

typedef CGAL::Delaunay_mesher_2<Tr, Criteria> Mesher;
typedef Tr::Vertex Vertex;
typedef Tr::Edge Edge;

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

struct Vertex_to_point {
  const Point_2& operator()(const Vertex& v) const
  {
    return v.point();
  }
};

struct Edge_to_segment {
  const Segment_2 operator()(const Edge& e) const
  {
    return Segment_2(e.first->vertex(Tr::cw (e.second))->point(),
                     e.first->vertex(Tr::ccw(e.second))->point());
  }
};

struct Background_color : public CGAL::Qt_widget_styled_layer
{
  Background_color (CGAL::Color c = CGAL::WHITE)
  {
    style()->setColor("Back ground color",
                      QColor(c.red(), c.green(), c.blue()));
  }

  void draw()
  {
    widget->setBackgroundColor(style()->getColor("Back ground color"));
  }
};

template <class Tr>
class Show_in_domain_faces : public CGAL::Qt_widget_styled_layer
{
  Tr *cdt;
public:
  Show_in_domain_faces(Tr *t, CGAL::Color c=CGAL::GREEN,
                       QObject* parent = 0, const char* name = 0)
    : Qt_widget_styled_layer(0, parent, name),
      cdt(t)
  {
    QColor color(c.red(), c.green(), c.blue());
    style()->setColor("Fill color", color);
    style()->setInt("Line width", 0);
    style()->setColor("Segment color", color);
  }

  typedef typename Tr::Finite_faces_iterator Face_iterator;

  void draw()
  {
    QColor old_fill_color = widget->fillColor();
    int old_line_width = widget->lineWidth();
    QColor old_line_color = widget->color();

    widget->setColor(style()->getColor("Segment color"));
    widget->setLineWidth(style()->getInt("Line width"));
    widget->setFillColor(style()->getColor("Fill color"));

    for(Face_iterator fit=cdt->finite_faces_begin();
        fit!=cdt->finite_faces_end();
        ++fit)
      if(fit->is_in_domain())
        *widget << cdt->triangle(fit);

    widget->setFillColor(old_fill_color);
    widget->setLineWidth(old_line_width);
    widget->setColor(old_line_color);
  }
};

template <class Mesher>
class Show_bad_faces : public CGAL::Qt_widget_layer
{
  Mesher *m;
  CGAL::Color color;
public:
  Show_bad_faces(Mesher *mesher,
                 CGAL::Color c=CGAL::Color(127,0,0),
                 QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name),
      m(mesher),
      color(c)
  {}

  void set_container(Mesher* mesher)
  {
    m = mesher;
  }

  typedef typename Mesher::Bad_faces_const_iterator Bad_faces_iterator;

  void draw()
  {
    if( m != 0 )
      {
        widget->lock();
        {
          QColor old_fill_color = widget->fillColor();
          int old_line_width = widget->lineWidth();

          *widget << CGAL::FillColor(color)
                  << CGAL::LineWidth(0);

          for(Bad_faces_iterator fit=m->bad_faces_begin();
                fit!=m->bad_faces_end();
              ++fit)
            *widget << Triangle_2((*fit)->vertex(0)->point(),
                                  (*fit)->vertex(1)->point(),
                                  (*fit)->vertex(2)->point());
          widget->setFillColor(old_fill_color);
          widget->setLineWidth(old_line_width);
        }
        widget->unlock();
      }
  }
};

class Meshing_debugging_layer : public CGAL::Qt_widget_layer
{
  Q_OBJECT
  Mesher* m;
  // bool point_on;
public:
  Meshing_debugging_layer(QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name), m(0)/* , point_on(false) */
  {
  }

  void set_container(Mesher* mesher)
  {
    m = mesher;
  }

  void draw()
  {
    if( m != 0 )
      {
        if( m->is_refinement_done() ) return;

        widget->lock();
        {
          QColor old_color = widget->color();
          int old_line_width = widget->lineWidth();
          CGAL::PointStyle old_point_style = widget->pointStyle();
          int old_point_size = widget->pointSize();

          if( m->is_edges_refinement_done() )
            {
              QColor old_fill_color = widget->fillColor();

              const Tr::Face_handle fh = m->next_bad_face();

              *widget << CGAL::FillColor(CGAL::RED)
                      << CGAL::LineWidth(0)
                      << Triangle_2(fh->vertex(0)->point(),
                                    fh->vertex(1)->point(),
                                    fh->vertex(2)->point());

              widget->setFillColor(old_fill_color);
            }
          else
            {
              Edge_to_segment edge_to_segment;

              *widget << CGAL::PURPLE
                      << CGAL::LineWidth(2)
                      << edge_to_segment(m->next_encroached_edge());
            }

          *widget << CGAL::BLACK
                  << CGAL::PointStyle(CGAL::PLUS)
                  << CGAL::PointSize(3)
                  << m->next_refinement_point()
                  << CGAL::PointSize(old_point_size)
                  << CGAL::PointStyle(old_point_style);

          widget->setLineWidth(old_line_width);
          widget->setColor(old_color);
        }
        widget->unlock();
      }
  }
};


class Follow_mouse : public CGAL::Qt_widget_layer
{
  QCursor oldcursor;
public:
  Follow_mouse(QObject* parent = 0, const char* name = 0)
    : Qt_widget_layer(parent, name)
  {
  }

  void mouseMoveEvent(QMouseEvent *e)
  {
    FT x=static_cast<FT>(widget->x_real(e->x()));
    FT y=static_cast<FT>(widget->y_real(e->y()));
    const Point_2 current_point = Point_2(x, y);
    widget->new_object(CGAL::make_object(current_point));
    previous_point = current_point;
    just_entered = false;
  }

  void enterEvent(QEvent*)
  {
    just_entered = true;
  }

  void leaveEvent(QEvent*)
  {
    just_entered = false;
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

  Point_2 previous_point;
  bool just_entered;
};

class Preferences : public QWidget
{
  Q_OBJECT
private:
  QGridLayout* layout;
  QTabWidget* tab;

public:
  Preferences(QWidget* parent = 0, const char* name = 0, bool modal = false)
    : QWidget(parent, name, modal)
  {
    layout = new QGridLayout(this, 1,1, 11, 6, "pref_layout");
    tab = new QTabWidget(this);
    layout->addWidget(tab, 0, 0);
  }

  template <typename LayersIterator>
  void setLayers(LayersIterator begin, LayersIterator end)
  {
    QWidget* w;
    while( (w = tab->currentPage()) != 0)
      {
	tab->removePage(w);
	delete w;
      }

    for(LayersIterator it = begin; it!=end; ++it)
      {
	QFrame* box = new QFrame(this);
	QHBoxLayout* hor = new QHBoxLayout(box, 5, 0);
	QVBoxLayout* layout = new QVBoxLayout(hor, 5);

	QCheckBox* check = new QCheckBox("&Activated", box);
	layout->addWidget(check);

	check->setChecked((*it)->is_active());

	connect(check, SIGNAL(stateChanged(int)),
		*it, SLOT(stateChanged(int)));
	connect(check, SIGNAL(stateChanged(int)),
		this, SLOT(is_changed()));

	//	QGroupBox* group = new QGroupBox("&Properties", box);
	CGAL::Qt_widget_style_editor* editor =
	  new CGAL::Qt_widget_style_editor((*it)->style(), box);
	editor->show();

	layout->addWidget(editor);
	layout->addItem(new QSpacerItem(0, 0,
					QSizePolicy::Minimum,
					QSizePolicy::Expanding));
	hor->addItem(new QSpacerItem(0, 0,
				     QSizePolicy::Expanding,
				     QSizePolicy::Minimum));

	connect(editor, SIGNAL(styleChanged()),
		this, SLOT(is_changed()));

	tab->addTab(box, (*it)->name());
      }
  }

private slots:
  void is_changed() { emit changed(); }

signals:
  void changed();
}; // end of class Preferences

class MyWindow : public QMainWindow
{
  Q_OBJECT

  typedef std::list<Point_2> Seeds;
public:
  MyWindow() : criteria(), mesher(0)
    {
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

      Background_color* background_color_layer = new Background_color();

      show_points =
        new Show_points_from_triangulation(&cdt,
                                           &Tr::finite_vertices_begin,
                                           &Tr::finite_vertices_end,
                                           CGAL::BLACK, 2,
                                           CGAL::DISC,
                                           this, "Show points");

      show_seeds = new Show_seeds(&seeds,
                                  &Seeds::begin,
                                  &Seeds::end,
                                  CGAL::BLUE,
                                  5,
                                  CGAL::CROSS,
                                  this, "Show seeds");
      show_triangulation =
        new CGAL::Qt_layer_show_triangulation<Tr>(&cdt,
                                                  CGAL::BLUE,1,
                                                  this,
                                                  "Show triangulation edges");
      show_in_domain =
        new Show_in_domain_faces<Tr>(&cdt, CGAL::GREEN,
                                     this, "Show in_domain faces");

      show_constraints =
        new CGAL::Qt_layer_show_triangulation_constraints<Tr>
        (&cdt, CGAL::RED, 1,
         this, "Show constrained edges");

      show_encroached_edges =
        new Show_encroached_edges(0,
                                  &Mesher::encroached_edges_begin,
                                  &Mesher::encroached_edges_end,
                                  CGAL::RED, 2,
                                  this, "Encroached edges layer");

      show_bad_faces =
        new Show_bad_faces<Mesher>(0, CGAL::Color(0,160,0),
                                   this, "Show bad faces");

      debug_layer = new Meshing_debugging_layer(this,
                                                "Debug layer");

      show_circles =
        new CGAL::Qt_layer_show_circles<Tr>(&cdt, CGAL::GRAY, 1,
                                            CGAL::WHITE, false,
                                            this, "Show circles");

      show_coordinates =
        new CGAL::Qt_widget_show_mouse_coordinates(*this,
                                                   this,
                                                   "Follow mouse coordinates");

      show_clusters = new Show_clusters<Mesher>(mesher,
                                                CGAL::BLACK,3,CGAL::RECT,
                                                CGAL::BLACK,CGAL::BLUE,2,
                                                this,
                                                "Show clusters");

      // layers order, first attached are "under" last attached
#ifdef CGAL_MESH_2_DEBUG_DRAW
      widget->attach(new CGAL::Debug_layer());
#endif
      widget->attach(background_color_layer);
      widget->attach(show_in_domain);
      widget->attach(show_bad_faces);
      widget->attach(show_triangulation);
      widget->attach(show_constraints);
      widget->attach(show_circles);
      widget->attach(show_encroached_edges);
      widget->attach(show_points);
      widget->attach(show_seeds);
      widget->attach(show_coordinates);
      widget->attach(debug_layer);
      widget->attach(show_clusters); // should be last

      show_circles->deactivate();
      show_encroached_edges->deactivate();
      show_bad_faces->deactivate();
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

      const int number_of_styled_layers = 6;
      CGAL::Qt_widget_styled_layer* styled_layers[number_of_styled_layers] =
        { background_color_layer,
          show_points,
          show_seeds,
          show_triangulation,
          show_constraints,
          show_in_domain};

      prefs = new Preferences(0, "Preferences", false);
      prefs->setCaption("Layers properties");
      prefs->setLayers(styled_layers, styled_layers + number_of_styled_layers);
      prefs->resize(300, 200);

      connect(prefs, SIGNAL(changed()),
	      widget, SLOT(redraw()));

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

      QToolButton *pbShowInDomain
        = new QToolButton(QPixmap( (const char**)marked_xpm ),
                          "Show faces in domain",
                          "Display faces that will be refined",
                          this, SLOT(fake_slot()),
                          toolbarLayers,
                          "show in domain");
      pbShowInDomain->setToggleButton(true);
      pbShowInDomain->setOn(true);
      connect(pbShowInDomain, SIGNAL(stateChanged(int)),
              show_in_domain, SLOT(stateChanged(int)));

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
      bgLayers->insert(pbShowInDomain);
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

      QSpinBox *sbStepLength =
        new QSpinBox(1, INT_MAX, 1, toolBarAdvanced);
      sbStepLength->setValue(1);
      step_length = 1;
      connect(sbStepLength, SIGNAL(valueChanged(int)),
              this, SLOT(updateStepLength(int)));

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

      pbShowCluster = new QPushButton("Show clusters", toolBarAdvanced);
      pbShowCluster->setToggleButton(true);
      pbShowCluster->setEnabled(false);
      connect(pbShowCluster, SIGNAL(stateChanged(int)),
              show_clusters, SLOT(stateChanged(int)));
      connect(pbShowCluster, SIGNAL(stateChanged(int)),
              widget, SLOT(redraw()));

      pbShowEncroachedEdges = new QPushButton("Show encroached edges",
                                              toolBarAdvanced);
      pbShowEncroachedEdges->setToggleButton(true);
      pbShowEncroachedEdges->setEnabled(false);
      connect(pbShowEncroachedEdges, SIGNAL(stateChanged(int)),
              show_encroached_edges, SLOT(stateChanged(int)));
      connect(pbShowEncroachedEdges, SIGNAL(stateChanged(int)),
              widget, SLOT(redraw()));

      pbShowBadFaces = new QPushButton("Show bad faces",
                                       toolBarAdvanced);
      pbShowBadFaces->setToggleButton(true);
      pbShowBadFaces->setEnabled(false);
      connect(pbShowBadFaces, SIGNAL(stateChanged(int)),
              show_bad_faces, SLOT(stateChanged(int)));
      connect(pbShowBadFaces, SIGNAL(stateChanged(int)),
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
      connect(this, SIGNAL(initializedMesher()),
              this, SLOT(after_initialized_mesher()));

      QPopupMenu *pmOptions = new QPopupMenu(this);
      menuBar()->insertItem("&Options", pmOptions);
      pmOptions->insertItem("Set &layer properties...", this,
                            SLOT(displayPreferences()));
      pmOptions->setCheckable(true);

      widget->set_window(-1.,1.,-1.,1.);
      widget->setMouseTracking(TRUE);
    };

  // compute bounds of the mesh
  void bounds(FT &xmin, FT &ymin,
              FT &xmax, FT &ymax)
    {
      Tr::Finite_vertices_iterator vi=cdt.finite_vertices_begin();
      xmin=xmax=vi->point().x();
      ymin=ymax=vi->point().y();
      vi++;
      while(vi != cdt.finite_vertices_end())
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
  void initializedMesher();

private slots:
  void after_initialized_mesher()
  {
    updatePointCounter();
    if( mesher !=0 )
      {
        init_status->setText("yes");
        pbShowCluster->setEnabled(true);
        pbShowEncroachedEdges->setEnabled(true);
        pbShowBadFaces->setEnabled(true);
      }
    else
      {
        init_status->setText("no");
        pbShowCluster->setOn(false);
        pbShowCluster->setEnabled(false);
        pbShowEncroachedEdges->setOn(false);
        pbShowEncroachedEdges->setEnabled(false);
        pbShowBadFaces->setOn(false);
        pbShowBadFaces->setEnabled(false);
      }
    show_clusters->change_mesher(mesher);
    show_encroached_edges->set_container(mesher);
    show_bad_faces->set_container(mesher);
    debug_layer->set_container(mesher);
  }

  void after_inserted_input()
  {
    delete mesher;
    mesher = 0;
    emit initializedMesher();
    mark_facets();
    nb_of_clusters_has_to_be_updated = true;
  }

  void mark_facets()
  {
    Mesher::mark_facets(cdt, seeds.begin(), seeds.end());
  }

  void displayPreferences()
  {
    prefs->show();
  }


public slots:

  void get_cgal_object(CGAL::Object obj)
    {
      Point_2 p;
      CGALPolygon poly;

      if(CGAL::assign(p,obj))
        if(follow_mouse->is_active())
          {
            typedef Tr::Face_handle Face_handle;
            Face_handle fh = cdt.locate(p);
            std::vector<Face_handle> faces_to_check;

            if(cdt.is_infinite(fh) && !follow_mouse->just_entered)
            { // make the line walk in opposite direction, if p is outside
              // the triangulation
              fh = cdt.locate(follow_mouse->previous_point);
              std::swap(follow_mouse->previous_point, p);
            }

            if(cdt.is_infinite(fh)) return;
            Segment_2 segment;
            if(follow_mouse->just_entered)
            {
              faces_to_check.push_back(fh);
              segment = Segment_2(p, p);
            }
            else
            {
              const Point_2& previous_point = follow_mouse->previous_point;
              Tr::Line_face_circulator 
                fc = cdt.line_walk(p,
                                   previous_point,
                                   fh),
                end(fc);
              do
              {
                faces_to_check.push_back(fc);
                ++fc;
              } while(!cdt.is_infinite(fc) &&
                      cdt.triangle(fc).has_on_unbounded_side(previous_point));

              segment = Segment_2(p, follow_mouse->previous_point);
            }

            criteria.set_local_size(true);
            criteria.set_segment(segment);

            std::vector<Face_handle> bad_faces;

            for(std::vector<Face_handle>::const_iterator 
                  fh_it = faces_to_check.begin(),
                  end = faces_to_check.end();
                fh_it != end; ++fh_it)
            {
              if( (*fh_it != NULL) && (!cdt.is_infinite(*fh_it)) && (*fh_it)->is_in_domain() )
              {
                Criteria::Quality q;
                if(criteria.is_bad_object().operator()(*fh_it, q) !=
                   CGAL::Mesh_2::NOT_BAD)
                  bad_faces.push_back(*fh_it);
              }
            }

            if( mesher!=0 )
            {
              mesher->set_criteria(criteria, false);
              mesher->set_bad_faces(bad_faces.begin(), bad_faces.end());
              while( mesher->step_by_step_refine_mesh() );
            }
          }
        else
          if(get_seed->is_active())
            {
              seeds.push_back(p);
              mark_facets();
            }
          else // get_point is active
            {
              cdt.insert(p);
              emit( insertedInput() );
            }
      else
        if (CGAL::assign(poly,obj))
          {
            for(CGALPolygon::Edge_const_iterator it=poly.edges_begin();
                it!=poly.edges_end();
                it++)
              cdt.insert((*it).source(),(*it).target());
            emit( insertedInput() );
          }
        else // obj should be a polygon or a point!
          CGAL_error();
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
      cdt.insert(bb1);
      cdt.insert(bb2);
      cdt.insert(bb3);
      cdt.insert(bb4);
      cdt.insert(bb1, bb2);
      cdt.insert(bb2, bb3);
      cdt.insert(bb3, bb4);
      cdt.insert(bb4, bb1);
      emit( insertedInput() );
      widget->redraw();
    }

  void updatePointCounter()
    {
      nb_of_points->setNum(static_cast<int>(cdt.number_of_vertices()));
      if(nb_of_clusters_has_to_be_updated &&
         mesher != 0)
        {
          nb_of_clusters_has_to_be_updated = false;
          nb_of_clusters->setNum(mesher->clusters().size());
        }
    }

  void refineMesh()
    {
      dumpTriangulation("last_input.edg");
      if( mesher == 0 )
        mesher = create_mesher();
      mesher->refine_mesh();
      emit initializedMesher();
      widget->redraw();
    }

  void conformMesh()
    {
      dumpTriangulation("last_input.edg");
      CGAL::make_conforming_Gabriel_2(cdt);
      mark_facets();
      delete mesher;
      mesher = 0;
      emit initializedMesher();
      updatePointCounter();
      widget->redraw();
    }

  void refineMeshStep()
    {
      int counter = step_length;
      if(mesher == 0)
        {
          mesher = create_mesher();
          mesher->init();
          emit initializedMesher();
          dumpTriangulation("last_input.edg");
        }
      while(counter>0)
        {
          --counter;
          if(!mesher->try_one_step_refine_mesh())
            {
              pbMeshTimer->setOn(false);
              counter = 0;
            }
        }
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

  void updateStepLength(int i)
    {
      step_length = i;
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
      cdt.clear();
      emit( insertedInput() );
      updatePointCounter();
      widget->clear_history();
      widget->redraw();
    }

  void clearSeeds_without_redraw()
    {
      seeds.clear();
      delete mesher;
      mesher = 0;
      emit initializedMesher();
      mark_facets();
    }

  void clearSeeds()
    {
      clearSeeds_without_redraw();
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
      std::ifstream f(s.ascii());
      if (!f) return;

      if(s.right(5) == ".poly")
        {
          clearSeeds_without_redraw();
          CGAL::read_triangle_poly_file(cdt, f, std::back_inserter(seeds));
        }
      else if(s.right(5) == ".data")
        {
          int nx, ny, niso, use_threshold;
          float threshold;

          std::ifstream ins(s.ascii());
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
          cdt.clear();
          //      std::copy(points.begin(), points.end(), std::back_inserter(*mesh));

          s.replace(QRegExp("\\.data$"), "_fault.data");

          std::ifstream ins2(s.ascii());
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
              cdt.insert_constraint(p, q);
              p = q;
              b = b + p.bbox();
            }
          }

          for(unsigned int k = 0; k < points.size(); k++)
            if (CGAL::do_overlap(b,points[k].bbox()))
              cdt.insert(points[k]);

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
          cdt.insert_constraint(bl, br);
          cdt.insert_constraint(br, tr);
          cdt.insert_constraint(tr, tl);
          cdt.insert_constraint(tl, bl);

          clearSeeds_without_redraw();
        }
      else
        {
          read_constraints(cdt, f);
          clearSeeds_without_redraw();
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
      // widget->set_window() calls widget->redraw()
      //      widget->redraw();
    }

  void saveTriangulation()
    {
      QString s( QFileDialog::getSaveFileName( "filename.edg",
                                               my_filters, this ) );
      if ( s.isEmpty() )
        return;
      std::ofstream of(s.ascii());
      if(s.right(5) == ".poly")
        CGAL::write_triangle_poly_file(cdt, of, seeds.begin(), seeds.end());
      else
        write_constraints(cdt, of);
    }

  void dumpTriangulation(QString s=QString("dump.edg"))
    {
      std::ofstream of(s.ascii());
      write_constraints(cdt, of);
    }

  inline
  void fake_slot()
    {
    }

  void setBound(const QString& bound)
    {
      criteria.set_bound(bound.toDouble());
      if( mesher != 0 )
        mesher->set_criteria(criteria);
    }

  void setSizeBound(const QString& size_bound)
    {
      criteria.set_size_bound(size_bound.toDouble());
      if ( mesher != 0 )
        mesher->set_criteria(criteria);
    }

  void setLocal(bool checked)
    {
      if( mesher == 0 )
        {
          mesher = create_mesher();
          mesher->init();
          emit initializedMesher();
        }
      if(checked == false)
      {
        criteria.set_local_size(false);
        mesher->set_criteria(criteria);
      }
      if(checked)
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
  Mesher* create_mesher()
  {
    Mesher* m = new Mesher(cdt, criteria);
    m->set_seeds(seeds.begin(), seeds.end());
    return m;
  }

private:
  static const QString my_filters;
  Criteria criteria;
  Tr cdt;
  Mesher* mesher;
  Seeds seeds;

  Preferences* prefs;
  // QPopupMenu *pmCriteria;
  // int menu_id;

  CGAL::Qt_widget* widget;
  CGAL::Qt_widget_get_point<K>* get_point;
  CGAL::Qt_widget_get_point<K>* get_seed;
  CGAL::Qt_widget_get_polygon<CGALPolygon>* get_polygon;
  Follow_mouse* follow_mouse;

  typedef CGAL::Show_points<Tr,
                            Tr::Finite_vertices_iterator,
                            Vertex_to_point>
    Show_points_from_triangulation;

  typedef CGAL::Show_points<Seeds, Seeds::const_iterator>
    Show_seeds;

  typedef CGAL::Show_segments<Mesher,
                              Mesher::Encroached_edges_const_iterator,
                              Edge_to_segment>
    Show_encroached_edges;

  Show_points_from_triangulation* show_points;
  Show_seeds* show_seeds;
  Show_encroached_edges* show_encroached_edges;
  Meshing_debugging_layer* debug_layer;
  Show_bad_faces<Mesher>* show_bad_faces;
  CGAL::Qt_layer_show_triangulation<Tr>* show_triangulation;
  CGAL::Qt_layer_show_triangulation_constraints<Tr>* show_constraints;
  CGAL::Qt_layer_show_circles<Tr>* show_circles;
  CGAL::Qt_widget_show_mouse_coordinates* show_coordinates;
  Show_in_domain_faces<Tr>* show_in_domain;

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
  QPushButton* pbShowCluster;
  QPushButton* pbShowEncroachedEdges;
  QPushButton* pbShowBadFaces;
  QToolBar *toolBarAdvanced;
  int timer_interval;
  int step_length;
};

const QString MyWindow::my_filters =
"Constrained edges (*.edg);;"
"Shewchuk Triangle .poly files (*.poly);;"
"All files (*)";

int main(int argc, char** argv)
{
  QApplication app( argc, argv );
  MyWindow* W = new MyWindow();
  app.setMainWidget(W);
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W->setIcon(cgal_icon);
#endif
  W->show();

  if( argc == 2 )
    W->openTriangulation(QString(argv[1]));

  return app.exec();
}

// moc_source_file: mesh_2_demo.cpp
#include "mesh_2_demo.moc"

// moc_source_file: Show_clusters.h
#include "Show_clusters.moc"

