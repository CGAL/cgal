// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// $URL$
// $Id$
//
// Authors : Monique Teillaud, Sylvain Pion, Pedro Machado, Radu Ursu
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (CGAL - Effective Computational Geometry for Curves and Surfaces)

#include <CGAL/basic.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>

#include <CGAL/IO/Qt_widget_circular_arc_2.h>
#include <CGAL/IO/Qt_widget_circular_arc_endpoint_2.h>

#include "Qt_widget_get_segment.h"
#include "Qt_widget_get_arc.h"
#include "sweeper.xpm"
#include "trash.xpm"
#include "get_arc.xpm"
#include "lines_icon.xpm"

#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Object.h>

typedef CGAL::Quotient<CGAL::MP_Float>                       NT;


typedef CGAL::Cartesian<NT>                                 Linear_k;

typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>          Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>       Circular_k;

typedef Circular_k::Line_arc_2                              Line_arc_2;
typedef Circular_k::Segment_2                               Segment;
typedef Circular_k::Circular_arc_2                          Circular_arc_2;
typedef Circular_k::Circular_arc_point_2                    Circular_arc_point_2;

typedef std::vector<CGAL::Object>                           ArcContainer;

const QString my_title_string("CGAL :: "
                              "Intersecting Circles, Line Arcs, Circular Arcs");

class Qt_layer_show_intersections
  : public CGAL::Qt_widget_layer
{
  ArcContainer _ac;
  bool show_intersections;

public:

  Qt_layer_show_intersections()
    : _ac(), show_intersections(true) {}

  void swap_show() { show_intersections = ! show_intersections; }

  const ArcContainer & arc_container() const { return _ac; }
        ArcContainer & arc_container()       { return _ac; }

  void draw()
  {
    if (!show_intersections)
      return;
    *widget << CGAL::RED;
    int l = _ac.size();

    for(int i=0; i<l; i++) {
      for(int j=i+1; j<l; j++) {
        std::vector< CGAL::Object > res;
        Circular_arc_2 ca1, ca2;
        Line_arc_2 la1, la2;
        if(assign(ca1, _ac[i])) {
          if(assign(ca2, _ac[j])) {
            CGAL::intersection(ca1, ca2, std::back_inserter(res));
          } else {
            CGAL::intersection(ca1, la2, std::back_inserter(res));
          }
        } else {
          if(assign(ca2, _ac[j])) {
            CGAL::intersection(la1, ca2, std::back_inserter(res));
          } else {
            CGAL::intersection(la1, la2, std::back_inserter(res));
          }
        }
        for(unsigned k=0; k<res.size(); k++) {
          std::pair<Circular_arc_point_2, unsigned> pair;
          Circular_arc_2 ca;
          Line_arc_2 la;
          if(assign(pair, res[k])) *widget << pair.first;
          else if(assign(ca, res[k])) *widget << ca;
          else if(assign(la, res[k])) {
            *widget << Segment(Circular_k::Point_2(to_double(la.source().x()),
                                                   to_double(la.source().y())),
                               Circular_k::Point_2(to_double(la.target().x()),
                                                   to_double(la.target().y())));
          }
        }
      }
    }
  }
};

struct Qt_layer_show_ch
  : public CGAL::Qt_widget_layer
{
  ArcContainer _ac;
  bool show_arcs;

public:

  Qt_layer_show_ch()
    : show_arcs(true) {}

  const ArcContainer & arc_container() const { return _ac; }

  ArcContainer & arc_container() { return _ac; }

  void swap_show() { show_arcs = ! show_arcs; }

  void draw()
  {
    if (!show_arcs)
      return;
    *widget << CGAL::BLUE;
    int l = _ac.size();
    for (int i=0; i<l; i++){
      Circular_arc_2 ca;
      Line_arc_2 la;
      if(assign(ca, _ac[i])) *widget << ca;
      else if(assign(la, _ac[i])) {
        *widget << Segment(Circular_k::Point_2(to_double(la.source().x()),
                                               to_double(la.source().y())),
                           Circular_k::Point_2(to_double(la.target().x()),
                                               to_double(la.target().y())));
      }
    }
  }
};

class MyWindow
  : public QMainWindow
{
  Q_OBJECT

public:

  MyWindow(int w, int h)
    : something_changed(true)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()), this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");

    QToolBar * layers_toolbar = new QToolBar("Tools", this,
                                             QMainWindow::Top, TRUE, "Tools");

    QToolButton * show_intersections_button =
                   new QToolButton(QPixmap((const char**)::sweeper_xpm ),
                                   "Showing Intersections",
                                   0,
                                   this,
                                   SLOT(show_intersections()),
                                   layers_toolbar,
                                   "Showing Intersections");

    widget->attach(&show_intersections_layer);
    connect(show_intersections_button, SIGNAL(stateChanged(int)),
            &show_intersections_layer, SLOT(stateChanged(int)));

    QToolButton * show_container_button =
                   new QToolButton(QPixmap((const char**)::get_arc),
                                   "Show Original arcs",
                                   0,
                                   this,
                                   SLOT(show_original_arcs()),
                                   layers_toolbar,
                                   "Show Original arcs");
    show_container_button->setToggleButton(true);

    show_container_button->toggle();
    connect(show_container_button, SIGNAL(stateChanged(int)),
            &testlayer, SLOT(stateChanged(int)));

    QToolButton * line_circle_button =
                       new QToolButton(QPixmap((const char**)::lines_icon),
                                       "Use Line",
                                       0,
                                       this,
                                       SLOT(change_line_circle()),
                                       layers_toolbar,
                                       "Use Circular arc");

    line_circle_button->setToggleButton(true);
    arc_circle = true;

    QToolButton * clear_button =
                       new QToolButton(QPixmap((const char**)trash),
                                       "Clear",
                                       0,
                                       this,
                                       SLOT(clear_container()),
                                       layers_toolbar,
                                       "Clear");

    connect(clear_button, SIGNAL(stateChanged(int)),
            &testlayer, SLOT(stateChanged(int)));

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor(CGAL::WHITE);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);

    // layers
    widget->attach(&testlayer);
    get_arc_layer = new CGAL::Qt_widget_get_arc<Circular_k>;
    get_segment_layer =  new CGAL::Qt_widget_get_segment<Circular_k>;
    widget->attach(get_arc_layer);
    connect(get_arc_layer, SIGNAL(new_object_time()), this, SLOT(get_arc()));
    connect(get_segment_layer, SIGNAL(new_object_time()),
	    this, SLOT(get_arc()));
  }

public slots:

  void new_instance()
  {
    widget->lock();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
    something_changed = true;
  }

  void get_arc()
  {
    CGAL::Object o;
    if (arc_circle) {
      Circular_arc_2 ca = get_arc_layer->get_circular_arc();
      o = make_object(ca);
    }
    else {
      Line_arc_2 la = get_segment_layer->get_line_arc();
      o = make_object(la);
    }
    arc_container().push_back(o);
    intersections_container().push_back(o);
    something_changed = true;
    widget->redraw();
  }

  void show_original_arcs()
  {
    testlayer.swap_show();
    something_changed = true;
    widget->redraw();
  }

  void show_intersections()
  {
    show_intersections_layer.swap_show();
    something_changed = true;
    widget->redraw();
  }

  void change_line_circle(){
    arc_circle = ! arc_circle;
    if(arc_circle){
      widget->detach(get_segment_layer);
      widget->attach(get_arc_layer);
    }
    else{
       widget->detach(get_arc_layer);
       widget->attach(get_segment_layer);
    }
    something_changed = true;
    widget->redraw();
  }

  void clear_container()
  {
    arc_container().clear();
    intersections_container().clear();
    something_changed = true;
    widget->redraw();
  }

private slots:

  void about()
  {
    QMessageBox::about(this, my_title_string, 
    "This is a demo of the CGAL's Circular_kernel_2.\n Particularly, the intersection functionality.");
  }

  void aboutQt()
  {
    QMessageBox::aboutQt(this, my_title_string);
  }

  void howto()
  {
    CGAL::Qt_help_window *help =
          new CGAL::Qt_help_window("help/index.html", ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void new_window()
  {
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->widget->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed = true;
  }

  void timer_done()
  {
    if (something_changed)
      widget->redraw();
    something_changed = false;
  }

private:
  bool arc_circle;

  const ArcContainer & arc_container() const
  { return testlayer.arc_container(); }

  ArcContainer & arc_container()
  { return testlayer.arc_container(); }

  const ArcContainer & intersections_container() const
  { return show_intersections_layer.arc_container(); }

  ArcContainer & intersections_container()
  { return show_intersections_layer.arc_container(); }

  CGAL::Qt_widget                         *widget;
  CGAL::Qt_widget_standard_toolbar        *stoolbar;
  bool                                    something_changed;
  Qt_layer_show_ch                        testlayer;
  Qt_layer_show_intersections             show_intersections_layer;
  CGAL::Qt_widget_get_segment<Circular_k> *get_segment_layer;
  CGAL::Qt_widget_get_arc<Circular_k>     *get_arc_layer;
};

#include "demo.moc"

int
main(int argc, char **argv)
{
  QApplication app(argc, argv);
  MyWindow widget(800, 700);
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
#endif
  widget.show();
  return app.exec();
}

