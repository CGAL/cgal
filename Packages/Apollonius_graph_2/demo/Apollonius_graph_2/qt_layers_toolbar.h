#ifndef QT_LAYERS_TOOLBAR_H
#define QT_LAYERS_TOOLBAR_H

#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qstring.h>
#include <qwhatsthis.h>

#include <CGAL/IO/Qt_layer_show_mouse_coordinates.h>

#include "qt_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/circle.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
//#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/notool.xpm>
//#include "removecircle.xpm"


class Layers_toolbar : public QObject
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *w, QMainWindow *mw, AG_2& ag)
    : nr_of_buttons(0) {

    showVD = new Voronoi_diagram_layer<AG_2>(ag);
    showDG = new Delaunay_graph_layer<AG_2>(ag);
    showNT = new Visible_sites_layer<AG_2>(ag);
    showTR = new Hidden_sites_layer<AG_2>(ag);

    // set the widget
    widget = w;
    window = mw;
    window->statusBar();

    widget->attach(showTR);
    widget->attach(showDG);
    widget->attach(showVD);
    widget->attach(showNT);

    maintoolbar = new QToolBar("tools", mw, QMainWindow::Top, TRUE, "Tools");

    but[0] = new QToolButton(QPixmap( (const char**)points_small_xpm ),
			     "Show weighted points", 
			     0, 
			     this, 
			     SLOT(show_sites()), 
			     maintoolbar, 
			     "Show weighted points");

    but[1] = new QToolButton(QPixmap( (const char**)voronoi_small_xpm ),
			     "Show Voronoi diagram", 
			     0, 
			     this, 
			     SLOT(show_apollonius_diagram()), 
			     maintoolbar, 
			     "Show Voronoi_diagram");

    
    but[2] = new QToolButton(QPixmap( (const char**)triangulation_small_xpm ),
			     "Show Delaunay graph", 
			     0, 
			     this, 
			     SLOT(show_apollonius_graph()), 
			     maintoolbar, 
			     "Show Delaunay graph");

    but[3] = new QToolButton(QPixmap( (const char**)point_small_xpm ),
			     "Insert point", 
			     0, 
			     this, 
			     SLOT(insert_point_mode()), 
			     maintoolbar, 
			     "Insert point");

    but[4] = new QToolButton(QPixmap( (const char**)circle_small_xpm ),
			     "Insert circle", 
			     0, 
			     this, 
			     SLOT(insert_circle_mode()), 
			     maintoolbar, 
			     "Insert circle");
    
#if 0
    but[5] = new QToolButton(QPixmap( (const char**)removecircle_xpm ),
			     "Remove weighted point", 
			     0, 
			     this, 
			     SLOT(remove_mode()), 
			     maintoolbar, 
			     "Remove weighted point");
#else
    but[5] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Remove weighted point", 
			     0, 
			     this, 
			     SLOT(remove_mode()), 
			     maintoolbar, 
			     "Remove weighted point");
#endif
    showDG->deactivate();

    nr_of_buttons = 6;
    for(int i = 0; i < nr_of_buttons; i++){
      but[i]->setToggleButton(TRUE);
    }

    but[0]->toggle();
    but[1]->toggle();
    //    but[2]->toggle();
    //    but[3]->toggle();
    but[4]->toggle();
    //    but[5]->toggle();
  }

  ~Layers_toolbar() {
    delete showVD;
    delete showDG;
    delete showTR;
    delete showNT;
  }

  inline QToolBar* toolbar() { return maintoolbar; };

signals:
  void new_object(CGAL::Object);
  void removeModeChanged(bool);
  void inputModeChanged(bool);
		
private slots:
  void show_sites() {
    if ( but[0]->isOn() ) {
      showNT->activate();
      showTR->activate();
    } else {
      showNT->deactivate();
      showTR->deactivate();
    }
    widget->redraw();
  }

  void show_apollonius_diagram() {
    if ( but[1]->isOn() ) {
      showVD->activate();
    } else {
      showVD->deactivate();
    }
    widget->redraw();
  }

  void show_apollonius_graph() {
    if ( but[2]->isOn() ) {
      //      widget->activate(showDG);
      showDG->activate();
    } else {
      //      widget->deactivate(showDG);
      showDG->deactivate();
    }
    widget->redraw();
  }

  void insert_point_mode() {
    if ( !but[3]->isOn() ) {
      but[3]->toggle();
      return;
    }

    but[4]->toggle();
    emit inputModeChanged( true );
  }

  void insert_circle_mode() {
    if ( !but[4]->isOn() ) {
      but[4]->toggle();
      return;
    }

    but[3]->toggle();
    emit inputModeChanged( false );
  }

  void remove_mode() {
    emit removeModeChanged( but[5]->isOn() );
  }

private:
  QToolBar		*maintoolbar;
  QToolButton		*but[10];
  CGAL::Qt_widget	*widget;
  QMainWindow		*window;	
  int			nr_of_buttons;

  Voronoi_diagram_layer<AG_2>             *showVD;
  Delaunay_graph_layer<AG_2>              *showDG;
  Visible_sites_layer<AG_2>               *showNT;
  Hidden_sites_layer<AG_2>                *showTR;

};//end class

#endif // QT_LAYERS_TOOLBAR_H
