#ifndef QT_LAYERS_TOOLBAR_H
#define QT_LAYERS_TOOLBAR_H

#include <CGAL/IO/Qt_widget.h>

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qstatusbar.h>
#include <qstring.h>
#include <qwhatsthis.h>

#include "qt_layers.h"

// icons
#include <CGAL/IO/pixmaps/points.xpm>
#include <CGAL/IO/pixmaps/point.xpm>
#include <CGAL/IO/pixmaps/line.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include <CGAL/IO/pixmaps/triangulation.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>
#include <CGAL/IO/pixmaps/notool.xpm>
#include <CGAL/IO/pixmaps/holddown.xpm>
//#include "remove.xpm"

typedef enum { SVD_POINT, SVD_SEGMENT, SVD_POLYGON } Input_mode;

class Layers_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Layers_toolbar(CGAL::Qt_widget *widget, SVD_2& svd,
		 const QString& label, QMainWindow* mainWindow,
		 QWidget* parent, bool newLine = FALSE,
		 const char* name = 0, WFlags f = 0 )
    : QToolBar(label, mainWindow, parent, newLine, name, f),
      nr_of_buttons(0), input_mode(SVD_SEGMENT)
  {
    showVD = new Voronoi_diagram_layer<SVD_2>(svd);
    showSI = new Sites_layer<SVD_2>(svd);
    showSK = new Skeleton_layer<SVD_2>(svd);

    // set the widget
    this->widget = widget;
    window = mainWindow;
    window->statusBar();

    widget->attach(showVD);
    widget->attach(showSK);
    widget->attach(showSI);


    but[0] = new QToolButton(QPixmap( (const char**)points_xpm ),
			     "Show sites", 
			     0, 
			     this, 
			     SLOT(show_sites()), 
			     this, 
			     "Show sites");

    but[1] = new QToolButton(QPixmap( (const char**)voronoi_xpm ),
			     "Show Voronoi diagram", 
			     0, 
			     this, 
			     SLOT(show_voronoi()), 
			     this, 
			     "Show Voronoi diagram");

    but[2] = new QToolButton(QPixmap( (const char**)triangulation_xpm ),
			     "Show skeleton", 
			     0, 
			     this, 
			     SLOT(show_skeleton()), 
			     this, 
			     "Show skeleton");
    
    but[3] = new QToolButton(QPixmap( (const char**)point_xpm ),
			     "Input points", 
			     0, 
			     this, 
			     SLOT(input_points_mode()), 
			     this, 
			     "Input points");

    
    but[4] = new QToolButton(QPixmap( (const char**)line_xpm ),
			     "Input segments", 
			     0, 
			     this, 
			     SLOT(input_segments_mode()), 
			     this, 
			     "Input segments");
    but[5] = new QToolButton(QPixmap( (const char**)polygon_xpm ),
			     "Input polygons", 
			     0, 
			     this, 
			     SLOT(input_polygons_mode()), 
			     this, 
			     "Input polygons");
#if 1
    but[6] = new QToolButton(QPixmap( (const char**)notool_xpm ),
			     "Remove site", 
			     0, 
			     this, 
			     SLOT(remove_mode()), 
			     this, 
			     "Remove site");
#else
    but[6] = new QToolButton(QPixmap( (const char**)remove_xpm ),
			     "Remove site", 
			     0, 
			     this, 
			     SLOT(remove_mode()), 
			     this, 
			     "Remove site");
#endif
    but[7] = new QToolButton(QPixmap( (const char**)holddown_xpm ),
			     "Snap mode", 
			     0, 
			     this, 
			     SLOT(snap_mode()), 
			     this, 
			     "Snap mode");

    showSK->deactivate();

    nr_of_buttons = 8;
    for(int i = 0; i < nr_of_buttons; i++){
      but[i]->setToggleButton(TRUE);
    }

    but[0]->toggle();
    but[1]->toggle();
    //    but[2]->toggle();
    //    but[3]->toggle();
    but[4]->toggle();
    //    but[5]->toggle();
    //    but[6]->toggle();
    //    but[7]->toggle();
  }

  ~Layers_toolbar() {
    delete showVD;
    delete showSI;
    delete showSK;
  }

  inline QToolBar* toolbar() { return this; };

signals:
  void new_object(CGAL::Object);
  void inputModeChanged(Input_mode);
  void insertModeChanged(bool);
  void snapModeChanged(bool);
  
private slots:
  void show_sites() {
    if ( but[0]->isOn() ) {
      showSI->activate();
    } else {
      showSI->deactivate();
    }
    widget->redraw();
  }

  void show_voronoi() {
    if ( !but[1]->isOn() ) {
      showVD->deactivate();
    } else {
      showVD->activate();
      if ( but[2]->isOn() ) {
	but[2]->toggle();
	showSK->deactivate();
      }
    }

    widget->redraw();
  }

  void show_skeleton() {
    if ( !but[2]->isOn() ) { 
      showSK->deactivate();
    } else {
      showSK->activate();
      if ( but[1]->isOn() ) {
	but[1]->toggle();
	showVD->deactivate();
      }
    }

    widget->redraw();
  }

  void input_points_mode() {
    if ( !but[3]->isOn() ) {
      but[3]->toggle();
      return;
    }
    if ( input_mode == SVD_SEGMENT ) { but[4]->toggle(); }
    if ( input_mode == SVD_POLYGON ) { but[5]->toggle(); }

    input_mode = SVD_POINT;
    emit inputModeChanged( SVD_POINT );
  }

  void input_segments_mode() {
    if ( !but[4]->isOn() ) {
      but[4]->toggle();
      return;
    }
    if ( input_mode == SVD_POINT ) { but[3]->toggle(); }
    if ( input_mode == SVD_POLYGON ) { but[5]->toggle(); }

    input_mode = SVD_SEGMENT;
    emit inputModeChanged( SVD_SEGMENT );
  }

  void input_polygons_mode() {
    if ( !but[5]->isOn() ) {
      but[5]->toggle();
      return;
    }
    if ( input_mode == SVD_POINT ) { but[3]->toggle(); }
    if ( input_mode == SVD_SEGMENT ) { but[4]->toggle(); }

    input_mode = SVD_POLYGON;
    emit inputModeChanged( SVD_POLYGON );
  }

  void remove_mode() {
    emit insertModeChanged( but[6]->isOn() );
  }

  void snap_mode() {
    emit snapModeChanged( but[7]->isOn() );
  }

private:
  QToolButton		*but[15];
  CGAL::Qt_widget	*widget;
  QMainWindow		*window;	
  int			nr_of_buttons;
  Input_mode            input_mode;

  Voronoi_diagram_layer<SVD_2>             *showVD;
  Skeleton_layer<SVD_2>                    *showSK;
  Sites_layer<SVD_2>                       *showSI;

};//end class

#endif // QT_LAYERS_TOOLBAR_H
