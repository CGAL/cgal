// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// file          : src/CGALQt/Qt_widget_standard_toolbar.C
// package       : Qt_widget (1.2.46)
// maintainer    : Laurent Rineau <rineau@clipper.ens.fr>
// author(s)     : Radu Ursu
// release       : $CGAL_Revision: CGAL-2.5-I-35 $
// release_date  : $CGAL_Date: 2002/10/18 $
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ======================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>

// icons
#include <CGAL/IO/pixmaps/zoom_in_rect.xpm>
#include <CGAL/IO/pixmaps/zoom_out.xpm>
#include <CGAL/IO/pixmaps/zoom_in.xpm>
#include <CGAL/IO/pixmaps/focus.xpm>
#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/back.xpm>
#include <CGAL/IO/pixmaps/forward.xpm>
#include <CGAL/IO/pixmaps/mouse_coord.xpm>

#include <qbuttongroup.h>
#include <qiconset.h>
#include <qmainwindow.h>
#include <qtoolbutton.h>
#include <qstring.h>
#include <qvariant.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_history.h>
#include <CGAL/IO/Qt_widget_focus.h>
#include <CGAL/IO/Qt_widget_zoomrect.h>
#include <CGAL/IO/Qt_widget_handtool.h>
#include <CGAL/IO/Qt_widget_show_mouse_coordinates.h>

namespace CGAL {
  Qt_widget_standard_toolbar::
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *parent,
			     const char* name) :
    QToolBar(parent, name),
    widget(w)
  {
    setLabel("Qt_widget standard toolbar");
    fill_toolbar(parent);
  };
  
  Qt_widget_standard_toolbar::
  Qt_widget_standard_toolbar(Qt_widget *w, QMainWindow *mw,
			     QWidget* parent,
			     bool newLine,
			     const char* name) :
    QToolBar("Qt_widget standard toolbar", mw, parent, newLine, name),
    widget(w)
  {
    fill_toolbar(mw);
  };
  
  void Qt_widget_standard_toolbar::fill_toolbar(QMainWindow *mw)
  {
    Qt_widget_focus* focuslayer = 
      new Qt_widget_focus(this, "focuslayer");
    Qt_widget_zoomrect* zoomrectlayer = 
      new Qt_widget_zoomrect(this, "zoomrectlayer");
    Qt_widget_handtool* handtoollayer = 
      new Qt_widget_handtool(this, "handtoollayer");
    Qt_widget_show_mouse_coordinates* showcoordlayer = 0; // created below

    widget->attach_standard(focuslayer);
    widget->attach_standard(zoomrectlayer);
    widget->attach_standard(handtoollayer);
    focuslayer->deactivate();
    zoomrectlayer->deactivate();
    handtoollayer->deactivate();

    if (mw)
      {
	mw->statusBar();
	showcoordlayer = 
	  new Qt_widget_show_mouse_coordinates(*mw, this,
					       "showcoordlayer");
	widget->attach_standard(showcoordlayer);
	showcoordlayer->does_eat_events = false;
	widget->setMouseTracking(true);
      }

    const QIconSet 
      back_pixmap(QPixmap(static_cast<const char **>(back_small_xpm) ),
		  QPixmap(static_cast<const char **>(back_xpm) ));

    const QIconSet 
      forward_pixmap(QPixmap(static_cast<const char **>(forward_small_xpm) ),
		     QPixmap(static_cast<const char **>(forward_xpm) ));

    QToolButton* backBt = new QToolButton(this, "History Back");
    backBt->setIconSet(back_pixmap);
    backBt->setTextLabel("History Back");

    QToolButton* forwardBt = new QToolButton(this, "History Forward");
    forwardBt->setIconSet(forward_pixmap);
    forwardBt->setTextLabel("History Forward");
  
    addSeparator();

    QToolButton* zoominBt =
      new QToolButton(QPixmap(zoomin_xpm ),
			       "Zoom in", 
			       0, 
			       this, 
			       SLOT(zoomin()), 
			       this, 
			       "Zoom in");
    zoominBt->setTextLabel("Scaling factor X2");

    QToolButton* zoomoutBt = 
      new QToolButton(QPixmap(zoomout_xpm ),
				"Zoom out", 
				0, 
				this, 
				SLOT(zoomout()), 
				this, 
				"Zoom out");
    zoomoutBt->setTextLabel("Scaling factor 1/2");

    addSeparator();

    QToolButton* zoomrectBt = new QToolButton(this, "focus on region");
    zoomrectBt->setPixmap(QPixmap(zoomin_rect_xpm ));
    zoomrectBt->setTextLabel("Focus on region");

    QToolButton* focusBt = new QToolButton(this, "focus");
    focusBt->setPixmap(QPixmap(focus_xpm ));
    focusBt->setTextLabel("Focus on point");

    QToolButton* handtoolBt = new QToolButton(this, "handtool");
    handtoolBt->setPixmap(QPixmap(hand_xpm ));
    handtoolBt->setTextLabel("Pan tool");

    addSeparator();

    QToolButton* showcoordBt = new QToolButton(this, "mouse");
    showcoordBt->setPixmap(QPixmap(mouse_coord_xpm) );
    showcoordBt->setTextLabel("Mouse Coordinates");

    button_group = new QButtonGroup(0, "My_group");
    // this button has no parent and is destroyed manually in the
    // destructor

    // below is the list of buttons in the group
    QToolButton* const button_group_list[] = { zoomrectBt,
					       focusBt,
					       handtoolBt };
    for(int i=0; i<3; ++i)
      {
	button_group_list[i]->setToggleButton(true);
	button_group->insert(button_group_list[i]);
      }
    button_group->setExclusive(false);
    connect(button_group, SIGNAL(clicked(int)),
	    this, SLOT(group_clicked(int)));

    showcoordBt->setToggleButton(true);
    showcoordBt->toggle();

    connect(zoomrectBt, SIGNAL(stateChanged(int)),
	    zoomrectlayer, SLOT(stateChanged(int)));
    connect(focusBt, SIGNAL(stateChanged(int)),
	    focuslayer, SLOT(stateChanged(int)));	
    connect(handtoolBt, SIGNAL(stateChanged(int)),
	    handtoollayer, SLOT(stateChanged(int)));
    connect(showcoordBt, SIGNAL(stateChanged(int)),
	    showcoordlayer, SLOT(stateChanged(int)));


    // history setting
    history = new Qt_widget_history(widget, "standard history");
    connect(backBt, SIGNAL(clicked()),
	    history, SLOT(backward()));
    connect(forwardBt, SIGNAL(clicked()),
	    history, SLOT(forward()));

    connect(history, SIGNAL(backwardAvaillable(bool)),
            backBt, SLOT(setEnabled(bool)));
    connect(history, SIGNAL(forwardAvaillable(bool)),
            forwardBt, SLOT(setEnabled(bool)));
    history->clear();
  };

  void Qt_widget_standard_toolbar::group_clicked(int i)
  {
    QButton* bt = button_group->find(i);
    
    if( bt->isOn() )
      for(int n = 0; n < button_group->count(); n++)
	if( n != i )
	  button_group->find(n)->setProperty("on", false);
  }

  void Qt_widget_standard_toolbar::zoomin()
  {
    widget->zoom(2); //
  }
  void Qt_widget_standard_toolbar::zoomout()
  {
    widget->zoom(0.5);
  }
  void Qt_widget_standard_toolbar::back()
  {
    history->backward();
  }
  void Qt_widget_standard_toolbar::forward()
  {
    history->forward();
  }
  void Qt_widget_standard_toolbar::clear_history()
  {
    history->clear();
  }

}//end namespace
#include "Qt_widget_standard_toolbar.moc"

#endif
