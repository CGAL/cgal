// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
// ----------------------------------------------------------------------
//
// file          : polygon_2_toolbar.C
// package       : Qt_widget
// author(s)     : Radu Ursu
// coordinator   : Laurent Rineau
//
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifdef CGAL_USE_QT

#include "polygon_2_toolbar.h"

// icons

#include <CGAL/IO/pixmaps/arrow.xpm>
#include <CGAL/IO/pixmaps/polygon.xpm>
#include <CGAL/IO/pixmaps/point.xpm>

#include <qiconset.h>

/* XPM */
static char *polygon2_xpm[] = {
/* columns rows colors chars-per-pixel */
"16 16 5 1",
"  c black",
". c #1F0670",
"X c #EDA1A1",
"o c gray100",
"O c None",
/* pixels */
"OOOOXXXXXXXXOOOO",
"OOOXX    XXXXOOO",
"OOXXX XX XXXXXOO",
"OXXXX XX XXXXXXO",
"XXXXX       XXXX",
"XXXXXXXX XX XXXX",
"XXXXXXXX XX XXXX",
"XXXXXXXX XX XXXX",
"XXXXXXXX       X",
"XX.....XXXX XX X",
"XXX....XXXX XX X",
"XXX....XXXX    X",
"OX.....XXXXXXXXO",
"O...XX.XXXXXXXOO",
"...XXXXXXXXXXOOO",
"..OOXXXXXXXXOOOO"
};

Polygon_toolbar::Polygon_toolbar(CGAL::Qt_widget *w, QMainWindow *mw) :
  QToolBar(mw, "NT")
{
    //when it is created, the toolbar has 0 buttons
    nr_of_buttons = 0;
    //set the widget
    widget = w;
    widget->attach(&getsimplepoly);
    widget->attach(&getpoly);
    widget->attach(&getpoint);
    getsimplepoly.deactivate();
    getpoly.deactivate();
    getpoint.deactivate();

    QIconSet set0(QPixmap( (const char**)arrow_small_xpm ),
                  QPixmap( (const char**)arrow_xpm ));
    QIconSet set1(QPixmap( (const char**)polygon_small_xpm ),
                  QPixmap( (const char**)polygon_xpm ));
    QIconSet set2(QPixmap( (const char**)point_small_xpm ),
                  QPixmap( (const char**)point_xpm ));
    QIconSet set3(QPixmap( (const char**)polygon2_xpm ),
                  QPixmap( (const char**)polygon2_xpm ));		
  but[0] = new QToolButton(this, "deactivate layer");
  but[0]->setIconSet(set0);
  but[0]->setTextLabel("Deactivate Layer");
  
  but[1] = new QToolButton(this, "spolygon");
  but[1]->setIconSet(set1);
  but[1]->setTextLabel("Input Simple Polygon");

  but[2] = new QToolButton(this, "polygon");
  but[2]->setIconSet(set3);
  but[2]->setTextLabel("Input Polygon");

  but[3] = new QToolButton(this, "point");
  but[3]->setIconSet(set2);
  but[3]->setTextLabel("Input Point");
  
  button_group = new QButtonGroup(0, "exclusive_group");
  button_group->insert(but[0]);
  button_group->insert(but[1]);
  button_group->insert(but[2]);
  button_group->insert(but[3]);
  button_group->setExclusive(true);

  but[0]->setToggleButton(true);
  but[1]->setToggleButton(true);
  but[2]->setToggleButton(true);
  but[3]->setToggleButton(true);
  
  connect(but[1], SIGNAL(stateChanged(int)),
        &getsimplepoly, SLOT(stateChanged(int)));
  connect(but[2], SIGNAL(stateChanged(int)),
        &getpoly, SLOT(stateChanged(int)));
  connect(but[3], SIGNAL(stateChanged(int)),
        &getpoint, SLOT(stateChanged(int)));
  nr_of_buttons = 4;
  };

#include "polygon_2_toolbar.moc"

#endif
