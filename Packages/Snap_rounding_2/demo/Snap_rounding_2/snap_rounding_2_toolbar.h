#ifndef CGAL_SNAP_ROUNDING_2_TOOLBAR_H
#define CGAL_SNAP_ROUNDING_2_TOOLBAR_H

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include "segment_input_layer_with_snapping.h"

#include <qobject.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qbuttongroup.h>
#include <qmainwindow.h>

#include <list>

class Tools_toolbar : public QToolBar
{
  Q_OBJECT
public:
  Tools_toolbar(CGAL::Qt_widget * w, QMainWindow * mw, 
    std::list<Segment_2> * l1);
  ~Tools_toolbar(){};

signals:
  void new_object(CGAL::Object);

private:
  QToolButton * but[10];
  QButtonGroup * button_group;
  CGAL::Qt_widget * widget;
  int nr_of_buttons;

  Segment_input_layer<Rep>         segment_layer;
};

#endif
