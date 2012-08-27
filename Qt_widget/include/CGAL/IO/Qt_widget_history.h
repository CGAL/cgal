// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Radu Ursu and Laurent Rineau

#ifndef CGAL_QT_WIDGET_HISTORY_H
#define CGAL_QT_WIDGET_HISTORY_H

#include <list>
#include <algorithm>
#include <CGAL/IO/Qt_widget.h>

namespace CGAL {

  class History_atom {
  public:
    History_atom() {};
    virtual ~History_atom() {};

    void save(const Qt_widget& widget){
      xmin = widget.x_min();
      ymin = widget.y_min();
      xmax = widget.x_max();
      ymax = widget.y_max();
    }
    
    void restore(Qt_widget& widget) const {
      widget.set_window(xmin, xmax, ymin, ymax);
    }
  private:
    double xmin, xmax, ymin, ymax;
  };

  class Qt_widget_history : public QObject {
    Q_OBJECT
  public:
    Qt_widget_history(Qt_widget* parent, const char* name = 0 );

  signals:
    void backwardAvaillable(bool);
    void forwardAvaillable(bool);

  public slots:
    void backward();
    void forward();
    
  private:
    struct Free {
      void operator()(History_atom* atom) const
      {
        delete atom;
      }
    };    
    
  public slots:
    void save();
    void clear() {
      std::for_each(history_list.begin(), history_list.end(), Free());
      history_list.clear();
      it = history_list.begin();
      emit backwardAvaillable(false);
      emit forwardAvaillable(false);
    }

  private:
    void restore(){
      disconnect( widget, SIGNAL(rangesChanged()), 
        this, SLOT(save()));
      (*it)->restore(*widget);
      connect(widget, SIGNAL(rangesChanged()), 
        this, SLOT(save()));
      widget->redraw();
    }

  private:
    std::list<History_atom*> history_list;
    std::list<History_atom*>::iterator it;
    Qt_widget* widget;
  };

} // namespace CGAL end

#endif // CGAL_QT_WIDGET_HISTORY_H
