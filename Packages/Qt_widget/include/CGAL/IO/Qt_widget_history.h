// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu & Laurent Rineau

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
