// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_history.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

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
      std::cerr << "save (" << xmin << ", " << xmax 
		<< ", " << ymin << ", " << ymax << ")" << std::endl;
    }
    
    void restore(Qt_widget& widget) const {
      widget.set_window(xmin, xmax, ymin, ymax);
      std::cerr << "restore (" << xmin << ", " << xmax 
		<< ", " << ymin << ", " << ymax << ")" << std::endl;
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
      emit(backwardAvaillable(false));
      emit(forwardAvaillable(false));
    }

  private:

    void restore()
      {
	disconnect(widget, SIGNAL(rangesChanged()),
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
