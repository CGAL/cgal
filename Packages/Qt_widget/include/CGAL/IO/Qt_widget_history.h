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
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_WIDGET_HISTORY_H
#define CGAL_QT_WIDGET_HISTORY_H

#include <vector>

namespace CGAL {
  class History_atom{
  public:
    History_atom(double xmin, double xmax, double ymin, double ymax)
      : m_x1(xmin), m_x2(xmax), m_y1(xmin), m_y2(xmax){};
    double x1(){return m_x1;}
    double x2(){return m_x2;}
    double y1(){return m_y1;}
    double y2(){return m_y2;}
  private:
    double  m_x1,
            m_x2,
            m_y1,
            m_y2;
  };

  class Qt_widget_history{
    public:
      Qt_widget_history() : nr_of_items(0), current_item(0) {};
      ~Qt_widget_history(){};

    bool back(){
      if(nr_of_items != 0){
        current_item--;
        return true;
      }
      return false;
    }
    bool forward(){
      if(current_item < nr_of_items){
        current_item++;
        return true;
      }
      return true;
    }
    History_atom* get_atom(){
      return &history_vector[current_item];
    }
    void add_to_history(double xmin, double xmax, double ymin, double ymax)
      {
        History_atom atom(xmin, xmax, ymin, ymax);
        history_vector.push_back(atom);
        nr_of_items++;
        current_item = history_vector.size();
      }
      
    private:
      
      std::vector<History_atom> history_vector;
      int current_item;
      int nr_of_items;
    protected:
  };
  
} // namespace CGAL end

#endif // CGAL_QT_WIDGET_HISTORY_H
