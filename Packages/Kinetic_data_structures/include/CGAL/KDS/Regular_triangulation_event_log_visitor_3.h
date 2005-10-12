// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_ELV_H
#define CGAL_KDS_REGULAR_TRIANGULATION_3_WATCHER_ELV_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Delaunay_triangulation_event_log_visitor_3.h>

CGAL_KDS_BEGIN_NAMESPACE

struct Regular_triangulation_event_log_visitor_3: 
  public Delaunay_triangulation_event_log_visitor_3  {
  typedef Delaunay_triangulation_event_log_visitor_3 P;
  Regular_triangulation_event_log_visitor_3(){}

  
  template <class Key, class Cell>
  void pre_move(Key k, Cell c){
    std::ostringstream out;
    out << "Moving " << k << " from ";
    write_cell(c, out);
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  template <class Key, class Cell>
  void post_move(Key k, Cell c){
    std::ostringstream out;
    out << "Moved " << k << " from ";
    write_cell(c, out);
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }
  
  template <class Key, class Cell>
  void pre_push(Key k, Cell c){
    std::ostringstream out;
    out << "Pushing " << k << " into ";
    write_cell(c, out);
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  template <class Vertex_handle>
  void post_push(Vertex_handle vh){
    std::ostringstream out;
    out << "Pushed " << vh->point();
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  template <class Vertex_handle>
  void pre_pop(Vertex_handle vh){
    std::ostringstream out;
    out << "Popping " << vh->point();
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }

  template <class Key, class Cell>
  void post_pop(Key k, Cell c){
    std::ostringstream out;
    out << "Popped " << k << " from ";
    write_cell(c, out);
    log_.push_back(out.str());
    CGAL_KDS_LOG(LOG_LOTS, "Logging: " << out.str() << std::endl);
  }


  template <class Cell>
  void write_cell(Cell c, std::ostream &out){
    out << "{";
    for (unsigned int i=0; i< 4; ++i){
      out << c->vertex(i)->point();
      if (i != 3){
	out << " ";
      }
    }
    out << "}";
  }
};

CGAL_KDS_END_NAMESPACE

#endif
