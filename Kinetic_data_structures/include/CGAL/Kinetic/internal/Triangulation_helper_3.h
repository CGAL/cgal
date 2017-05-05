// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_TRIANGULATION_HELPER_H
#define CGAL_KINETIC_TRIANGULATION_HELPER_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/internal/triangulation_helpers_3.h>
namespace CGAL { namespace Kinetic { namespace internal {

template <class Tri>
class Triangulation_helper_3: public Tri
{

public:
  typedef typename Tri::Edge Edge;
  typedef typename Tri::Facet Facet;
  typedef typename Tri::Cell_handle Cell_handle;
  typedef typename Tri::Vertex_handle Vertex_handle;
  typedef typename Tri::Cell::Edge_label Edge_label;
  typedef typename Tri::Cell::Facet_label Facet_label;
  typedef typename Tri::Cell_circulator Cell_circulator;

  using Tri::incident_cells;

  Triangulation_helper_3(typename Tri::Geom_traits gt= Tri::Geom_traits()): Tri(gt){}

  bool has_degree_3(const Edge &e) const
  {
    return CGAL::Kinetic::internal::has_degree_3(*this, e);
  };

  bool hd3(const Facet &f) const
  {
    for (unsigned int i=0; i<3 ;++i) {
      bool d3 = has_degree_3(edge(f, i));
      CGAL_assertion(has_degree_3(opposite(edge(f,i)))== d3);
      if (d3) return true;
    }
    return false;
  }

  bool has_degree_4(const Vertex_handle vh) const
  {
    return CGAL::Kinetic::internal::has_degree_4(*this, vh);
  }

  bool has_degree_3(const Facet &f) const
  {
    CGAL_precondition_code(Edge e0= edge(f,0));
    CGAL_precondition_code(Edge e1= edge(f,1));
    CGAL_precondition_code(Edge e2= edge(f,2));
    CGAL_precondition(vertex(e0,1)== vertex(e1,0));
    CGAL_precondition(vertex(e1,1)== vertex(e2,0));
    CGAL_precondition(vertex(e2,1)== vertex(e0,0));
    CGAL_precondition(vertex(e0,0) != f.first->vertex(f.second));
    CGAL_precondition(vertex(e1,0) != f.first->vertex(f.second));
    CGAL_precondition(vertex(e2,0) != f.first->vertex(f.second));
    CGAL_precondition(hd3(f)== hd3(opposite(f)));
    return hd3(f);
  }

  Edge edge(const Facet &f, unsigned int i) const
  {
    CGAL_precondition( i <3);
    int f0= (i+1)%3;
    int f1= (i+2)%3;
    int i0= (f.second+1+f0)%4;
    int i1= (f.second+1+f1)%4;
    if (i1== f.second) {
      i1= (i1+1)%4;
    }
    Edge ret(f.first, i0, i1);
    CGAL_postcondition(i0 != f.second);
    CGAL_postcondition(i1 != f.second);
    return ret;
  }

  Vertex_handle other_vertex(const Facet &f, const Edge &e) const
  {
    return CGAL::Kinetic::internal::other_vertex(*this, f, e);
  }

  Facet_label label(const Facet &f) const
  {
    CGAL_assertion(f.first->facet_label(f.second) == opposite(f).first->facet_label(opposite(f).second));
    return f.first->facet_label(f.second);
  }
  Edge_label label(const Edge &f) const
  {
    Edge_label ret= f.first->edge_label(f.second, f.third);

    CGAL_postcondition_code(Cell_circulator cc= incident_cells(f));
    CGAL_postcondition_code(Cell_circulator ce= cc);
    CGAL_postcondition_code(do {
			      )
			    CGAL_postcondition_code( Edge ec= edge_in_cell(f, cc));
			    CGAL_postcondition(ec.first->edge_label(ec.second, ec.third)==ret);
			    CGAL_postcondition_code(++cc);
			    CGAL_postcondition_code(
						    } while (cc != ce));
    return ret;
  }
  void set_oriented_label(const Facet &f, Facet_label l) const
  {
    f.first->set_facet_label(f.second, l);

  }
  void set_label(const Facet &f, Facet_label l) const
  {
    set_oriented_label(f, l);
    set_oriented_label(opposite(f), l);
    CGAL_postcondition(label(f)==l);
  }

  Edge cross(Edge &e) const
  {
    int a=-1, b=-1;
    for (int i=0; i<4; ++i) {
      if (i != e.second && i != e.third) {
	if (a==-1) {
	  a=i;
	}
	else {
	  b=i;
	  break;
	}
      }
    }
    return Edge(e.first, a,b);
  }

  Edge edge_around_vertex(const Cell_handle cell, const Vertex_handle vh, unsigned int i) const
  {
    CGAL_assertion(i<3);
    int vi= cell->index(vh);
    Edge ret(cell, vi, (vi+i+1)%4);
    //CGAL_assertion(is_edge(ret.first, ret.second, ret.third));
    return ret;
  }

  //! Get a facet around the vh
  /*!
    The facet i should be opposite edge i.
  */
  Facet facet_around_vertex(const Cell_handle cell, const Vertex_handle vh, unsigned int i) const
  {
    int vi= cell->index(vh);
    Facet ret(cell, (vi+i+1)%4);
#ifndef NDEBUG
    Edge e= edge_around_vertex(cell, vh, i);

    for (int j=0; j<3; ++j) {
      Vertex_handle v2= vertex(ret, j);
      if (v2 != vh) {
	CGAL_assertion(v2 != cell->vertex(e.second) &&
		       v2 != cell->vertex(e.third));
      }
    }
#endif
    return ret;
  }

  void set_label(const Edge &e, Edge_label l) const
  {
    CGAL::Kinetic::internal::set_edge_label(*this, e, l);
  }

  void clear_cell_labels(Cell_handle h) const
  {
    for (unsigned int i=0; i<4; ++i) {
      h->set_facet_label(i, Facet_label());
      for (unsigned int j=0; j<i; ++j) {
	h->set_edge_label(i,j, Edge_label());
      }
    }

  }

  Vertex_handle vertex(const Facet &f, unsigned int i) const
  {
    return vertex_of_facet(f, i);
  }

  Vertex_handle vertex(const Edge &f, unsigned int i) const
  {
    //hi_there<Edge>(f);
    //hi_there(f);
    //typedef typename Edge::first_type::value_type::Vertex_handle Q;
    //Q q;
    return vertex_of_edge(f, i);
  }

  Facet opposite(const Facet &e) const
  {
    return Facet(e.first->neighbor(e.second),
		 Tri::mirror_index(e.first, e.second));// update
  }
  Edge opposite(const Edge &e) const
  {
    return Edge(e.first, e.third, e.second);
  }
  bool equal(const Edge &e0, const Edge &e1) const
  {
    Vertex_handle e0a= e0.first->vertex(e0.second);
    Vertex_handle e0b= e0.first->vertex(e0.third);
    Vertex_handle e1a= e1.first->vertex(e1.second);
    Vertex_handle e1b= e1.first->vertex(e1.third);
    bool ret=( (e0a==e1a && e0b==e1b) || (e0a==e1b && e0b== e1a));
    /*if (verbose) {
      std::cout << "Comparing ";
      rwite_edge(e0);
      std::cout << " and ";
      write_edge(e1);
      std::cout << " and getting " << ret << std::endl;
      }*/
    return ret;
  }

  template <class Stream>
  void write_facet(const Facet &f, Stream &out) const
  {
    internal::write_facet(f, out);
  }
  template <class Stream>
  void write_edge(const Edge &e, Stream &out) const
  {
    std::vector<typename Tri::Point> pts;
    pts.push_back(vertex(e,0)->point());
    pts.push_back(vertex(e,1)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << "]";
    /*if (label(e) != Edge_label::null()){
      out << " " << label(e);
      }*/
  }

  template <class Stream>
  void write_cell(const Cell_handle h, Stream &out) const
  {
    out << "[";
    for (unsigned int i=0; i< 4; ++i) {
      out << h->vertex(i)->point();
      if (i !=3) out << " ";
    }
    out << "]";
  }

  template <class Stream>
  void write_labeled_facet(const Facet &f, Stream &out) const
  {
    std::vector<typename Tri::Point> pts;
    pts.push_back(vertex(f,0)->point());
    pts.push_back(vertex(f,1)->point());
    pts.push_back(vertex(f,2)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << ", " << pts[2] << "]";
    if (label(f).is_valid()) {
      out << " " << label(f);
    }
  }

  template <class Stream>
  void write_labeled_edge(const Edge &e, Stream &out) const
  {
    std::vector<typename Tri::Point> pts;
    pts.push_back(vertex(e,0)->point());
    pts.push_back(vertex(e,1)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << "]";
    if (label(e).is_valid() ) {
      out << " " << label(e);
    }
  }

  Edge edge_in_cell(const Edge &e, const Cell_handle c) const
  {
    return CGAL::Kinetic::internal::edge_in_cell(e, c);
  }
protected:

};

} } } //namespace CGAL::Kinetic::internal
#endif
