// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_COMPUTE_ANCHOR_3_H
#define CGAL_COMPUTE_ANCHOR_3_H

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_simplex_3.h>

namespace CGAL {

template < class RegularTriangulation_3>
class Compute_anchor_3
{
public:
  typedef RegularTriangulation_3                       Regular_triangulation;
  typedef typename Regular_triangulation::Geom_traits  Geom_traits;
  typedef typename Geom_traits::Weighted_point         Weighted_point;
  
  typedef typename RegularTriangulation_3::Vertex_handle Vertex_handle;
  typedef typename RegularTriangulation_3::Cell_handle   Cell_handle;
  typedef typename RegularTriangulation_3::Facet         Facet;
  typedef typename RegularTriangulation_3::Edge          Edge;
  
  typedef typename RegularTriangulation_3::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename RegularTriangulation_3::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename RegularTriangulation_3::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename RegularTriangulation_3::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename RegularTriangulation_3::Facet_circulator         Facet_circulator;
  typedef Triangulation_simplex_3<RegularTriangulation_3>           Simplex;

  typedef typename std::list<Simplex>::iterator                     Simplex_iterator;

  Compute_anchor_3(const RegularTriangulation_3 &reg) : reg(reg) {
  }

  Simplex anchor_del( const Vertex_handle v ) {
    equiv_anchors.clear();
    return v;
  }
  Simplex anchor_del( const Edge &e ) {
    return compute_anchor_del(e);
  }
  Simplex anchor_del( const Facet &f ) {
    return compute_anchor_del(f);
  }
  Simplex anchor_del( const Cell_handle ch ) {
    return compute_anchor_del(ch);
  }
  Simplex anchor_del( const Simplex &s ) {
    int dim = s.dimension();
    if (dim == 0) {
      Vertex_handle vh = s;
      return anchor_del(vh);
    } else if (dim == 1) {
      Edge e = s;
      return anchor_del(e);
    } else if (dim == 2) {
      Facet f = s;
      return anchor_del(f);
    } else if (dim == 3) {
      Cell_handle ch = s;
      return anchor_del(ch);
    }
    CGAL_error();
    return Simplex();
  }
  Simplex anchor_vor( const Vertex_handle v ) {
    return compute_anchor_vor(v);
  }
  Simplex anchor_vor( const Edge &e ) {
    return compute_anchor_vor(e);
  }
  Simplex anchor_vor( const Facet &f ) {
    return compute_anchor_vor(f);
  }
  Simplex anchor_vor( const Cell_handle ch ) {
    equiv_anchors.clear();
    for (int i=0; i<4; i++) {
      Cell_handle ch2 = ch->neighbor(i);
      if (!reg.is_infinite(ch2)) {
        Sign side = test_anchor(ch,ch2);
        CGAL_assertion(test_anchor(ch2,ch)==side);
        if (side==ZERO) {
          equiv_anchors.push_back(ch2);
        } else {
          CGAL_assertion(side==POSITIVE);
        }
      }
    }
    return ch;
  }
  Simplex anchor_vor( const Simplex &s ) {
    int dim = s.dimension();
    if (dim == 0) {
      return anchor_vor(Vertex_handle(s));
    } else if (dim == 1) {
      return anchor_vor(Edge(s));
    } else if (dim == 2) {
      return anchor_vor(Facet(s));
    } else if (dim == 3) {
      return anchor_vor(Cell_handle(s));
    }
    return Simplex();
  }

  bool is_degenerate() {
    return !equiv_anchors.empty();
  }
  Simplex_iterator equivalent_anchors_begin() {
    return equiv_anchors.begin();
  }
  Simplex_iterator equivalent_anchors_end() {
    return equiv_anchors.end();
  }
private:
  ///////////////////////////////
  // Anchor functions
  ///////////////////////////////
  Simplex compute_anchor_del( Edge const &e );
  Simplex compute_anchor_del( Facet const &f );
  Simplex compute_anchor_del( Cell_handle const ch );
  
  Simplex compute_anchor_vor( Vertex_handle const v );
  Simplex compute_anchor_vor( Edge const &e );
  Simplex compute_anchor_vor( Facet const &f );

  // Test whether the anchor of edge (wp1,wp2) and wp2 are equal
  Sign test_anchor(Weighted_point &wp1, Weighted_point &wp2) {
    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2);
  }
  Sign test_anchor(Weighted_point const& wp1, Weighted_point const& wp2, 
                   Weighted_point const& wp3) {
    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2, wp3);
  }
  Sign test_anchor(Weighted_point const& wp1, Weighted_point const& wp2, 
                   Weighted_point const& wp3, Weighted_point const& wp4) {
    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2, wp3, wp4);
  }
  // Test whether the anchor of e and anchor of e.first->vertex(i) are equal
  Sign test_anchor(Edge e, int i) {
    CGAL_assertion(!reg.is_infinite(e));
    CGAL_assertion(e.second == i || e.third == i);
    Weighted_point wp1, wp2;
    Cell_handle ch = e.first;
    return test_anchor(
      ch->vertex(e.second+e.third-i)->point(),
      ch->vertex(i)->point());
  }
  // Test whether the anchor of f and anchor of the edge f - f.first->vertex(i)
  // are equal
  Sign test_anchor(Facet f, int i){
    CGAL_assertion(!reg.is_infinite(f));
    CGAL_assertion(f.second != i);
    CGAL_assertion(0<=f.second && f.second<4);
    CGAL_assertion(0<=i && i<4);

    Weighted_point wp1, wp2, wp3;
    Cell_handle ch = f.first;
    wp3 = ch->vertex(i)->point();
    switch ((4+i-f.second)&3) {
      case 1: CGAL_assertion (((f.second+1)&3) == i);
        wp1 = ch->vertex((f.second+2)&3)->point();
        wp2 = ch->vertex((f.second+3)&3)->point();
        break;
      case 2: CGAL_assertion (((f.second+2)&3) == i);
        wp1 = ch->vertex((f.second+1)&3)->point();
        wp2 = ch->vertex((f.second+3)&3)->point();
        break;
      case 3: CGAL_assertion (((f.second+3)&3) == i);
        wp1 = ch->vertex((f.second+1)&3)->point();
        wp2 = ch->vertex((f.second+2)&3)->point();
        break;
      default:
        CGAL_error();
    }

    return 
      reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(wp1, wp2, wp3);
  }


  // Test whether the anchor of ch and anchor of the facet (ch,i) are equal
  Sign test_anchor(Cell_handle ch, int i) {
    CGAL_assertion(!reg.is_infinite(ch));

    return reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(
      ch->vertex((i+1)&3)->point(),
      ch->vertex((i+2)&3)->point(),
      ch->vertex((i+3)&3)->point(),
      ch->vertex(i)->point());                
  }
  Sign test_anchor(Cell_handle ch, Cell_handle ch2) {
    CGAL_assertion(!reg.is_infinite(ch));
    CGAL_assertion(!reg.is_infinite(ch2));

    int index = ch2->index(ch);
    return reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(
      ch->vertex(0)->point(),
      ch->vertex(1)->point(),
      ch->vertex(2)->point(),
      ch->vertex(3)->point(),
      ch2->vertex(index)->point());                
  }
  Sign test_anchor(
    Weighted_point const& wp1, Weighted_point const& wp2,
    Weighted_point const& wp3, Weighted_point const& wp4,
    Weighted_point const& wp5) {
    return reg.geom_traits().in_smallest_orthogonal_sphere_3_object()(
      wp1, wp2, wp3, wp4, wp5);                
  }

  const Regular_triangulation &reg;
  std::list<Simplex> equiv_anchors;
};


// compute_anchor_del

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Edge const &e) {
  CGAL_assertion(!reg.is_infinite(e));
  equiv_anchors.clear();
  Sign result = test_anchor(e, e.second);
  if (result==NEGATIVE) {
    return Simplex(e.first->vertex(e.third));
  } else if (result==ZERO) {
    equiv_anchors.push_back(Simplex(e.first->vertex(e.third)));
  }
  result = test_anchor(e, e.third);
  if (result==NEGATIVE) {
    equiv_anchors.clear();
    return Simplex(e.first->vertex(e.second));
  } else if (result==ZERO) {
    equiv_anchors.push_back(Simplex(e.first->vertex(e.third)));
  } 
  return Simplex(e);
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Facet const &f) {
  CGAL_assertion(!reg.is_infinite(f));
  equiv_anchors.clear();

  int i;
  Sign result;
  bool contains_center = true;
  for (i=1; (i<4) && contains_center; i++) {
    result = test_anchor(f, (f.second+i)&3);
    contains_center = (result != NEGATIVE);
    if (result == ZERO) {
      //equiv_anchors.push_back(Edge(f.first, f.second,(f.second+i)&3));
      equiv_anchors.push_back(
        Edge(f.first, (f.second+(i==1?2:1))&3, (f.second+(i==3?2:3))&3));
    }
  }
  
  if (contains_center) {
    return Simplex(f);
  } else {
    i--;
    Edge e; e.first = f.first; 
    if (i==1) e.second = ((f.second+2)&3); 
    else e.second = ((f.second+1)&3);
    if (i==3) e.third = ((f.second+2)&3); 
    else e.third = ((f.second+3)&3);

    Simplex s = anchor_del(e);
    if (s.dimension() == 1) {
      equiv_anchors.clear();
      return s;
    } else {
      // The anchor is the anchor of the other edge adjacent to tmp
      CGAL_assertion(s.dimension() == 0);
      Vertex_handle vh=s;
      e.second = e.first->index(vh);
      e.third = (f.second+i)&3;

      s = anchor_del(e);
      equiv_anchors.clear();
      return s;
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_del(Cell_handle const ch) {
  CGAL_assertion(!reg.is_infinite(ch));
  equiv_anchors.clear();
  
  Simplex s;
  bool contains_center = true;
  Sign result;
  for (int i=0; (i<4) && contains_center; i++) {
    result = test_anchor(ch, i);
    if (result == NEGATIVE) {
      contains_center = false;
      s = anchor_del(Facet(ch,i));
    } else if (result == ZERO) {
      equiv_anchors.push_back(Facet(ch,i));
    }
  }

  if (contains_center) {
    return Simplex(ch);
  } else {
    Simplex tmp;

    bool found=true;
    while (true) {
      if (s.dimension() == 1) {
        // Test two adjacent facets
        Edge e=s;
        int ind1 = ch->index(e.first->vertex(e.second));
        int ind2 = ch->index(e.first->vertex(e.third));
        for (int i=0; (i<4) && found; i++) {
          if ((i != ind1) && (i != ind2)) {
            tmp = anchor_del(Facet(ch, i));
            found = (s == tmp);
          }
        }
      } else if (s.dimension() == 0) {
        // Test adjacent edges
        Vertex_handle vh=s;
        int index = ch->index(vh);
        for (int i=0; (i<4) && found; i++) {
          if (i != index) {
            tmp = anchor_del(Edge(ch, index, i));
            found = (s == tmp);
          }
        }
      } else {
        CGAL_assertion(s.dimension() == 2);
      }
      if (found) {
        equiv_anchors.clear();
  return s;
      }
      found = true;
      s = tmp;
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::
compute_anchor_vor (Vertex_handle const v) {
  CGAL_assertion(!reg.is_infinite(v));
  CGAL_assertion(reg.is_vertex(v));

  equiv_anchors.clear();

  Sign side;
  bool contains_center=true;
  Simplex s;

  std::list<Vertex_handle> adj_vertices;
  typename std::list<Vertex_handle>::iterator adj_vertex;
  reg.incident_vertices(v, std::back_inserter(adj_vertices));
  
  for (adj_vertex = adj_vertices.begin(); 
       (adj_vertex != adj_vertices.end()) && contains_center;
       adj_vertex++) {
    if (!reg.is_infinite(*adj_vertex)) {
      side = test_anchor(v->point(),(*adj_vertex)->point());
      if (side == NEGATIVE) { 
        contains_center = false;
      } else if (side == ZERO) {
        Edge e;
        if (!reg.is_edge (v, *adj_vertex, e.first, e.second, e.third)) {
          CGAL_error();
        }
        equiv_anchors.push_back(Simplex(e));
      }
    }
  }

  if (contains_center) {
    return Simplex(v);
  } else {
    adj_vertex--;
    Edge e;
    if (!reg.is_edge(v, *adj_vertex, e.first, e.second, e.third)) {
      CGAL_error();
    }
    s = anchor_vor(e);
    Simplex tmp;
    while (true) {
      bool found = true;
      if (s.dimension() == 2) {
        // s lies on a Voronoi edge
        Facet f=s;
        int index = f.first->index(v);
        for (int i=1; (i<4) && found; i++) {
          if (((f.second+i)&3) != index) {
            tmp = anchor_vor(Edge(f.first, index, (f.second+i)&3));
            found = (tmp == s);
          }
        }
      } else if (s.dimension() == 3) {
        //  s lies on a Voronoi vertex
        Cell_handle ch=s;
        CGAL_assertion(ch != Cell_handle());
        int index = ch->index(v);
        for (int i=1; (i<4) && (found); i++) {
          tmp = anchor_vor(Facet(ch, (index+i)&3));
          found = (tmp == s);
        }
      } else {
        CGAL_assertion(s.dimension() == 1);
      }
      if (found) {
        equiv_anchors.clear();
        if (s.dimension() == 1) {
          Edge e = s;
          Vertex_handle v_other = e.first->vertex(e.second+e.third-e.first->index(v));
          for (adj_vertex = adj_vertices.begin(); 
               adj_vertex != adj_vertices.end();
               adj_vertex++) {
            if ((v_other != (*adj_vertex)) && (!reg.is_infinite(*adj_vertex))) {
              CGAL_assertion(!reg.is_infinite(v));
              CGAL_assertion(!reg.is_infinite(v_other));
              CGAL_assertion(!reg.is_infinite(*adj_vertex));
              side = test_anchor(v->point(), v_other->point(),
                                 (*adj_vertex)->point());
              if (side==ZERO) {
                Edge e2;
                if (!reg.is_edge(v, *adj_vertex, e2.first, e2.second, e2.third)) {
                  CGAL_error();
                }
                equiv_anchors.push_back(e2);
              }
            }
          }          
        }
        return s;
      }
      s = tmp;
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::compute_anchor_vor (Edge const &e) {
  CGAL_assertion(!reg.is_infinite(e));

  equiv_anchors.clear();

  Vertex_handle v0 = e.first->vertex(e.second);
  Vertex_handle v1 = e.first->vertex(e.third);

  bool contains_center = true;
  Sign side;
  Simplex s;
  
  Facet_circulator fcir, fstart;
  fstart = fcir = reg.incident_facets(e);

  do {
    if (!reg.is_infinite(*fcir)) {
      int i = 6 - (*fcir).second -
        (*fcir).first->index(v0) - (*fcir).first->index(v1);
      side = test_anchor(*fcir, i);
      if (side == NEGATIVE) {
        contains_center = false;
        s = anchor_vor(Facet(*fcir));
      } else if (side == ZERO) {
        equiv_anchors.push_back(Facet(*fcir));
      }
    }
  } while (++fcir != fstart);
  
  if (contains_center) {
    s = Simplex(e);
    return s;
  } else {
    Simplex tmp;
    while (true) {
      bool found = true;
      if (s.dimension() == 3) {
        CGAL_assertion(s.dimension() == 3);
        Cell_handle ch=s;
        int index0 = ch->index(v0);
        int index1 = ch->index(v1);
        for (int i=0; (i<4) && found; i++) {
          if ((i != index0) && (i != index1)) {
            if (!reg.is_infinite(Facet(ch,i))) {
              side = test_anchor(ch,6-index0-index1-i);
              if (side != POSITIVE) {
                tmp = anchor_vor(Facet(ch,i));
                found = (s==tmp);
              }
            }
          }
        }
      } else {
        CGAL_assertion(s.dimension() == 2);
      }
      if (found) {
        equiv_anchors.clear();
        if (s.dimension() == 2) {
          // Check whether facet is degenerate (a line segment):
          Facet f = s;
          int index = 6 - f.second - f.first->index(v0) - f.first->index(v1);
          fstart = fcir = reg.incident_facets(e);
          do {
            if (!reg.is_infinite(*fcir)) {
              int index2 = 6 - (*fcir).second
                             - (*fcir).first->index(v0) 
                             - (*fcir).first->index(v1);
            if (!(f.first->vertex(index) == (*fcir).first->vertex(index2))) {
                side = test_anchor(v0->point(), v1->point(),
                                   f.first->vertex(index)->point(), 
                                   (*fcir).first->vertex(index2)->point());
                if (side == ZERO) {
                  equiv_anchors.push_back(Facet(*fcir));
                }
              }
            }
          } while (++fcir != fstart);
        }
        return s;
      }
      s = tmp;
    }
  }
}

template < class RegularTriangulation3 >
typename Compute_anchor_3<RegularTriangulation3>::Simplex
Compute_anchor_3<RegularTriangulation3>::compute_anchor_vor (Facet const &f) {
  CGAL_assertion(!reg.is_infinite(f));
  equiv_anchors.clear();
  
  Sign side;
  
  CGAL_assertion(f.first != Cell_handle());
  if (!reg.is_infinite(f.first)) {
    side = test_anchor(f.first, f.second);
    if (side==NEGATIVE) {
      return Simplex(f.first);
    } else if (side == ZERO) {
      equiv_anchors.push_back(f.first);
    }
  }

  Cell_handle neighbor = f.first->neighbor(f.second);
  CGAL_assertion(neighbor != Cell_handle());
  if (!reg.is_infinite(neighbor)) {
    int n_index = neighbor->index(f.first);
    side = test_anchor(neighbor, n_index);
    if (side==NEGATIVE) {
      CGAL_assertion(equiv_anchors.empty());
      return Simplex(neighbor);
    } else if (side == ZERO) {
      equiv_anchors.push_back(neighbor);
    }
  }

  return Simplex(f);
}


} //namespace CGAL

#endif // CGAL_COMPUTE_ANCHOR_3_H
