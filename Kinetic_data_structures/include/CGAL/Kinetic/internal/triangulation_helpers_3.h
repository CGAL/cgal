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

#ifndef CGAL_KINETIC_INTERNAL_TRIANGULATION_HELPERS_3_H
#define CGAL_KINETIC_INTERNAL_TRIANGULATION_HELPERS_3_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/utility.h>
#include <vector>

namespace CGAL { namespace Kinetic { namespace internal {

/*template <class C>
typename C::first_type hi_there(C){
  return NULL;
  }*/

template <class C>
typename C::first_type::value_type::Vertex_handle vertex_of_facet(const C& f,
unsigned int i)
{
    CGAL_precondition( i<3);
    int o;
    if (f.second%2==1) o = (f.second+i+1)%4;
    else {
        o= (f.second+(2-i)+1)%4;
    }
    return f.first->vertex(o);
}


template <class C>
typename C::first_type::value_type::Vertex_handle vertex_of_edge(const C& f,
unsigned int i)
{
    CGAL_precondition(i<2);
    if (i==0) {
        return f.first->vertex(f.second);
    }
    else {
        return f.first->vertex(f.third);
    }
}


template <class F>
void write_facet(F f, std::ostream &out)
{
    std::vector<typename F::first_type::value_type::Point> pts;
    pts.push_back(vertex_of_facet(f,0)->point());
    pts.push_back(vertex_of_facet(f,1)->point());
    pts.push_back(vertex_of_facet(f,2)->point());
    std::sort(pts.begin(), pts.end());
    out << "{" << pts[0] << ", " << pts[1] << ", " << pts[2] << "}";
}


template <class E>
void write_edge(const E &e, std::ostream &out)
{
    std::vector<typename E::first_type::value_type::Point> pts;
    pts.push_back(vertex_of_edge(e,0)->point());
    pts.push_back(vertex_of_edge(e,1)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << "]";
}


template <class T>
bool has_degree_3(const T&t, const typename T::Edge &e)
{
    typename T::Cell_circulator ccir= t.incident_cells(e), ecir= ccir;
    int degree=0;
    do {
        ++degree;
        if (degree >3) break;
        ++ccir;
    } while (ccir != ecir);
    bool ret=( degree==3);
    return ret;
}

template <class Tr>
bool has_degree_4(const Tr &t, const typename Tr::Vertex_handle vh)
{
    typename Tr::Cell_handle h= vh->cell();
    int vi= h->index(vh);
    typename Tr::Vertex_handle ovh;
    bool ret=true;
    for (int i=0; i< 4; ++i) {
        if (i==vi) continue;
        if (ovh== typename Tr::Vertex_handle()) {
            ovh= t.mirror_vertex(h, i);           //h->mirror_vertex(i);
        }
        else {
            if (ovh != t.mirror_vertex(h, i)) {
                ret=false;
                break;
            }
        }
    }
#ifndef NDEBUG
    std::vector<typename Tr::Vertex_handle> vhs;
    t.incident_vertices(vh, std::back_insert_iterator<std::vector<typename Tr::Vertex_handle> >(vhs));
    if (vhs.size()==4) CGAL_postcondition(ret);
    else CGAL_postcondition(!ret);
#endif
    return ret;
}


template <class Stream, class Ch>
void write_cell(const Ch h, Stream &out)
{
    out << "[";
    for (unsigned int i=0; i< 4; ++i) {
        out << h->vertex(i)->point();
        if (i !=3) out << " ";
    }
    out << "]";
}


template <class Edge, class Cell_handle>
Edge edge_in_cell(const Edge &e, const Cell_handle c)
{
    if (e.first == c) return e;
    typename Cell_handle::value_type::Vertex_handle p0= e.first->vertex(e.second);
    typename Cell_handle::value_type::Vertex_handle p1= e.first->vertex(e.third);
    return Edge(c, c->index(p0), c->index(p1));
}


template <class Edge>
typename Edge::first_type::value_type::Edge_label edge_label(const Edge &f);

template <class Tr>
void set_edge_label(const Tr &tr, const typename Tr::Edge &e, typename Tr::Cell::Edge_label l)
{
    typename Tr::Cell_circulator cc= tr.incident_cells(e), ce= cc;
    do {
        typename Tr::Edge ec= edge_in_cell(e, cc);
        ec.first->set_edge_label(ec.second, ec.third, l);
//CGAL_postcondition(label(Edge(ec))==l);
        ++cc;
    } while (cc != ce);
    CGAL_postcondition(edge_label(e)==l);
}


template <class Facet, class Edge>
typename Facet::first_type::value_type::Vertex_handle other_vertex(const Facet &f, const Edge &e)
{
//CGAL_precondition(e.first == f.first);
    for (int i=0; i< 4; ++i) {
        if (i== f.second) continue;
        if (f.first->vertex(i) == e.first->vertex(e.second)) continue;
        if (f.first->vertex(i) == e.first->vertex(e.third)) continue;
        return f.first->vertex(i);
    }
    CGAL_postcondition(0);
    return typename Facet::first_type::value_type::Vertex_handle();
}


template <class Facet>
typename Facet::first_type::value_type::Triangulation_data_structure::Edge facet_edge(const Facet &f, unsigned int i)
{
    CGAL_precondition( i <3);
    int f0= (i+1)%3;
    int f1= (i+2)%3;
    int i0= (f.second+1+f0)%4;
    int i1= (f.second+1+f1)%4;
    if (i1== f.second) {
        i1= (i1+1)%4;
    }
    typename Facet::first_type::value_type::Triangulation_data_structure::Edge ret(f.first, i0, i1);
    CGAL_postcondition(i0 != f.second);
    CGAL_postcondition(i1 != f.second);
    return ret;
}


template <class Tr>
bool has_degree_3(const Tr &t, const typename Tr::Facet &f)
{
/*CGAL_precondition_code(Edge e0= edge(f,0));
CGAL_precondition_code(Edge e1= edge(f,1));
CGAL_precondition_code(Edge e2= edge(f,2));
CGAL_precondition(vertex(e0,1)== vertex(e1,0));
CGAL_precondition(vertex(e1,1)== vertex(e2,0));
CGAL_precondition(vertex(e2,1)== vertex(e0,0));
CGAL_precondition(vertex(e0,0) != f.first->vertex(f.second));
CGAL_precondition(vertex(e1,0) != f.first->vertex(f.second));
CGAL_precondition(vertex(e2,0) != f.first->vertex(f.second));
CGAL_precondition(hd3(f)== hd3(opposite(f)));*/
    for (unsigned int i=0; i<3 ;++i) {
        bool d3 = has_degree_3(t, facet_edge(f, i));
        CGAL_assertion(has_degree_3(t, opposite_edge(facet_edge(f,i)))== d3);
        if (d3) return true;
    }
    return false;
}


template <class Facet>
typename Facet::first_type::value_type::Facet_label facet_label(const Facet &f)
{
    return f.first->facet_label(f.second);
}


template <class Edge>
typename Edge::first_type::value_type::Edge_label edge_label(const Edge &f)
{
    return f.first->edge_label(f.second, f.third);

/*CGAL_postcondition_code(Cell_circulator cc= incident_cells(f));
CGAL_postcondition_code(Cell_circulator ce= cc);
CGAL_postcondition_code(do {)
            CGAL_postcondition_code( Edge ec= edge_in_cell(f, cc));
            CGAL_postcondition(ec.first->edge_label(ec.second, ec.third)==ret);
            CGAL_postcondition_code(++cc);
            CGAL_postcondition_code(} while (cc != ce));
            return ret;*/
}


template <class Facet, class Label>
void set_oriented_facet_label(const Facet &f, Label l)
{
    f.first->set_facet_label(f.second, l);
}


template <class Triangulation, class Facet, class Label>
void set_facet_label(const Triangulation &t, const Facet &f, Label l)
{
    set_oriented_facet_label(f, l);
    set_oriented_facet_label(opposite_facet(t, f), l);
    CGAL_postcondition(facet_label(f)==l);
}


template <class Edge>
Edge cross_edge(Edge &e)
{
    int a=-1, b;
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


template <class Cell_handle, class Vertex_handle>
typename Cell_handle::first_type::value_type::Triangulation::Edge edge_around_vertex(const Cell_handle cell,
const Vertex_handle vh,
unsigned int i)
{
    CGAL_assertion(i<3);
    int vi= cell->index(vh);
    typename Cell_handle::first_type::value_type::Triangulation::Edge ret(cell, vi, (vi+i+1)%4);
//CGAL_assertion(is_edge(ret.first, ret.second, ret.third));
    return ret;
}


//! Get a facet around the vh
/*!
  The facet i should be opposite edge i.
*/
template <class Cell_handle, class Vertex_handle>
typename Cell_handle::first_type::value_type::Triangulation::Facet facet_around_vertex(const Cell_handle cell,
const Vertex_handle vh,
unsigned int i)
{
    int vi= cell->index(vh);
    typename Cell_handle::first_type::value_type::Triangulation::Facet ret(cell, (vi+i+1)%4);
#ifndef NDEBUG
    typename Cell_handle::first_type::value_type::Triangulation::Edge e= edge_around_vertex(cell, vh, i);

    for (int j=0; j<3; ++j) {
        typename Cell_handle::first_type::value_type::Triangulation::Vertex_handle v2= vertex(ret, j);
        if (v2 != vh) {
            CGAL_assertion(v2 != cell->vertex(e.second) &&
                v2 != cell->vertex(e.third));
        }
    }
#endif
    return ret;
}


template <class Cell_handle>
void clear_cell_labels(Cell_handle h)
{
    for (unsigned int i=0; i<4; ++i) {
        h->set_facet_label(i, typename Cell_handle::value_type::Facet_label());
        for (unsigned int j=0; j<i; ++j) {
            h->set_edge_label(i,j, typename Cell_handle::value_type::Edge_label());
        }
    }
}


/*  Vertex_handle vertex(const Facet &f, unsigned int i) const {
    return vertex_of_facet(f, i);
  }

  Vertex_handle vertex(const Edge &f, unsigned int i) const {
    //hi_there<Edge>(f);
    //hi_there(f);
    typedef typename Edge::first_type::value_type::Vertex_handle Q;
    Q q;
    return vertex_of_edge(f, i);
    }*/

template <class Triangulation, class Facet>
Facet opposite_facet(const Triangulation &t, const Facet &e)
{
    return Facet(e.first->neighbor(e.second),
        t.mirror_index(e.first, e.second));
}


template <class Edge>
Edge opposite_edge(const Edge &e)
{
    return Edge(e.first, e.third, e.second);
}


template <class Edge>
bool equal_edge(const Edge &e0, const Edge &e1)
{
    typename Edge::first_type::value_type::Triangulation::Vertex_handle e0a= e0.first->vertex(e0.second);
    typename Edge::first_type::value_type::Triangulation::Vertex_handle e0b= e0.first->vertex(e0.third);
    typename Edge::first_type::value_type::Triangulation::Vertex_handle e1a= e1.first->vertex(e1.second);
    typename Edge::first_type::value_type::Triangulation::Vertex_handle e1b= e1.first->vertex(e1.third);
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


template <class Edge, class Stream>
void write_edge(const Edge &e, Stream &out)
{
    std::vector<typename Edge::first_type::value_type::Geom_traits::Point_3> pts;
    pts.push_back(vertex(e,0)->point());
    pts.push_back(vertex(e,1)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << "]";
/*if (label(e) != Edge_label::null()){
  out << " " << label(e);
  }*/
}


/*
  template <class Stream>
  void write_labeled_facet(const Facet &f, Stream &out) const {
    std::vector<typename Tri::Geom_traits::Point_3> pts;
    pts.push_back(vertex(f,0)->point());
    pts.push_back(vertex(f,1)->point());
    pts.push_back(vertex(f,2)->point());
    std::sort(pts.begin(), pts.end());
    out << "[" << pts[0] << ", " << pts[1] << ", " << pts[2] << "]";
    if (label(f)){
      out << " " << label(f);
}
}

template <class Stream>
void write_labeled_edge(const Edge &e, Stream &out) const {
std::vector<typename Tri::Geom_traits::Point_3> pts;
pts.push_back(vertex(e,0)->point());
pts.push_back(vertex(e,1)->point());
std::sort(pts.begin(), pts.end());
out << "[" << pts[0] << ", " << pts[1] << "]";
if (label(e) ){
out << " " << label(e);
}
}
*/

} } } //namespace CGAL::Kinetic::internal
#endif
