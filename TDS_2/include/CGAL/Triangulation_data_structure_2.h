// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_TRIANGULATION_DATA_STRUCTURE_2_H
#define CGAL_TRIANGULATION_DATA_STRUCTURE_2_H

#include <CGAL/license/TDS_2.h>


#include <CGAL/config.h>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <vector>
#include <algorithm>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
 
#include <CGAL/Compact_container.h>

#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Triangulation_ds_iterators_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>

namespace CGAL { 

template < class Vb = Triangulation_ds_vertex_base_2<>, 
           class Fb = Triangulation_ds_face_base_2<> >
class Triangulation_data_structure_2 
  :public Triangulation_cw_ccw_2
{
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;

  typedef typename Vb::template Rebind_TDS<Tds>::Other  Vertex_base;
  typedef typename Fb::template Rebind_TDS<Tds>::Other  Face_base;

  friend class Triangulation_ds_edge_iterator_2<Tds,false>;
  friend class Triangulation_ds_edge_iterator_2<Tds,true>;
  friend class Triangulation_ds_face_circulator_2<Tds>;
  friend class Triangulation_ds_edge_circulator_2<Tds>;
  friend class Triangulation_ds_vertex_circulator_2<Tds>;

public:
  // Tools to change the Vertex and Face types of the TDS.
  template < typename Vb2 >
  struct Rebind_vertex {
    typedef Triangulation_data_structure_2<Vb2, Fb>  Other;
  };

  template < typename Fb2 >
  struct Rebind_face {
    typedef Triangulation_data_structure_2<Vb, Fb2>  Other;
  };

  typedef Vertex_base                                Vertex;
  typedef Face_base                                  Face;
  
  typedef Compact_container<Face>                    Face_range;
  typedef Compact_container<Vertex>                  Vertex_range;

  typedef typename Face_range::size_type             size_type;
  typedef typename Face_range::difference_type       difference_type;

  typedef typename Face_range::iterator              Face_iterator;
  typedef typename Vertex_range::iterator            Vertex_iterator;

  typedef Triangulation_ds_edge_iterator_2<Tds>      Edge_iterator;
  typedef Triangulation_ds_edge_iterator_2<Tds,false> Halfedge_iterator;

  typedef Triangulation_ds_face_circulator_2<Tds>    Face_circulator;
  typedef Triangulation_ds_vertex_circulator_2<Tds>  Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_2<Tds>    Edge_circulator;

  typedef Vertex_iterator                            Vertex_handle;
  typedef Face_iterator                              Face_handle;

  typedef std::pair<Face_handle,int> Edge;

  typedef std::list<Edge> List_edges;

protected:
  int _dimension;
  Face_range   _faces;
  Vertex_range _vertices;

  //CREATORS - DESTRUCTORS
public:
  Triangulation_data_structure_2(); 
  Triangulation_data_structure_2(const Tds &tds);
  ~Triangulation_data_structure_2();
  Tds& operator= (const Tds &tds);
  void swap(Tds &tds);

  //ACCESS FUNCTIONS
  // We need the const_cast<>s because TDS is not const-correct.
  Face_range& faces()             { return _faces;}
  Face_range& faces() const 
    { return  const_cast<Tds*>(this)->_faces;}
  Vertex_range& vertices()         {return _vertices;}
  Vertex_range& vertices() const
    {return  const_cast<Tds*>(this)->_vertices;}

  int  dimension() const { return _dimension;  }
  size_type number_of_vertices() const {return vertices().size();}
  size_type number_of_faces() const ;
  size_type number_of_edges() const;
  size_type number_of_full_dim_faces() const; //number of faces stored by tds
  
  // TEST FEATURES
  bool is_vertex(Vertex_handle v) const;
  bool is_edge(Face_handle fh, int i) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb, 
	       Face_handle& fr,  int& i) const;
  bool is_face(Face_handle fh) const;
  bool is_face(Vertex_handle v1, 
	       Vertex_handle v2, 
	       Vertex_handle v3) const;
  bool is_face(Vertex_handle v1, 
	       Vertex_handle v2, 
	       Vertex_handle v3,
	       Face_handle& fr) const;

  // ITERATORS AND CIRCULATORS
public:
// The face_iterator_base_begin  gives the possibility to iterate over all
// faces in the container  independently of the dimension.
  // public for the need of file_ouput() of Constrained triangulation
  // should be made private later

  Face_iterator face_iterator_base_begin() const    {
    return faces().begin();
  }
  Face_iterator face_iterator_base_end() const    {
    return faces().end();
  }

public:
  Face_iterator faces_begin() const {
    if (dimension() < 2) return faces_end();
    return faces().begin();
  }
    
  Face_iterator faces_end() const {
    return faces().end();
  }

  Vertex_iterator vertices_begin() const  {
    return vertices().begin();
  }

  Vertex_iterator vertices_end() const {
    return vertices().end();
  }
  
  Edge_iterator edges_begin() const {
    return Edge_iterator(this);
  }

  Edge_iterator edges_end() const {
    return Edge_iterator(this,1);
  }
  
  Halfedge_iterator halfedges_begin() const {
    return Halfedge_iterator(this);
  }

  Halfedge_iterator halfedges_end() const {
    return Halfedge_iterator(this,1);
  }
  
  Face_circulator incident_faces(Vertex_handle v, 
				 Face_handle f =  Face_handle()) const{
    return Face_circulator(v,f);
  }
  Vertex_circulator incident_vertices(Vertex_handle v, 
				      Face_handle f = Face_handle()) const
  {    
    return Vertex_circulator(v,f);  
  }

  Edge_circulator incident_edges(Vertex_handle v, 
				 Face_handle f = Face_handle()) const{
    return Edge_circulator(v,f);
  }

  size_type degree(Vertex_handle v) const {
    int count = 0;
    Vertex_circulator vc = incident_vertices(v), done(vc);
    if ( ! vc.is_empty()) {
      do { 
	count += 1;
      } while (++vc != done);
    }
    return count;
  }

  
  Vertex_handle
  mirror_vertex(Face_handle f, int i) const
  {
    CGAL_triangulation_precondition ( f->neighbor(i) != Face_handle()
				    && f->dimension() >= 1);
  return f->neighbor(i)->vertex(mirror_index(f,i));
  }

  int
  mirror_index(Face_handle f, int i) const
  {
    // return the index of opposite vertex in neighbor(i);
    CGAL_triangulation_precondition (f->neighbor(i) != Face_handle() &&
				     f->dimension() >= 1);
    if (f->dimension() == 1) {
      CGAL_assertion(i<=1);
      const int j = f->neighbor(i)->index(f->vertex((i==0) ? 1 : 0));
      CGAL_assertion(j<=1);
      return (j==0) ? 1 : 0;
    }
    return ccw( f->neighbor(i)->index(f->vertex(ccw(i))));
  }

  Edge 
  mirror_edge(const Edge e) const 
  {
    CGAL_triangulation_precondition(e.first->neighbor(e.second) != Face_handle()
                                    && e.first->dimension() >= 1);
    return Edge(e.first->neighbor(e.second),
                mirror_index(e.first,  e.second));
  }

  // MODIFY
  void flip(Face_handle f, int i);
 
  Vertex_handle insert_first();
  Vertex_handle insert_second();
  Vertex_handle insert_in_face(Face_handle f);
  Vertex_handle insert_in_edge(Face_handle f, int i);
  Vertex_handle insert_dim_up(Vertex_handle w = Vertex_handle(), 
			      bool orient=true);

  void remove_degree_3(Vertex_handle v, Face_handle f = Face_handle());
  void remove_1D(Vertex_handle v); 
   
  void remove_second(Vertex_handle v);
  void remove_first(Vertex_handle v);
  void remove_dim_down(Vertex_handle v);
  void dim_down(Face_handle f, int i);

  Vertex_handle star_hole(List_edges& hole);
  void    star_hole(Vertex_handle v, List_edges& hole);
  void    make_hole(Vertex_handle v, List_edges& hole);

//   template< class EdgeIt>
//   Vertex_handle star_hole(EdgeIt edge_begin,EdgeIt edge_end);
 
//   template< class EdgeIt>
//   void  star_hole(Vertex_handle v, EdgeIt edge_begin,  EdgeIt edge_end);

//   template< class EdgeIt, class FaceIt>
//   Vertex_handle star_hole(EdgeIt edge_begin, 
// 		    EdgeIt edge_end,
// 		    FaceIt face_begin,
// 		    FaceIt face_end);
 
//   template< class EdgeIt, class FaceIt>
//   void  star_hole(Vertex_handle v,
// 		  EdgeIt edge_begin, 
// 		  EdgeIt edge_end,
// 		  FaceIt face_begin,
// 		  FaceIt face_end);
  
  Vertex_handle create_vertex(const Vertex &v = Vertex());
  Vertex_handle create_vertex(Vertex_handle v); //calls copy constructor 
  Face_handle create_face(const Face& f = Face());
  Face_handle create_face(Face_handle f); //calls copy constructor 

 Face_handle create_face(Face_handle f1, int i1, 
			  Face_handle f2, int i2, 
			  Face_handle f3, int i3);
  Face_handle create_face(Face_handle f1, int i1, 
			  Face_handle f2, int i2);
  Face_handle create_face(Face_handle f1, int i1, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, 
			  Vertex_handle v2, 
			  Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, 
			  Vertex_handle v2, 
			  Vertex_handle v3,
			  Face_handle f1, 
			  Face_handle f2, 
			  Face_handle f3);

  void set_adjacency(Face_handle f0, int i0, Face_handle f1, int i1) const;

  void delete_face(Face_handle);
  void delete_vertex(Vertex_handle);

  // split and join operations
 protected:
  Vertex_handle join_vertices(Face_handle f, int i, Vertex_handle v);

  typedef
  boost::tuples::tuple<Vertex_handle,Vertex_handle,Face_handle,Face_handle>
  Fourtuple;

 public:
  Fourtuple split_vertex(Vertex_handle v, Face_handle f1, Face_handle g1);

  inline Vertex_handle join_vertices(Face_handle f, int i) {
    return join_vertices(f, i, f->vertex( ccw(i) ));
  }

  inline Vertex_handle join_vertices(Edge e) {
    return join_vertices(e.first, e.second);
  }

  inline Vertex_handle join_vertices(Edge_iterator eit) {
    return join_vertices(*eit);
  }

  inline Vertex_handle join_vertices(Edge_circulator ec) {
    return join_vertices(*ec);
  }

  // insert_degree_2 and remove_degree_2 operations
  Vertex_handle insert_degree_2(Face_handle f, int i);
  void remove_degree_2(Vertex_handle v);

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;
  
  // HELPING
private:
  typedef std::pair<Vertex_handle,Vertex_handle> Vh_pair;
public:
  void  set_adjacency(Face_handle fh, 
		      int ih, 
		      std::map< Vh_pair, Edge>& edge_map);
  void reorient_faces();
private:
  bool dim_down_precondition(Face_handle f, int i);

public:
  void clear();

  template <class TDS_src>
  Vertex_handle copy_tds(const TDS_src &tds, typename TDS_src::Vertex_handle);

  template <class TDS_src>
  Vertex_handle copy_tds(const TDS_src &tds)
  {
    return copy_tds(tds, typename TDS_src::Vertex_handle());
  }

  template <class TDS_src,class ConvertVertex,class ConvertFace>
  Vertex_handle copy_tds(const TDS_src&, typename TDS_src::Vertex_handle,const ConvertVertex&,const ConvertFace&);

  Vertex_handle collapse_edge(Edge e)
  {
    std::cout << "before collapse"<<std::endl;
    Face_handle fh = e.first;
    int i = e.second;
    Vertex_handle vh = fh->vertex(cw(i));
    Vertex_handle wh = fh->vertex(ccw(i));
    Face_handle left = fh->neighbor(cw(i));
    Face_handle right = fh->neighbor(ccw(i));
    Face_handle nh = fh->neighbor(i);
    int li = left->index(fh);
    int ri = right->index(fh);
    int ni = nh->index(fh);
    left->set_neighbor(li, right);
    right->set_neighbor(ri,left);
    left->set_vertex(ccw(li), vh);
    vh->set_face(right);
    right->vertex(ccw(ri))->set_face(right);

    left = nh->neighbor(ccw(ni));
    right = nh->neighbor(cw(ni));
    li = left->index(nh);
    ri = right->index(nh);
    left->set_neighbor(li, right);
    right->set_neighbor(ri,left);
    left->set_vertex(cw(li), vh);
    right->vertex(cw(ri))->set_face(right);
    delete_face(fh);
    delete_face(nh);
    delete_vertex(wh);
    std::cout << "after collapse"<<std::endl;
    return vh;
  }


  // I/O
  Vertex_handle file_input(std::istream& is, bool skip_first=false);
  void file_output(std::ostream& os,
		   Vertex_handle v = Vertex_handle(),
		   bool skip_first=false) const;
  Vertex_handle off_file_input(std::istream& is, bool verbose=false);
  void  vrml_output(std::ostream& os,
		    Vertex_handle v = Vertex_handle(),
		    bool skip_first=false) const;

  // SETTING (had to make them public for use in remove from Triangulations)
  void set_dimension (int n) {_dimension = n ;}

  // template members definition
public:

  /************* START OF MODIFICATIONS ***************/

  template< class FaceIt >
  Vertex_handle insert_in_hole(FaceIt face_begin, FaceIt face_end) 
  {
    Vertex_handle newv = create_vertex();
    insert_in_hole(newv, face_begin, face_end);
    return newv;
  }


  template< class FaceIt >
  void insert_in_hole(Vertex_handle v, FaceIt face_begin, FaceIt face_end) 
  {

    CGAL_triangulation_precondition(dimension() == 2);

    std::vector<Face_handle>  new_faces;
    std::vector<Edge>         bdry_edges;

    Face_handle fh = *face_begin;
    int ii = 0;
    bool found_boundary = false;
    do {
      if (std::find(face_begin, face_end, fh->neighbor(ii)) == face_end) {
        bdry_edges.push_back(Edge(fh, ii));
        found_boundary = true;
      } else {
        int newi = fh->neighbor(ii)->index(fh->vertex(ccw(ii)));
        fh = fh->neighbor(ii);
        ii = newi;
      }
    } while(!found_boundary);
    // Now we have found ONE edge on the boundary. 
    // From that one edge we must walk on the boundary 
    // of the hole until we've covered the whole thing.

    bool complete_walk = false;
    do {
      Face_handle nh = fh->neighbor(ccw(ii));
      if (std::find(face_begin, face_end, nh) == face_end) {
        ii = ccw(ii);
        Edge new_edge(fh, ii);
        if (std::find(bdry_edges.begin(), bdry_edges.end(), new_edge) == bdry_edges.end()) {
          bdry_edges.push_back(Edge(fh, ii));
        } else {
          complete_walk = true;
        }
      } else {
        int newi = cw(nh->index(fh->vertex(ii)));
        fh = nh;
        ii = newi;
      }
    } while (!complete_walk);
    // At this point, bdry_edges contains the edges that define
    // the boundary of the hole with a specific ordering: for any
    // two consecutive edges in the vector e1 = (f1, i1), 
    // e2 = (f2, i2) it holds that 
    //      f1->vertex(cw(i1)) == f2->vertex(ccw(i2))

    for (unsigned int jj = 0; jj < bdry_edges.size(); jj++) {
      Face_handle fh = bdry_edges[jj].first;
      int idx = bdry_edges[jj].second;
      
      Vertex_handle v1 = fh->vertex(ccw(idx));
      Vertex_handle v2 = fh->vertex(cw(idx));

      Face_handle nf = fh->neighbor(idx);
      int jdx = mirror_index(fh, idx);

      Face_handle new_f = create_face(v, v1, v2);
      v1->set_face(new_f);
      set_adjacency(new_f, 0, nf, jdx);
      new_faces.push_back(new_f);
    }
    // At this point we have created all the new faces of the triangulation,
    // and we have set adjacency relationships with the faces on the border
    // of the hole.

    for (unsigned int i = 0; i < new_faces.size() - 1; i++) {
      set_adjacency(new_faces[i], 1, new_faces[i+1], 2);
    }
    set_adjacency(new_faces[0], 2, new_faces[new_faces.size()-1], 1);
    // Now we have also set adjacency relationships between the new faces.

    for (FaceIt it = face_begin; it != face_end; it++) {
      delete_face(*it);
    }
    // The old faces that were in conflict are now deleted.

    v->set_face(new_faces[0]);
    // Set the pointer of the new vertex to one of the new faces.
  }



  /************* END OF MODIFICATIONS ***************/


  template< class EdgeIt>
  Vertex_handle star_hole(EdgeIt edge_begin, EdgeIt edge_end)
  // creates a new vertex 
  // and stars from it
  // the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
  {
     Vertex_handle newv = create_vertex();
     star_hole(newv, edge_begin, edge_end);
     return newv;
  }
 
  template< class EdgeIt>
  void  star_hole(Vertex_handle v, EdgeIt edge_begin,  EdgeIt edge_end)
  // uses vertex v
  // to star the hole described by the range [edge_begin,edge_end[
  // the triangulation is assumed to have dim=2
  // the hole is supposed to be ccw oriented
  { 
    std::list<Face_handle> empty_list;
    star_hole(v, 
	      edge_begin, 
	      edge_end, 
	      empty_list.begin(),
	      empty_list.end());
    return;    
  }


  template< class EdgeIt, class FaceIt>
  Vertex_handle star_hole(EdgeIt edge_begin, 
		    EdgeIt edge_end,
		    FaceIt face_begin,
		    FaceIt face_end)
  // creates a new vertex 
  // and stars from it
  // the hole described by the range [edge_begin,edge_end[
    // reusing the faces in the range [face_begin,face_end[
  // the triangulation is assumed to have dim=2
  // the hole is supposed to be ccw oriented
  {
    Vertex_handle newv = create_vertex();
    star_hole(newv, edge_begin, edge_end, face_begin, face_end);
    return newv;
  }
 
  template< class EdgeIt, class FaceIt>
  void  star_hole(Vertex_handle newv,
		  EdgeIt edge_begin, 
		  EdgeIt edge_end,
		  FaceIt face_begin,
		  FaceIt face_end)
    // uses vertex v
    // to star the hole described by the range [edge_begin,edge_end[
    // reusing the faces in the range [face_begin,face_end[
    // the triangulation is assumed to have dim=2
    // hole is supposed to be ccw oriented
  {
    CGAL_triangulation_precondition(dimension() == 2);
    EdgeIt eit = edge_begin;
    FaceIt fit = face_begin;

    Face_handle fn = (*eit).first;
    int in = (*eit).second;
    fn->vertex(cw(in))->set_face(fn);
    Face_handle first_f =  reset_or_create_face(fn, in , newv, fit, face_end);
    Face_handle previous_f=first_f, next_f;
    ++eit; 

    for( ; eit != edge_end ; eit++) {
      fn = (*eit).first;
      in = (*eit).second;
      fn->vertex(cw(in))->set_face(fn);
      next_f = reset_or_create_face(fn, in , newv, fit, face_end);
      set_adjacency(next_f, 1, previous_f, 0);
      previous_f=next_f;
    }

    set_adjacency(next_f, 0, first_f, 1);
    newv->set_face(first_f);
    return;    
  }

private:
  template< class FaceIt>
  Face_handle  reset_or_create_face(Face_handle fn, 
			      int in, 
			      Vertex_handle v,
			      FaceIt& fit,
			      const FaceIt& face_end)
  {
    if (fit == face_end) return create_face(fn, in, v);
    (*fit)->set_vertices(fn->vertex(cw(in)), fn->vertex(ccw(in)), v);
    (*fit)->set_neighbors(Face_handle(),Face_handle(),fn);
    fn->set_neighbor(in, *fit);
    return *fit++;    
  }

};

//for backward compatibility
template < class Gt , class Vb, class Fb>
class Triangulation_default_data_structure_2 
  : public Triangulation_data_structure_2<Vb,Fb>
{
public:
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;
  typedef Triangulation_default_data_structure_2<Gt,Vb,Fb> Tdds;
  typedef Gt                           Geom_traits; 

  Triangulation_default_data_structure_2(const Geom_traits& = Geom_traits())
    : Tds() {}
 
  Triangulation_default_data_structure_2(const Tdds &tdds)
    : Tds(tdds) {}
};

//for backward compatibility
template <class Vb, class Fb>
class Triangulation_data_structure_using_list_2
  :public Triangulation_data_structure_2<Vb, Fb>
{
public:
  typedef Triangulation_data_structure_2<Vb,Fb>  Tds;
  typedef Triangulation_data_structure_using_list_2<Vb,Fb>  Tdsul;

  Triangulation_data_structure_using_list_2(): Tds() {} 
  Triangulation_data_structure_using_list_2(const Tdsul &tdsul)
    : Tds(tdsul) {}
};

 
template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
Triangulation_data_structure_2() 
  : _dimension(-2)
{ }

template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
Triangulation_data_structure_2(const Tds &tds)
{
  copy_tds(tds);
}

template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb> ::
~Triangulation_data_structure_2()
{
  clear();
}

//assignement  
template < class Vb, class Fb>
Triangulation_data_structure_2<Vb,Fb>&
Triangulation_data_structure_2<Vb,Fb> ::
operator= (const Tds &tds)
{
  copy_tds(tds);
  return *this;
}  

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
clear()
{
  faces().clear();
  vertices().clear();
  set_dimension(-2);
  return;
}

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
swap(Tds &tds)
{
  CGAL_triangulation_expensive_precondition(tds.is_valid() && is_valid());
  std::swap(_dimension, tds._dimension);
  faces().swap(tds.faces());
  vertices().swap(tds.vertices());
  return;
}

//ACCESS FUNCTIONS
template <class Vb, class Fb>
inline 
typename Triangulation_data_structure_2<Vb,Fb>::size_type
Triangulation_data_structure_2<Vb,Fb> ::
number_of_faces() const 
{
  if (dimension() < 2) return 0;
  return faces().size();
}

template <class Vb, class Fb>
inline 
typename Triangulation_data_structure_2<Vb,Fb>::size_type
Triangulation_data_structure_2<Vb,Fb>::
number_of_edges() const
{
  switch (dimension()) {
  case 1:  return number_of_vertices();
  case 2:  return 3*number_of_faces()/2;
  default: return 0;
  }
}
      
template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::size_type
Triangulation_data_structure_2<Vb,Fb>::
number_of_full_dim_faces() const
{
  return faces().size();
}

template <class Vb, class Fb>
inline bool
Triangulation_data_structure_2<Vb,Fb>::
is_vertex(Vertex_handle v) const
{
  Vertex_iterator vit = vertices_begin();
  while (vit != vertices_end() && v != vit)
        ++vit;
  return v == vit;
}

template <class Vb, class Fb>
inline bool
Triangulation_data_structure_2<Vb,Fb>::
is_edge(Face_handle fh, int i) const
{
  if ( dimension() == 0 )  return false;
  if ( dimension() == 1 && i != 2) return false;
  if (i > 2) return false;
  Face_iterator fit = face_iterator_base_begin();
  while (fit != face_iterator_base_end() && fh != fit ) ++fit;
  return fh == fit;
}

template <class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_edge(Vertex_handle va, Vertex_handle vb) const
// returns true (false) if the line segment ab is (is not) an edge of t
//It is assumed that va is a vertex of t
{
  Vertex_circulator vc = incident_vertices(va), done(vc);
  if ( vc == 0) return false;
  do {
    if( vb == vc ) {return true;} 
  } while (++vc != done);
  return false;
}
 

template <class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_edge(Vertex_handle va, Vertex_handle vb, 
	Face_handle &fr,  int & i) const
// assume va is a vertex of t
// returns true (false) if the line segment ab is (is not) an edge of t
// if true is returned (fr,i) is the edge ab
// with face fr on the right of a->b
{
  Face_handle fc = va->face(); 
  Face_handle start = fc;
  if (fc == 0) return false;
  int inda, indb;
  do {
    inda=fc->index(va);
    indb = (dimension() == 2 ? cw(inda) : 1-inda);
    if(fc->vertex(indb) == vb) {
      fr=fc;
      i = 3 - inda - indb; //works in dim 1 or 2
      return true;
    }
    fc=fc->neighbor(indb); //turns ccw around va
  } while (fc != start);
  return false;
}

template <class Vb, class Fb>
inline bool 
Triangulation_data_structure_2<Vb,Fb>::
is_face(Face_handle fh) const
{
  if (dimension() < 2)  return false;
  Face_iterator fit = faces_begin();
  while (fit != faces_end() && fh != fit ) ++fit;
  return fh == fit;
}

template <class Vb, class Fb>
inline bool 
Triangulation_data_structure_2<Vb,Fb>::
is_face(Vertex_handle v1, 
	Vertex_handle v2, 
	Vertex_handle v3) const
{
  Face_handle f;
  return is_face(v1,v2,v3,f);
}

template <class Vb, class Fb>
bool 
Triangulation_data_structure_2<Vb,Fb>::
is_face(Vertex_handle v1, 
	Vertex_handle v2, 
	Vertex_handle v3,
	Face_handle &f) const
{
  if (dimension() != 2) return false;
  int i;
  bool b = is_edge(v1,v2,f,i);
  if (!b) return false;
  else if (v3== f->vertex(i)) return true;
  f = f-> neighbor(i);
  int ind1= f->index(v1);
  int ind2= f->index(v2);
  if (v3 == f->vertex(3-ind1-ind2)) { return true;}
  return false;  
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
flip(Face_handle f, int i)
{
  CGAL_triangulation_precondition( dimension()==2);
  Face_handle n  = f->neighbor(i);
  int ni = mirror_index(f,i); //ni = n->index(f);
    
  Vertex_handle  v_cw = f->vertex(cw(i));
  Vertex_handle  v_ccw = f->vertex(ccw(i));

  // bl == bottom left, tr == top right
  Face_handle tr = f->neighbor(ccw(i));
  int tri =  mirror_index(f,ccw(i));  
  Face_handle bl = n->neighbor(ccw(ni));
  int bli =  mirror_index(n,ccw(ni)); 
      
  f->set_vertex(cw(i), n->vertex(ni));
  n->set_vertex(cw(ni), f->vertex(i));
    
  // update the neighborhood relations
  set_adjacency(f, i, bl, bli);
  set_adjacency(f, ccw(i), n, ccw(ni));
  set_adjacency(n, ni, tr, tri);

  if(v_cw->face() == f) {
    v_cw->set_face(n);
  }
    
  if(v_ccw->face() == n) {
    v_ccw->set_face(f);
  }
}
  
template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_first( )
{
  CGAL_triangulation_precondition( number_of_vertices() == 0 &&
				   dimension()==-2 );
  return insert_dim_up();
}

template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle 
Triangulation_data_structure_2<Vb,Fb>::
insert_second()
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 &&
				   dimension()==-1 );
  return insert_dim_up();

}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_in_face(Face_handle f)
  // New vertex will replace f->vertex(0) in face f
{
  CGAL_triangulation_precondition( f != Face_handle() && dimension()== 2);
  Vertex_handle  v = create_vertex();

  Vertex_handle v0 = f->vertex(0);
  Vertex_handle v2 = f->vertex(2);
  Vertex_handle v1 = f->vertex(1);
    
  Face_handle n1 = f->neighbor(1);
  Face_handle n2 = f->neighbor(2);
    
  Face_handle f1 = create_face(v0, v, v2, f, n1, Face_handle());
  Face_handle f2 = create_face(v0, v1, v, f, Face_handle(), n2);

  set_adjacency(f1, 2, f2, 1);
  if (n1 != Face_handle()) {
    int i1 = mirror_index(f,1); //int i1 = n1->index(f);
    n1->set_neighbor(i1,f1);
  }
  if (n2 != Face_handle()) {
    int i2 = mirror_index(f,2);//int i2 = n2->index(f);
    n2->set_neighbor(i2,f2);}

  f->set_vertex(0, v);
  f->set_neighbor(1, f1);
  f->set_neighbor(2, f2);

  if( v0->face() == f  ) {  v0->set_face(f2); }
  v->set_face(f);

  return v;
}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_in_edge(Face_handle f, int i)
  //insert in the edge opposite to vertex i of face f
{
  CGAL_triangulation_precondition(f != Face_handle() && dimension() >= 1); 
  if (dimension() == 1) {CGAL_triangulation_precondition(i == 2);}
  if (dimension() == 2) {CGAL_triangulation_precondition(i == 0 || 
							 i == 1 || 
							 i == 2);}
  Vertex_handle v;
  if (dimension() == 1) {
    v = create_vertex();
    Face_handle ff = f->neighbor(0);
    Vertex_handle vv = f->vertex(1);
    Face_handle g = create_face(v,vv,Vertex_handle(),ff, f, Face_handle());
    f->set_vertex(1,v);f->set_neighbor(0,g);
    ff->set_neighbor(1,g);
    v->set_face(g);
    vv->set_face(ff);
  }

    else { //dimension() ==2
    Face_handle n = f->neighbor(i);
    int in = mirror_index(f,i); //n->index(f);
    v = insert_in_face(f);
    flip(n,in); 
    }

  return v;
}


template <  class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_dim_up(Vertex_handle w,  bool orient)
{
  // the following function insert 
  // a vertex  v which is outside the affine  hull of Tds
  // The triangulation will be starred from  v and w 
  // ( geometrically w=  // the infinite vertex )
  // w=NULL for first and second insertions
  // orient governs the orientation of the resulting triangulation

  Vertex_handle v = create_vertex();
  set_dimension( dimension() + 1);
  Face_handle f1;
  Face_handle f2;

  const int dim = dimension(); //it is the resulting dimension
    
  switch (dim) { 
  case -1:
    f1 = create_face(v,Vertex_handle(),Vertex_handle());
    v->set_face(f1);
    break;
  case 0 :
    f1 = face_iterator_base_begin();
    f2 = create_face(v,Vertex_handle(),Vertex_handle());
    set_adjacency(f1, 0, f2, 0);
    v->set_face(f2);
    break;
  case 1 :
  case 2 :
    {
      std::list<Face_handle> faces_list;
      Face_iterator ib= face_iterator_base_begin(); 
      Face_iterator ib_end = face_iterator_base_end();
      for (; ib != ib_end ; ++ib){
	faces_list.push_back( ib);
      }
      
      std::list<Face_handle>  to_delete;
      typename std::list<Face_handle>::iterator lfit = faces_list.begin();
      Face_handle f, g;

      for ( ; lfit != faces_list.end() ; ++lfit) {
	f = * lfit;
	g = create_face(f); //calls copy constructor of face
	f->set_vertex(dim,v);
	g->set_vertex(dim,w);
	set_adjacency(f, dim, g, dim);
	if (f->has_vertex(w)) to_delete.push_back(g); // flat face to delete
      }

      lfit = faces_list.begin();
      for ( ; lfit != faces_list.end() ; ++lfit) {
	f = * lfit;
	g = f->neighbor(dim);
	for(int j = 0; j < dim ; ++j) {
	  g->set_neighbor(j, f->neighbor(j)->neighbor(dim));
	}
      }

      // couldn't unify the code for reorientation mater
      lfit = faces_list.begin() ; 
      if (dim == 1){
	if (orient) {
	  (*lfit)->reorient(); ++lfit ;  (*lfit)->neighbor(1)->reorient();
	}
	else {
	  (*lfit)->neighbor(1)->reorient(); ++lfit ; (*lfit)->reorient(); 
	}
      }
      else { // dimension == 2
	for( ;lfit  != faces_list.end(); ++lfit ) {
	  if (orient) {(*lfit)->neighbor(2)->reorient();}
	  else { (*lfit)->reorient();}
	}
      }

      lfit = to_delete.begin();
      int i1, i2;
      for ( ;lfit  != to_delete.end(); ++lfit){
	f = *lfit ;
	int j ;
	if (f->vertex(0) == w) {j=0;}
	else {j=1;}
	f1= f->neighbor(dim); i1= mirror_index(f,dim); //f1->index(f);
	f2= f->neighbor(j); i2= mirror_index(f,j); //f2->index(f);
	set_adjacency(f1, i1, f2, i2);
	delete_face(f);
      }
    
      v->set_face( *(faces_list.begin()));
    }
    break;
  default:
    CGAL_triangulation_assertion(false);
    break;  }
  return v;
}


template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
remove_degree_3(Vertex_handle v, Face_handle f)
// remove a vertex of degree 3
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(degree(v) == 3);

  if (f == Face_handle()) {f= v->face();}
  else { CGAL_triangulation_assertion( f->has_vertex(v));}
      
  int i = f->index(v);
  Face_handle left = f->neighbor(cw(i));
  int li = mirror_index(f,cw(i)); 
  Face_handle right = f->neighbor(ccw(i));
  int ri = mirror_index(f,ccw(i)); 

  Face_handle ll, rr;
  Vertex_handle q = left->vertex(li);
  CGAL_triangulation_assertion( left->vertex(li) == right->vertex(ri));
    
  ll = left->neighbor(cw(li));
  if(ll != Face_handle()) {
    int lli = mirror_index(left,cw(li)); 
    ll->set_neighbor(lli, f);
  } 
  f->set_neighbor(cw(i), ll);
  if (f->vertex(ccw(i))->face() == left) f->vertex(ccw(i))->set_face(f);    
        
  rr = right->neighbor(ccw(ri));
  if(rr != Face_handle()) {
    int rri =  mirror_index(right,ccw(ri)); //rr->index(right);
    rr->set_neighbor(rri, f);
  } 
  f->set_neighbor(ccw(i), rr);
  if (f->vertex(cw(i))->face() == right) f->vertex(cw(i))->set_face(f);  
        
  f->set_vertex(i, q);
  if (q->face() == right || q->face() == left) {
    q->set_face(f);
  }
  delete_face(right);
  delete_face(left);
        
  delete_vertex(v);
} 

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
dim_down(Face_handle f, int i)
{
  CGAL_triangulation_expensive_precondition( is_valid() );
  CGAL_triangulation_precondition( dimension() == 2 );
  CGAL_triangulation_precondition( number_of_vertices() > 3 );
  CGAL_triangulation_precondition( degree( f->vertex(i) ) == 
                                   number_of_vertices()-1 );

  Vertex_handle v = f->vertex(i);
  std::list<Face_handle > to_delete;
  std::list<Face_handle> to_downgrade;
  Face_iterator ib = face_iterator_base_begin();
  for( ; ib != face_iterator_base_end(); ++ib ){
    if ( ! ib->has_vertex(v) ) { to_delete.push_back(ib);}
    else { to_downgrade.push_back(ib);}
  }

  typename std::list<Face_handle>::iterator lfit = to_downgrade.begin();
  int j;
  for( ; lfit !=  to_downgrade.end() ; ++lfit) {
    Face_handle fs = *lfit; j = fs->index(v);
    if (j == 0) fs->cw_permute();
    else if(j == 1) fs->ccw_permute();
    fs->set_vertex(2, Vertex_handle());
    fs->set_neighbor(2, Face_handle());
    fs->vertex(0)->set_face(fs);
  }
  lfit = to_delete.begin();
  for( ; lfit !=  to_delete.end() ; ++lfit) {
    delete_face(*lfit);
  }
  set_dimension(dimension() -1);
  Face_handle n0 = f->neighbor(0);
  //Face_handle n1 = f->neighbor(1);
  //Vertex_handle v0 = f->vertex(0);
  Vertex_handle v1 = f->vertex(1);
  f->set_vertex(1, v);
  Face_handle fl = create_face(v, v1, Vertex_handle(),
	                       n0, f, Face_handle());
  f->set_neighbor(0, fl);
  n0->set_neighbor(1, fl);
  v->set_face(f);
}
  
template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
remove_dim_down(Vertex_handle v)
{
  Face_handle f;
  switch( dimension()){
  case -1: 
    delete_face(v->face());
    break;
  case 0:
    f = v->face();
    f->neighbor(0)->set_neighbor(0,Face_handle());
    delete_face(v->face());
    break;
  case 1:
  case 2:
//  CGAL_triangulation_precondition ( 
//           (dimension() == 1 &&  number_of_vertices() == 3) ||
//           (dimension() == 2 && number_of_vertices() > 3) );
    // the faces incident to v are down graded one dimension
    // the other faces are deleted
    std::list<Face_handle > to_delete;
    std::list<Face_handle > to_downgrade;
    Face_iterator ib = face_iterator_base_begin();
    for( ; ib != face_iterator_base_end(); ++ib ){
      if ( ! ib->has_vertex(v) ) { to_delete.push_back(ib);}
      else { to_downgrade.push_back(ib);}
    }

    typename std::list<Face_handle>::iterator lfit = to_downgrade.begin();
    int j;
    for( ; lfit !=  to_downgrade.end() ; ++lfit) {
      f = *lfit; j = f->index(v);
      if (dimension() == 1) {
	if (j == 0) 	f->reorient();
	f->set_vertex(1,Vertex_handle());
	f->set_neighbor(1, Face_handle());
      }
      else { //dimension() == 2
	if (j == 0) f->cw_permute();
	else if(j == 1) f->ccw_permute();
	f->set_vertex(2, Vertex_handle());
	f->set_neighbor(2, Face_handle());
      }
      f->vertex(0)->set_face(f);
    }

    lfit = to_delete.begin();
    for( ; lfit !=  to_delete.end() ; ++lfit) {
      delete_face(*lfit);
    }
  }  
  delete_vertex(v);
  set_dimension(dimension() -1);
  return;
}

template <  class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::  
remove_1D(Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 1 &&
				   number_of_vertices() > 3);
  Face_handle f = v->face();
  int i = f->index(v);
  if (i==0) {f = f->neighbor(1);}
  CGAL_triangulation_assertion( f->index(v) == 1);
  Face_handle g= f->neighbor(0);
  f->set_vertex(1, g->vertex(1));
  set_adjacency(f, 0, g->neighbor(0), 1);
  g->vertex(1)->set_face(f);
  delete_face(g);
  delete_vertex(v);
  return;
}



template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
remove_second(Vertex_handle v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 2 &&
 				  dimension() == 0);
  remove_dim_down(v);
  return;
}

    
template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
remove_first(Vertex_handle v)
{
  CGAL_triangulation_precondition(number_of_vertices()== 1 && 
 				  dimension() == -1);
  remove_dim_down(v);
  return; 
}

template <class Vb, class Fb>
inline
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
star_hole(List_edges& hole)
{
  Vertex_handle newv = create_vertex();
  star_hole(newv, hole);
  return newv;
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
star_hole(Vertex_handle newv, List_edges& hole)
  // star the hole represented by hole around newv
  // the triangulation is assumed to have dim=2
  // hole is supposed to be ccw oriented
{
   
  star_hole(newv, hole.begin(), hole.end());
  return;	    
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
make_hole(Vertex_handle v, List_edges& hole)
  // delete the faces incident to v and v
  // and return the dscription of the hole in hole
{
 CGAL_triangulation_precondition(dimension() == 2);
 std::list<Face_handle> to_delete;  

 Face_handle  f, fn;
 int i =0, in =0;
 Vertex_handle  vv;

 Face_circulator fc = incident_faces(v);
 Face_circulator done(fc);
 do {
   f = fc ;
   i = f->index(v);
   fn = f->neighbor(i);
   in = mirror_index(f,i); //fn->index(f);
   vv = f->vertex(cw(i));
   if( vv->face()==  f) vv->set_face(fn);
   vv = fc->vertex(ccw(i));
   if( vv->face()== f) vv->set_face(fn);
   fn->set_neighbor(in, Face_handle());
   hole.push_back(Edge(fn,in));
   to_delete.push_back(f);
 }
  while(++fc != done);

  while (! to_delete.empty()){
    delete_face(to_delete.front());
    to_delete.pop_front();
  }
  delete_vertex(v);
  return;
}


template <class Vb, class Fb>
inline
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
create_vertex(const Vertex &v)
{
  return vertices().insert(v);
}

template <class Vb, class Fb>
inline
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
create_vertex(Vertex_handle vh)
{
  return vertices().insert(*vh);
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(const Face& f)
{
  return faces().insert(f);
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face( Face_handle fh)
{
  return create_face(*fh);
}


template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, 
	    Face_handle f2, int i2, 
	    Face_handle f3, int i3)
{
  Face_handle newf = faces().emplace(f1->vertex(cw(i1)),
					      f2->vertex(cw(i2)),
					      f3->vertex(cw(i3)),
					      f2, f3, f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  f3->set_neighbor(i3,newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, Face_handle f2, int i2)
{
  Face_handle newf = faces().emplace(f1->vertex(cw(i1)),
					      f2->vertex(cw(i2)),
					      f2->vertex(ccw(i2)),
					      f2, Face_handle(), f1);
  f1->set_neighbor(i1,newf);
  f2->set_neighbor(i2,newf);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Face_handle f1, int i1, Vertex_handle v)
{
  Face_handle newf = create_face();
  newf->set_vertices(f1->vertex(cw(i1)), f1->vertex(ccw(i1)), v);
  newf->set_neighbors(Face_handle(), Face_handle(), f1);
  f1->set_neighbor(i1,newf);
  return newf;
}


template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  Face_handle newf = faces().emplace(v1, v2, v3);
  return newf;
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Face_handle
Triangulation_data_structure_2<Vb,Fb>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
	    Face_handle f1, Face_handle f2, Face_handle f3)
{
  Face_handle newf = faces().emplace(v1, v2, v3, f1, f2, f3);

  return(newf);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
set_adjacency(Face_handle f0, int i0, Face_handle f1, int i1) const
{
  CGAL_triangulation_assertion(i0 >= 0 && i0 <= dimension());
  CGAL_triangulation_assertion(i1 >= 0 && i1 <= dimension());
  CGAL_triangulation_assertion(f0 != f1);
  f0->set_neighbor(i0,f1);
  f1->set_neighbor(i1,f0);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
delete_face(Face_handle f)
{
  CGAL_triangulation_expensive_precondition( dimension() != 2 || is_face(f));
  CGAL_triangulation_expensive_precondition( dimension() != 1 || is_edge(f,2));
  CGAL_triangulation_expensive_precondition( dimension() != 0 ||
					     is_vertex(f->vertex(0)) );
  faces().erase(f);
}

template <class Vb, class Fb>
inline void
Triangulation_data_structure_2<Vb,Fb>::
delete_vertex(Vertex_handle v)
{
  CGAL_triangulation_expensive_precondition( is_vertex(v) );
  vertices().erase(v);
}

// split and join operations

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Fourtuple
Triangulation_data_structure_2<Vb,Fb>::
split_vertex(Vertex_handle v, Face_handle f1, Face_handle g1)
{
  /*
  // The following method preforms a split operation of the vertex v
  // using the faces f1 and g1. The split operation is shown
  // below.
  // The names of the variables in the method correspond to the
  // quantities in the drawings below
  //
  // The configuration before the split:
  //
  //                  cw(i1)   v3   ccw(i2)
  //                     *-----*-----*
  //                    / \    |    / \
  //                   /   \ f1|f2 /   \
  //                  /     \  |  /     \
  //                 /       \ | /       \
  //                /         \|/v        \
  //               *-----------*-----------*
  //                \         /|\         /
  //                 \       / | \       /
  //                  \     /  |  \     /
  //                   \   / g2|g1 \   /
  //                    \ /    |    \ /
  //                     *-----*-----*
  //                 ccw(j2)   v4   cw(j1)
  //
  //
  // The configuration after the split:
  //
  //
  //               cw(i1)      v3     ccw(i2)
  //                 *---------*---------*
  //                / \       / \       / \
  //               /   \  f1 /   \  f2 /   \
  //              /     \   /  f  \   /     \
  //             /       \ /     v2\ /       \
  //            *---------*---------*---------*
  //             \       / \v1     / \       /
  //              \     /   \  g  /   \     /
  //               \   /  g2 \   /  g1 \   /
  //                \ /       \ /       \ /
  //                 *---------*---------*
  //              ccw(j2)      v4      cw(j1)
  //
  */

  CGAL_triangulation_expensive_precondition( is_valid() );

  CGAL_triangulation_precondition( dimension() == 2 );
  CGAL_triangulation_precondition( f1 != Face_handle() && f1->has_vertex(v) );
  CGAL_triangulation_precondition( g1 != Face_handle() && g1->has_vertex(v) );

  // 1. first we read some information that we will need
  int i1 = f1->index(v);
  int j1 = g1->index(v);
  Face_handle f2 = f1->neighbor( cw(i1) );
  Face_handle g2 = g1->neighbor( cw(j1) );

  int i2 = f2->index(v);
  int j2 = g2->index(v);

  Vertex_handle v3 = f1->vertex( ccw(i1) );
  Vertex_handle v4 = g1->vertex( ccw(j1) );

  // lst is the list of faces adjecent to v stored in
  // counterclockwise order from g2 to f1) inclusive.
  // the list idx contains the indices of v in the
  // faces in lst.
  std::list<Face_handle> lst;
  std::list<int>         idx;

  Face_circulator fc(v, g1);
  Face_handle ff(fc);
  while ( ff != f2 ) {
    lst.push_back( ff );
    idx.push_back( ff->index(v) );
    fc++;
    ff = Face_handle(fc);
  }
  lst.push_back( ff );
  idx.push_back( ff->index(v) );

  // 2. we create the new vertices and the two new faces
  Vertex_handle v1 = v;
  Vertex_handle v2 = create_vertex();
  Face_handle f = create_face(v1, v2, v3);
  Face_handle g = create_face(v2, v1, v4);

  // 3. we update the adjacency information for the new vertices and
  //    the new faces
  f->set_neighbor(0, f2);
  f->set_neighbor(1, f1);
  f->set_neighbor(2, g);
  g->set_neighbor(0, g2);
  g->set_neighbor(1, g1);
  g->set_neighbor(2, f);
  v1->set_face(f);
  v2->set_face(g);

  // 4. update the vertex for the faces f2 through g1 in
  //    counterclockwise order
  typename std::list<Face_handle>::iterator fit = lst.begin();
  typename std::list<int>::iterator         iit = idx.begin();
  for (; fit != lst.end(); ++fit, ++iit) {
    (*fit)->set_vertex(*iit, v2);
  }

  lst.clear();
  idx.clear();

  // 5. make f and g the new neighbors of f1, f2 and g1, g2
  //    respectively.
  f1->set_neighbor(  cw(i1), f );
  f2->set_neighbor( ccw(i2), f );
  g1->set_neighbor(  cw(j1), g );
  g2->set_neighbor( ccw(j2), g );

  CGAL_triangulation_expensive_postcondition( is_valid() );

  // 6. return the new stuff
  return Fourtuple(v1, v2, f, g);
}

template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
join_vertices(Face_handle f, int i, Vertex_handle v)
{
  CGAL_triangulation_expensive_precondition( is_valid() );
  CGAL_triangulation_precondition( f != Face_handle() );
  CGAL_triangulation_precondition( i >= 0 && i <= 2 );

  // this methods does the "join"-operation and preserves
  // the vertex v among the two vertices that define the edge (f, i) 

  Vertex_handle v1 = f->vertex( ccw(i) );
  Vertex_handle v2 = f->vertex( cw(i)  );

  CGAL_triangulation_precondition( v == v1 || v == v2 );

  if ( v == v2 ) {
    return join_vertices(f->neighbor(i), mirror_index(f,i), v);
  }

  size_type deg2 = degree(v2);

  CGAL_triangulation_precondition( deg2 >= 3 );

  if ( deg2 == 3 ) {
    remove_degree_3(v2, f->neighbor(ccw(i)));
    return v1;
  }
  
  /*
  // The following drawing corrsponds to the variables
  // used in this part...
  // The vertex v1 is returned...
  //
  //      itl       i=v3      itr
  //       *---------*---------*
  //        \       / \       /
  //         \  tl /   \  tr /
  //          \   /  f  \   /
  //           \ /       \ /
  //  v1=ccw(i) *---------*  cw(i)=v2
  //           / \       / \
  //          /   \  g  /   \
  //         /  bl \   /  br \
  //        /       \ /	      \
  //       *---------*---------*
  //      ibl       j=v4      ibr
  //                                                           
  // The situation after the "join"-operation is as follows:
  //
  //                 i
  //           *-----*-----*
  //            \    |    /
  //             \ tl|tr /
  //              \  |  /
  //               \ | /
  //                \|/
  //                 *  v1
  //                /|\
  //               / | \
  //              /  |	\
  //             / bl|br \
  //            /    |	  \
  //           *-----*-----*
  //
  */

  // first we register all the needed info
  Face_handle g = f->neighbor(i);
  int j = mirror_index(f,i);

  Face_handle tl = f->neighbor( cw(i)  );
  Face_handle tr = f->neighbor( ccw(i) );

  int itl = mirror_index(f, cw(i)  );
  int itr = mirror_index(f, ccw(i) );

  Face_handle bl = g->neighbor( ccw(j) );
  Face_handle br = g->neighbor( cw(j)  );

  int ibl = mirror_index(g, ccw(j) );
  int ibr = mirror_index(g, cw(j)  );

  // we need to store the faces adjacent to v2 as well as the
  // indices of v2 w.r.t. these faces, so that afterwards we can set 
  // v1 to be the vertex for these faces
  std::vector<Face_handle> star_faces_of_v2;
  std::vector<int> star_indices_of_v2;
  Face_circulator fc_start(v2);
  Face_circulator fc = fc_start;

  do {
    Face_handle ff(fc);
    star_faces_of_v2.push_back(ff);
    star_indices_of_v2.push_back(ff->index(v2));
    ++fc;
  } while ( fc != fc_start );

  CGAL_triangulation_assertion(
    static_cast<size_type>(star_faces_of_v2.size()) == deg2 );

  // from this point and on we modify the values

  // first set the neighbors
  set_adjacency(tl, itl, tr, itr);
  set_adjacency(bl, ibl, br, ibr);

  // make sure that all the faces containing v2 as a vertex, now
  // contain v1
  for (unsigned int k = 0; k < star_faces_of_v2.size(); k++) {
    int id = star_indices_of_v2[k];
    CGAL_triangulation_assertion( star_faces_of_v2[k]->vertex(id) == v2 );
    star_faces_of_v2[k]->set_vertex( id, v1 );
  }

  // then make sure that all the vertices have correct pointers to 
  // faces
  Vertex_handle v3 = f->vertex(i);
  Vertex_handle v4 = g->vertex(j);
  if ( v3->face() == f )  v3->set_face(tr);
  if ( v4->face() == g )  v4->set_face(br);
  if ( v1->face() == f || v1->face() == g ) v1->set_face(tl);


#if ! defined(CGAL_TRIANGULATION_NO_ASSERTIONS) && ! defined(CGAL_NO_ASSERTIONS)
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit) {
    int id;
    CGAL_triangulation_assertion( !fit->has_vertex(v2, id) );
  }
#endif

  // memory management
  star_faces_of_v2.clear();
  star_indices_of_v2.clear();

  delete_face(f);
  delete_face(g);

  delete_vertex(v2);

  CGAL_triangulation_expensive_postcondition( is_valid() );

  return v1;
}

// insert_degree_2 and remove_degree_2 operations
template <class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
insert_degree_2(Face_handle f, int i)
{
  /*
  // This method basically does the following transformation
  // The remove_degree_2 method performs the same operation in the
  // opposite direction
  //
  //
  //                                                *
  //                 i                             / \
  //                 *                            /   \
  //                / \                          /  f  \
  //               /   \                        / _____	\
  //              /  f  \                      / /  f1 \ \
  //             /       \                     |/   v   \|
  //  v0=ccw(i) *---------* v1=cw(i)  ===>  v0 *----*----* v1
  //             \       /                     |\   f2  /|
  //              \  g  /                      \ \_____/ /
  //               \   /                        \       /
  //                \ /                          \  g  /
  //                 *                            \   /
  //                 j                             \ /
  //                                                *
  //
  */

  Face_handle g = f->neighbor(i);
  int j = mirror_index(f,i);

  Vertex_handle  v = create_vertex();

  Vertex_handle v0 = f->vertex( ccw(i) );
  Vertex_handle v1 = f->vertex( cw(i)  );

  Face_handle f_undef;

  Face_handle f1 = create_face(v0, v, v1, f_undef, f, f_undef);
  Face_handle f2 = create_face(v0, v1, v, f_undef, f_undef, g);

  set_adjacency(f1, 0, f2, 0);
  set_adjacency(f1, 2, f2, 1);

  f->set_neighbor(i, f1);
  g->set_neighbor(j, f2);

  v->set_face(f1);

  return v;
}

template <class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
remove_degree_2(Vertex_handle v)
{
  CGAL_precondition( degree(v) == 2 );

  Face_handle f1 = v->face();
  int i = f1->index(v);

  Face_handle f2 = f1->neighbor( ccw(i) );
  int j = f2->index(v);

  Face_handle ff1 = f1->neighbor( i );
  Face_handle ff2 = f2->neighbor( j );

  int id1 = mirror_index(f1,i);
  int id2 = mirror_index(f2,j);

  set_adjacency(ff1, id1, ff2, id2);

  Vertex_handle v1 = f1->vertex( ccw(i) );
  //    if ( v1->face() == f1 || v1->face() == f2 ) {
  v1->set_face(ff1);
  //    }

  Vertex_handle v2 = f1->vertex( cw(i) );
  //    if ( v2->face() == f1 || v2->face() == f2 ) {
  v2->set_face(ff2);
  //    }

  delete_face(f1);
  delete_face(f2);

  delete_vertex(v);
}

// CHECKING
template <  class Vb, class Fb>
bool
Triangulation_data_structure_2<Vb,Fb>::
is_valid(bool verbose, int level) const
{
  if(number_of_vertices() == 0){ 
    return (dimension() == -2);
  }

      
  bool result = (dimension()>= -1);
  CGAL_triangulation_assertion(result);

  //count and test the validity of the faces (for positive dimensions)
  Face_iterator ib = face_iterator_base_begin(); 
  Face_iterator ib_end = face_iterator_base_end();
  size_type count_stored_faces =0;
  for ( ; ib != ib_end ; ++ib){
    count_stored_faces += 1;
    if (dimension()>= 0) {
      result = result && ib->is_valid(verbose,level);
      CGAL_triangulation_assertion(result);
    }
  }
  
  result = result && (count_stored_faces == number_of_full_dim_faces());
  CGAL_triangulation_assertion(
		 count_stored_faces == number_of_full_dim_faces());
 
  // vertex count
  size_type vertex_count = 0;
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end();
      ++vit) {
    CGAL_triangulation_assertion( vit->face() != Face_handle());
    result = result && vit->is_valid(verbose,level);
    CGAL_triangulation_assertion( result );
    ++vertex_count;
  }
  result = result && (number_of_vertices() == vertex_count);
  CGAL_triangulation_assertion( number_of_vertices() == vertex_count );
    
  //edge count
  size_type edge_count = 0;
  for(Edge_iterator eit = edges_begin(); eit != edges_end(); ++eit) { 
    ++edge_count;
  }

  // face count
  size_type face_count = 0;
  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit) {
    ++face_count;
  }
        
  switch(dimension()) {
  case -1: 
    result = result && vertex_count == 1 && face_count == 0
      && edge_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 0:
    result = result && vertex_count == 2 && face_count == 0
      && edge_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 1:
    result = result &&  edge_count == vertex_count;
    CGAL_triangulation_assertion(result);
    result = result &&  face_count == 0;
    CGAL_triangulation_assertion(result);
    break;
  case 2:
    result = result &&  edge_count == 3*face_count/2 ;
    CGAL_triangulation_assertion(edge_count == 3*face_count/2);
    break;
  default:
    result = false;
    CGAL_triangulation_assertion(result);
  }
  return result;
}

template <class Vb, class Fb>
template <class TDS_src,class ConvertVertex,class ConvertFace>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
copy_tds(const TDS_src& tds_src,
        typename TDS_src::Vertex_handle vert,
        const ConvertVertex& convert_vertex,
        const ConvertFace& convert_face)
{
  if (vert != typename TDS_src::Vertex_handle()) 
    CGAL_triangulation_precondition( tds_src.is_vertex(vert));

  clear();
  size_type n = tds_src.number_of_vertices();
  set_dimension(tds_src.dimension());

  // Number of pointers to cell/vertex to copy per cell.
  int dim = (std::max)(1, dimension() + 1);
 
  if(n == 0) {return Vertex_handle();}
  
  //initializes maps
  Unique_hash_map<typename TDS_src::Vertex_handle,Vertex_handle> vmap;
  Unique_hash_map<typename TDS_src::Face_handle,Face_handle> fmap;

  // create vertices
  typename TDS_src::Vertex_iterator vit1 = tds_src.vertices_begin();
  for( ; vit1 != tds_src.vertices_end(); ++vit1) {
    Vertex_handle vh = create_vertex( convert_vertex(*vit1) );
    vmap[vit1] = vh;
    convert_vertex(*vit1, *vh);
  }

  //create faces 
  typename TDS_src::Face_iterator fit1 = tds_src.faces().begin();
  for( ; fit1 != tds_src.faces_end(); ++fit1) {
    Face_handle fh = create_face( convert_face(*fit1) );
    fmap[fit1] = fh;
    convert_face(*fit1, *fh);
  }

  //link vertices to a cell 
  vit1 = tds_src.vertices_begin();
  for ( ; vit1 != tds_src.vertices_end(); vit1++) {
    vmap[vit1]->set_face(fmap[vit1->face()]);
  }

  //update vertices and neighbor pointers
  fit1 = tds_src.faces().begin();
  for ( ; fit1 != tds_src.faces_end(); ++fit1) {
      for (int j = 0; j < dim ; ++j) {
	fmap[fit1]->set_vertex(j, vmap[fit1->vertex(j)] );
	fmap[fit1]->set_neighbor(j, fmap[fit1->neighbor(j)]);
      }
    }
   
  // remove the post condition because it is false when copying the
  // TDS of a regular triangulation because of hidden vertices
  // CGAL_triangulation_postcondition( is_valid() );
  return (vert == typename TDS_src::Vertex_handle())  ? Vertex_handle() : vmap[vert];
}

//utilities for copy_tds
namespace internal { namespace TDS_2{
  template <class Vertex_src,class Vertex_tgt>
  struct Default_vertex_converter
  {
    Vertex_tgt operator()(const Vertex_src& src) const {
      return Vertex_src( src.point() );
    }
    
    void operator()(const Vertex_src&,Vertex_tgt&) const {}
  };

  template <class Face_src,class Face_tgt>
  struct Default_face_converter
  {
    Face_tgt operator()(const Face_src& /*src*/) const {
      return Face_tgt();
    } 
    
    void operator()(const Face_src&,Face_tgt&) const {}
  };
  
  template <class Vertex>
  struct Default_vertex_converter<Vertex,Vertex>
  {
    const Vertex& operator()(const Vertex& src) const {
      return src;
    }
    
    void operator()(const Vertex&,Vertex&) const {}
  };
  
  template <class Face>
  struct Default_face_converter<Face,Face>{
    const Face& operator()(const Face& src) const {
      return src;
    } 
    
    void operator()(const Face&,Face&) const {}
  };
} } //namespace internal::TDS_2

template <  class Vb, class Fb>
template < class TDS_src>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
copy_tds(const TDS_src &src, typename TDS_src::Vertex_handle vh)
  // return the vertex corresponding to vh in the new tds
{
  if (this == &src) return Vertex_handle();
  internal::TDS_2::Default_vertex_converter<typename TDS_src::Vertex,Vertex> setv;
  internal::TDS_2::Default_face_converter<typename TDS_src::Face,Face>  setf;
  return copy_tds(src,vh,setv,setf);
}

template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
file_output( std::ostream& os, Vertex_handle v, bool skip_first) const
{
  // ouput to a file
  // if non NULL, v is the vertex to be output first
  // if skip_first is true, the point in the first vertex is not output
  // (it may be for instance the infinite vertex of the triangulation)
  
  size_type n = number_of_vertices();
  size_type m = number_of_full_dim_faces();
  if(is_ascii(os))  os << n << ' ' << m << ' ' << dimension() << std::endl;
  else     os << n << m << dimension();
  if (n==0) return;

  Unique_hash_map<Vertex_handle,int> V;
  Unique_hash_map<Face_handle,int> F;


  // first vertex 
  int inum = 0;
  if ( v != Vertex_handle()) {
    V[v] = inum++;
    if( ! skip_first){
      // os << v->point();
      os << *v ;
    if(is_ascii(os))  os << std::endl;
    }
  }
  
  // other vertices
  for( Vertex_iterator vit= vertices_begin(); vit != vertices_end() ; ++vit) {
    if ( v != vit ) {
	V[vit] = inum++;
	// os << vit->point();
	os << *vit;
	if(is_ascii(os)) os << "\n";
    }
  }
  if(is_ascii(os)) os << "\n";

  // vertices of the faces
  inum = 0;
  int dim = (dimension() == -1 ? 1 :  dimension() + 1);
  for( Face_iterator ib = face_iterator_base_begin();
       ib != face_iterator_base_end(); ++ib) {
    F[ib] = inum++;
    for(int j = 0; j < dim ; ++j) {
      os << V[ib->vertex(j)];
      if(is_ascii(os)) os << " ";
    }
    os << *ib ;
    if(is_ascii(os)) os << "\n";
  }
  if(is_ascii(os)) os << "\n";
    
  // neighbor pointers of the  faces
  for( Face_iterator it = face_iterator_base_begin();
       it != face_iterator_base_end(); ++it) {
    for(int j = 0; j < dimension()+1; ++j){
      os << F[it->neighbor(j)];
      if(is_ascii(os))  os << " ";
    }
    if(is_ascii(os)) os << "\n";
  }

  return ;
}


template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
file_input( std::istream& is, bool skip_first)
{
  //input from file
  //return a pointer to the first input vertex
  // if skip_first is true, a first vertex is added (infinite_vertex)
  //set this  first vertex as infinite_Vertex
  if(number_of_vertices() != 0)    clear();
  
  size_type n, m;
  int d;
  is >> n >> m >> d;

  if (n==0){ return Vertex_handle();}

  set_dimension(d);

  std::vector<Vertex_handle > V(n);
  std::vector<Face_handle> F(m);

  // read vertices
  size_type i = 0;
  if(skip_first){
    V[0] = create_vertex();
    ++i;
  }
  for( ; i < n; ++i) {
    V[i] = create_vertex();
    is >> *(V[i]);
  }
  
  // Creation of the faces
  int index;
  int dim = (dimension() == -1 ? 1 :  dimension() + 1);
  {
    for(i = 0; i < m; ++i) {
      F[i] = create_face() ;
      for(int j = 0; j < dim ; ++j){
	is >> index;
	F[i]->set_vertex(j, V[index]);
	// The face pointer of vertices is set too often,
	// but otherwise we had to use a further map
	V[index]->set_face(F[i]);
      }
      // read in non combinatorial info of the face
      is >> *(F[i]) ;
    }
  }

  // Setting the neighbor pointers 
  {
    for(i = 0; i < m; ++i) {
      for(int j = 0; j < dimension()+1; ++j){
	is >> index;
	F[i]->set_neighbor(j, F[index]);
      }
    }
  }
  
  return V[0];
}


template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
vrml_output( std::ostream& os, Vertex_handle v, bool skip_infinite) const
{
  // ouput to a vrml file style
  // Point are assumed to be 3d points with a stream oprator <<
  // if non NULL, v is the vertex to be output first
  // if skip_inf is true, the point in the first vertex is not output
  // and the faces incident to v are not output
  // (it may be for instance the infinite vertex of the terrain)
  os << "#VRML V2.0 utf8" << std::endl;
  os << "Shape {" << std::endl;
  os << "\tgeometry IndexedFaceSet {" << std::endl;
  os << "\t\tcoord Coordinate {" << std::endl;
  os << "\t\t\tpoint [" << std::endl;

  Unique_hash_map<Vertex_handle,int> vmap;

  Vertex_iterator vit;
  Face_iterator fit;

  //first vertex
  int inum = 0;
  if ( v != Vertex_handle()) {
    vmap[v] = inum++;
    if( ! skip_infinite)  os << "\t\t\t\t" << *v << std::endl;
  }

  //other vertices
  for( vit= vertices_begin(); vit != vertices_end() ; ++vit) {
    if ( v != vit) {
      vmap[vit] = inum++;
      os << "\t\t\t\t" << *vit << std::endl;
    }
  }

   os << "\t\t\t]" << std::endl;
   os << "\t\t}" << std::endl;
   os << "\t\tcoordIndex [" << std::endl;

   // faces
   for(fit= faces_begin(); fit != faces_end(); ++fit) {
     if (!skip_infinite || !fit->has_vertex(v)) {
   	os << "\t\t\t";
	os << vmap[(*fit).vertex(0)] << ", ";
	os << vmap[(*fit).vertex(1)] << ", ";
	os << vmap[(*fit).vertex(2)] << ", ";
	os << "-1, " << std::endl;  
     }
   }
   os << "\t\t]" << std::endl;
   os << "\t}" << std::endl;
   os << "}" << std::endl;
   return;
}

template < class Vb, class Fb>
typename Triangulation_data_structure_2<Vb,Fb>::Vertex_handle
Triangulation_data_structure_2<Vb,Fb>::
off_file_input( std::istream& is, bool verbose)
{
  // input from an OFF file
  // assume a dimension 2 triangulation
  // create an infinite-vertex and  infinite faces with the
  // boundary edges if any.
  // return the infinite vertex if created
  Vertex_handle vinf;
  File_scanner_OFF scanner(is, verbose);
  if (! is) {
    if (scanner.verbose()) {
         std::cerr << " " << std::endl;
	 std::cerr << "TDS::off_file_input" << std::endl;
	 std::cerr << " input error: file format is not OFF." << std::endl;
    }
    return vinf;
  }

  if(number_of_vertices() != 0)    clear();
  int dim = 2;
  set_dimension(dim);

  std::vector<Vertex_handle > vvh(scanner.size_of_vertices());
  std::map<Vh_pair, Edge> edge_map;
  typedef typename Vb::Point   Point;

  // read vertices
  std::size_t i;
  for ( i = 0; i < scanner.size_of_vertices(); i++) {
    Point p;
    file_scan_vertex( scanner, p);
    vvh[i] = create_vertex();
    vvh[i]->set_point(p);
    scanner.skip_to_next_vertex( i);
  }
  if ( ! is ) {
    is.clear( std::ios::badbit);
    return vinf;
  }
  //vinf = vvh[0];

  // create the facets
  for ( i = 0; i < scanner.size_of_facets(); i++) {
    Face_handle fh = create_face();
    std::size_t no;
    scanner.scan_facet( no, i);
    if( ! is || no != 3) {
      if ( scanner.verbose()) {
	std::cerr << " " << std::endl;
	std::cerr << "TDS::off_file_input" << std::endl;
	std::cerr << "facet " << i << "does not have  3 vertices." 
		  << std::endl;
      }
      is.clear( std::ios::badbit);
      return vinf;
    }

    for ( std::size_t j = 0; j < no; ++j) {
      std::size_t index;
      scanner.scan_facet_vertex_index( index, i);
      fh->set_vertex(j, vvh[index]);
      vvh[index]->set_face(fh);
    }

    for (std::size_t ih  = 0; ih < no; ++ih) {
	set_adjacency(fh, ih, edge_map);
    }
  }

  // deal with  boundaries
  if ( !edge_map.empty()) {
    vinf = create_vertex();
    std::map<Vh_pair, Edge> inf_edge_map;
   while (!edge_map.empty()) {
     Face_handle fh = edge_map.begin()->second.first;
     int ih = edge_map.begin()->second.second;
     Face_handle fn = create_face( vinf, 
				   fh->vertex(cw(ih)), 
				   fh->vertex(ccw(ih)));
     vinf->set_face(fn);
     set_adjacency(fn, 0, fh, ih);
     set_adjacency(fn, 1, inf_edge_map);
     set_adjacency(fn, 2, inf_edge_map);
     edge_map.erase(edge_map.begin());
   }
   CGAL_triangulation_assertion(inf_edge_map.empty());
  }
  
  
  // coherent orientation
  reorient_faces();
  return vinf;
}


template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
set_adjacency(Face_handle fh, 
	      int ih, 
	      std::map< Vh_pair, Edge>& edge_map)
{
  // set adjacency to (fh,ih) using the the map edge_map
  // or insert (fh,ih) in edge map
  Vertex_handle vhcw  =  fh->vertex(cw(ih));
  Vertex_handle vhccw =  fh->vertex(ccw(ih)); 
  Vh_pair  vhp =  vhcw < vhccw ?  
                  std::make_pair(vhcw, vhccw) 
                : std::make_pair(vhccw, vhcw) ;
  typename std::map<Vh_pair, Edge>::iterator emapit = edge_map.find(vhp);
  if (emapit == edge_map.end()) {// not found, insert edge
    edge_map.insert(std::make_pair(vhp, Edge(fh,ih)));
  }
  else { //found set adjacency and erase
    Edge e = emapit->second;
    set_adjacency( fh,ih, e.first, e.second);
    edge_map.erase(emapit);
  } 
}



template < class Vb, class Fb>
void
Triangulation_data_structure_2<Vb,Fb>::
reorient_faces()
{
  // reorient the faces of a triangulation 
  // needed for example in off_file_input
  // because the genus is not known, the number of faces 
  std::set<Face_handle> oriented_set;
  std::stack<Face_handle>  st;
  Face_iterator fit = faces_begin();
  std::ptrdiff_t nf  = std::distance(faces_begin(),faces_end());

  while (0 != nf) {
    while ( !oriented_set.insert(fit).second ){
      ++fit; // find a germ for  non oriented components 
    }
    // orient component
    --nf;
    st.push(fit);
    while ( ! st.empty()) {
      Face_handle fh = st.top();
      st.pop();
      for(int ih = 0 ; ih < 3 ; ++ih){
	Face_handle fn = fh->neighbor(ih);
	if (oriented_set.insert(fn).second){
	  int in = fn->index(fh);
	  if (fn->vertex(cw(in)) != fh->vertex(ccw(ih))) fn->reorient();
          --nf;
	  st.push(fn);
	}
      }
    }

  }
  return;
}
	  

template <  class Vb, class Fb>
std::istream&
operator>>(std::istream& is,  
	   Triangulation_data_structure_2<Vb,Fb>& tds) 
{
  tds.file_input(is);
  return is;
}


template <  class Vb, class Fb>
std::ostream&
operator<<(std::ostream& os, 
	   const Triangulation_data_structure_2<Vb,Fb>  &tds) 
{
   tds.file_output(os);
   return os;
}


} //namespace CGAL 

#endif //CGAL_TRIANGULATION_DATA_STRUCTURE_2_H
