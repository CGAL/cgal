// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Andreas Fabri, Mariette Yvinec

#ifndef CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Polyline_constraint_hierarchy_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

// Tr the base triangulation class 
// Tr has to be Constrained or Constrained_Delaunay

template < class Tr >
class Polyline_constrained_triangulation_2  
  : public Tr
{


template<class CDT>
class Face_container
{
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Face_handle Face_handle;
private:
  typedef boost::tuple<Vertex_handle, Vertex_handle, Vertex_handle> TFace; 
  std::list<TFace> faces;
  bool init;
  CDT& cdt;

public:
  typedef Face_handle value_type;
  typedef Face_handle& reference;
  typedef const Face_handle& const_reference;

  Face_container(CDT& cdt_ ) : init(false), cdt(cdt_) {}

  void push_back(Face_handle fh)
  {
    faces.push_back(boost::make_tuple(fh->vertex(0),
                                      fh->vertex(1),
                                      fh->vertex(2)));

  }

  template <class OutputIterator>
  void
  write_faces(OutputIterator out)
  {
    while(! faces.empty()){
      TFace tf = faces.back();
      faces.pop_back();
      Face_handle fh;
      if(cdt.is_face(boost::get<0>(tf), boost::get<1>(tf), boost::get<2>(tf), fh)){
	*out = fh;
	out++;
      }
    } 
  }
};




public:
  typedef Tr                                   Triangulation;
  typedef typename Tr::Intersection_tag        Intersection_tag;
  typedef Polyline_constrained_triangulation_2<Tr> Self;
  typedef Tr                                   Base;

  typedef typename Triangulation::Edge             Edge;
  typedef typename Triangulation::Vertex           Vertex;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Face_handle      Face_handle;
  typedef typename Triangulation::Face_circulator  Face_circulator ;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Vertex_circulator  Vertex_circulator;
  typedef typename Triangulation::Locate_type      Locate_type;
  typedef typename Triangulation::Line_face_circulator Line_face_circulator;
  typedef typename Triangulation::Geom_traits      Geom_traits;
  typedef typename Geom_traits::Point_2            Point;
  typedef typename Geom_traits::Segment_2          Segment;
  typedef typename Triangulation::Constraint       Constraint;
  typedef typename Triangulation::size_type        size_type;

  typedef typename Triangulation::List_edges       List_edges;
  typedef typename Triangulation::List_faces       List_faces;
  typedef typename Triangulation::List_vertices    List_vertices;
  typedef typename Triangulation::List_constraints List_constraints;

  typedef Polyline_constraint_hierarchy_2<Vertex_handle, Point> Polyline_constraint_hierarchy;
  typedef Tag_true                                Polyline_constraint_hierarchy_tag;

  // for user interface with the constraint hierarchy
  typedef typename Polyline_constraint_hierarchy::H_vertex_it    
                                            Vertices_in_constraint_iterator;
  typedef typename Polyline_constraint_hierarchy::H_context      Context;
  typedef typename Polyline_constraint_hierarchy::H_context_iterator  Context_iterator;
  typedef typename Polyline_constraint_hierarchy::H_c_iterator   Constraint_iterator;
  typedef typename Polyline_constraint_hierarchy::H_sc_iterator  Subconstraint_iterator;
 
  typedef typename Polyline_constraint_hierarchy::H_vertex_list* Constraint_id;   
                                            
  //for backward compatibility
 typedef Vertices_in_constraint_iterator     Vertices_in_constraint;
  
  using Triangulation::geom_traits;
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::number_of_vertices;
  using Triangulation::vertices_begin;
  using Triangulation::locate;
  using Triangulation::includes_edge;
  using Triangulation::mark_constraint;
  using Triangulation::find_intersected_faces;
  using Triangulation::triangulate_hole;

protected:
  Polyline_constraint_hierarchy hierarchy;
 
public:
  Polyline_constraint_hierarchy& hierarchy_ref()
  {
    return hierarchy;
  }

  Polyline_constrained_triangulation_2(const Geom_traits& gt=Geom_traits()) 
    : Triangulation(gt) { }

  Polyline_constrained_triangulation_2(const Self& ctp)
    : Triangulation()    { copy_triangulation(ctp);}

  virtual ~Polyline_constrained_triangulation_2() {}

  Polyline_constrained_triangulation_2 & operator=(const Self& ctp);

  Polyline_constrained_triangulation_2(List_constraints& lc, 
				   const Geom_traits& gt=Geom_traits())
      : Triangulation(gt)
  {
    typename List_constraints::iterator lcit=lc.begin();
    for( ;lcit != lc.end(); lcit++) {
      insert_constraint( (*lcit).first, (*lcit).second);
    }
     CGAL_triangulation_postcondition( this->is_valid() );
  }

  template<class InputIterator>
  Polyline_constrained_triangulation_2(InputIterator first,
				   InputIterator last,
				   const Geom_traits& gt=Geom_traits() )
     : Triangulation(gt)
  {
    while( first != last){
      insert_constraint((*first).first, (*first).second);
      ++first;
    }
      CGAL_triangulation_postcondition( this->is_valid() );
  }

    //Helping
  void clear() { Base::clear(); hierarchy.clear();}
  void copy_triangulation(const Polyline_constrained_triangulation_2 &ctp);
  void swap(Polyline_constrained_triangulation_2 &ctp);

  // INSERTION
  Vertex_handle insert(const Point& a, 
		       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  Constraint_id insert_constraint(const Point& a, const Point& b);

  Vertices_in_constraint_iterator
  insert_vertex_in_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos, 
			      Vertex_handle vh)
  {
    return insert_vertex_in_constraint(cid, pos, vh, Emptyset_iterator());
  }

  Vertices_in_constraint_iterator
 remove_vertex_from_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos)
  {
    return remove_vertex_from_constraint(cid, pos, Emptyset_iterator());
  }


 template <class OutputIterator>
  Vertices_in_constraint_iterator
  remove_vertex_from_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos, 
				OutputIterator out)
  {
    if(pos == vertices_in_constraint_begin(cid)){
      ++pos;
      Constraint_id aux = hierarchy.split2(cid,pos);
      remove_constraint(aux, out);
      return pos;
    }

    Vertices_in_constraint_iterator it = vertices_in_constraint_end(cid);
    it--;
    if(pos == it){
      --pos;
      Constraint_id aux = hierarchy.split(cid, pos);
      remove_constraint(aux, out);
      return vertices_in_constraint_end(cid);
    }

    Vertices_in_constraint_iterator pp = pos;
    --pp;
    Vertex_handle a = *pp;
    pp = pos;
    ++pp;
    Vertex_handle b = *pp;
    --it;
    Vertices_in_constraint_iterator beg = vertices_in_constraint_begin(cid);
    ++beg;
    Face_container<Self> fc(*this);

    Constraint_id head = 0, tail = 0;
    if(pos != beg){
      // split off head
      --pos;
      head = hierarchy.split2(cid, pos);
      ++pos;
    }
    if(pos != it){
      // split off tail
      ++pos;
      tail = hierarchy.split(cid,pos);
    }
  
    Constraint_id aux = insert_constraint(a, b, std::back_inserter(fc));
    pos = vertices_in_constraint_end(aux);
    --pos;
    --pos; // this is not necessarily == vertices_in_constraint_begin(aux);
    hierarchy.swap(cid, aux);
    remove_constraint(aux, std::back_inserter(fc));

    if(head){
      hierarchy.concatenate2(head, cid);
    }

    if(tail){
      hierarchy.concatenate(cid, tail);
    }
    fc.write_faces(out);
    ++pos; // we went one too far back because the last vertex gets removed by concatenate
    return pos;
  }

  // Inserts vh before pos
  // Returns an iterator pointing on the newly inserted vertex
  // Writes the modified faces to out
  template <class OutputIterator>
  Vertices_in_constraint_iterator
  insert_vertex_in_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos, 
			      Vertex_handle vh, OutputIterator out)
  {
    // Insertion before the first vertex
    if(pos == vertices_in_constraint_begin(cid)){
      //std::cout << "insertion before first vertex" << std::endl;
      Constraint_id head = insert_constraint(vh, *pos, out);
      hierarchy.concatenate2(head, cid);
      return vertices_in_constraint_begin(cid);
    } 

    // Insertion after the last vertex
    if(pos == vertices_in_constraint_end(cid)){
      //std::cout << "insertion after last vertex" << std::endl;
      pos--;
      Constraint_id tail = insert_constraint(*pos, vh, out);
      pos = vertices_in_constraint_end(tail);
      --pos;
      hierarchy.concatenate(cid, tail);
      return pos;
    }
    Vertex_handle b = *pos;
    --pos;    
    Vertex_handle a = *pos;
    ++pos;
    Face_container<Self> fc(*this);
    Vertices_in_constraint_iterator beg = vertices_in_constraint_begin(cid), vcit;
    ++beg;
    vcit = beg;
    ++beg;
    // If the constraint consists only of a segment, and we want to insert
    // in the middle 
    if((pos == vcit) && (beg == vertices_in_constraint_end(cid))){
      //std::cout << "insertion in constraint which is a segment" << std::endl;
      Constraint_id aux1 = insert_constraint(a, vh, std::back_inserter(fc));
      Constraint_id aux2 = insert_constraint(vh, b, std::back_inserter(fc));
      pos = vertices_in_constraint_begin(aux2);
      concatenate(aux1, aux2);
      hierarchy.swap(cid, aux1);
      remove_constraint(aux1, std::back_inserter(fc));
      fc.write_faces(out);
      return pos;
      
    }
    Constraint_id head = 0, tail = 0;
    Vertices_in_constraint_iterator bit = vertices_in_constraint_begin(cid);
    Vertices_in_constraint_iterator pred = pos;
    --pred;
    ++bit;
    if(pos != bit){
      //std::cout << "split head" << std::endl;
      head = split(cid, pred);
      std::swap(head,cid); // split2 does the job
      pred = vertices_in_constraint_begin(cid);
      pos = pred;
      ++pos;
    }
    Vertices_in_constraint_iterator eit = vertices_in_constraint_end(cid);
    --eit;
    if(pos != eit){
      //std::cout << "split tail" << std::endl;
      tail = split(cid, pos);    
    }
    
    // make the new constraint
    Constraint_id aux1 = insert_constraint(a, vh, std::back_inserter(fc));
    Constraint_id aux2 = insert_constraint(vh, b, std::back_inserter(fc));
    pos = vertices_in_constraint_begin(aux2);
    concatenate(aux1, aux2);

    if(head){
      //std::cout << "concatenate head" << std::endl;
      remove_constraint(cid, std::back_inserter(fc));
      hierarchy.concatenate(head, aux1);
    } else {
      hierarchy.swap(cid, aux1);
      remove_constraint(aux1, std::back_inserter(fc));
      head = cid;
    }

    if(tail){
      //std::cout << "concatenate tail" << std::endl;
      concatenate(head, tail);
    }
    fc.write_faces(out);
    return pos;
  }

  template < class InputIterator, class OutputIterator>
  Constraint_id insert_polyline(InputIterator first, InputIterator last, OutputIterator out)
  {
    Face_handle hint;
    Face_container<Self> fc(*this);
    std::vector<Vertex_handle> vertices;
    for(;first!= last; first++){
      Vertex_handle vh = insert(*first, hint);
      hint = vh->face();
      // no duplicates
      if(vertices.empty() || (vertices.back() != vh)){
	vertices.push_back(vh);
      }
    }
    int n = vertices.size();
    if(n == 1){
      return NULL;
    }
    Constraint_id ca = hierarchy.insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1], std::back_inserter(fc)); 

    if(n>2){
      for(int j=1; j<n-1; j++){
	hierarchy.append_constraint(ca, vertices[j], vertices[j+1]);
	insert_subconstraint(vertices[j], vertices[j+1], std::back_inserter(fc));
      }
    }
    for(Vertices_in_constraint_iterator vcit = vertices_in_constraint_begin(ca);
	vcit != vertices_in_constraint_end(ca);
	vcit++){
      insert_incident_faces(vcit, out);
    }
    vertices_in_constraint_begin(ca)->fixed = true;
    Vertices_in_constraint_iterator end = vertices_in_constraint_end(ca);
    --end;
    end->fixed = true;
    fc.write_faces(out);
    
    return ca;
  }


  Constraint_id insert_constraint(Vertex_handle va, Vertex_handle vb)
  {
  // protects against inserting twice the same constraint
  Constraint_id cid = hierarchy.insert_constraint(va, vb);
  if (va != vb && (cid != NULL) )  insert_subconstraint(va,vb); 
  
  return cid;
  }

  template < class InputIterator>
  Constraint_id insert_polyline(InputIterator first, InputIterator last)
  {
  Face_handle hint;
    std::vector<Vertex_handle> vertices;
    for(;first!= last; first++){
      Vertex_handle vh = insert(*first, hint);
      hint = vh->face();
      // no duplicates
      if(vertices.empty() || (vertices.back() != vh)){
	vertices.push_back(vh);
      }
    }
    int n = vertices.size();
    if(n == 1){
      return NULL;
    }
    Constraint_id ca = hierarchy.insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1]); 

    if(n>2){
      for(int j=1; j<n-1; j++){
	hierarchy.append_constraint(ca, vertices[j], vertices[j+1]);
	insert_subconstraint(vertices[j], vertices[j+1]);
      }
    }
 
    vertices_in_constraint_begin(ca)->fixed = true;
    Vertices_in_constraint_iterator end = vertices_in_constraint_end(ca);
    --end;
    end->fixed = true;
    
    return ca;
  }


template <class OutputIterator>
typename Polyline_constrained_triangulation_2<Tr>::Constraint_id
insert_constraint(Vertex_handle va, Vertex_handle vb, OutputIterator out)
{
  // protects against inserting twice the same constraint
  Constraint_id cid = hierarchy.insert_constraint(va, vb);
  if (va != vb && (cid != NULL) )  insert_subconstraint(va,vb,out); 
  
    for(Vertices_in_constraint_iterator vcit = vertices_in_constraint_begin(cid);
	vcit != vertices_in_constraint_end(cid);
	vcit++){
      insert_incident_faces(vcit, out);
    }
  return cid;
}



  Vertex_handle push_back(const Point& a);
  Constraint_id push_back(const Constraint& c);
  
  //for backward compatibility
  // not const Point&, because otherwise VC6/7 messes it up with 
  // the insert that takes an iterator range
  Constraint_id insert(Point a, Point b) { insert_constraint(a, b);}
  Constraint_id insert(Vertex_handle va, Vertex_handle  vb) {insert_constraint(va,vb);}

  virtual Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  No_intersection_tag);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Exact_intersections_tag);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Exact_predicates_tag);
 
  // REMOVAL

  template <class OutputIterator>
  void remove_constraint(Constraint_id cid, OutputIterator out)
  {
    std::list<Vertex_handle> vertices(hierarchy.vertices_in_constraint_begin(cid),
				      hierarchy.vertices_in_constraint_end(cid));

    hierarchy.remove_constraint(cid);
    for(typename std::list<Vertex_handle>::iterator it = vertices.begin(), 
	  succ = it; 
	++succ != vertices.end(); 
	++it){
      if(! is_subconstraint(*it, *succ)){ // this checks whether other constraints pass
	Face_handle fh;
	int i;
	bool b = Triangulation::is_edge(*it, *succ, fh, i);
	CGAL_triangulation_assertion(b);
	Triangulation::remove_constrained_edge(fh,i, out); // this does also flipping if necessary.
      }
    }
  }
  void remove_constraint(Constraint_id cid)
  {
    remove_constraint(cid, Emptyset_iterator());
  }


  void simplify(Vertices_in_constraint_iterator u,
                Vertices_in_constraint_iterator v,
                Vertices_in_constraint_iterator w)

  {
    hierarchy.simplify(u,v,w);
    
    Triangulation::remove_incident_constraints(v->vertex);
  
    Triangulation::remove(v->vertex);
  
    Triangulation::insert_constraint(u->vertex,w->vertex);
  }

  // CONCATENATE AND SPLIT

  // concatenates two constraints
  Constraint_id
  concatenate(Constraint_id first, Constraint_id second);

  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint 
  Constraint_id
  split(Constraint_id first, Vertices_in_constraint_iterator vcit);
  
  // Query of the constraint hierarchy
  Constraint_iterator constraints_begin() const;
  Constraint_iterator constraints_end()   const;
  Subconstraint_iterator subconstraints_begin() const;
  Subconstraint_iterator subconstraints_end() const;
  Context   context(Vertex_handle va, Vertex_handle vb); 

  bool is_subconstraint(Vertex_handle va, 
			Vertex_handle vb);
  size_type number_of_enclosing_constraints(Vertex_handle va, 
				      Vertex_handle vb);
  Context_iterator   contexts_begin(Vertex_handle va, 
				    Vertex_handle vb);
  Context_iterator   contexts_end(Vertex_handle va, 
				  Vertex_handle vb);

  Vertices_in_constraint_iterator vertices_in_constraint_begin(Constraint_id cid) const;
  Vertices_in_constraint_iterator vertices_in_constraint_end(Constraint_id cid) const ;

  size_type number_of_constraints() {
    return static_cast<size_type> (hierarchy.number_of_constraints());}
  size_type number_of_subconstraints(){
    return static_cast<size_type> (hierarchy.number_of_subconstraints());}

  // public member, used by Mesh_2::Refine_edges
  void split_constraint(Vertex_handle v1, Vertex_handle v2,
                        Vertex_handle va) {
    hierarchy.split_constraint(v1,v2,va);
  }

protected:
  template <class OutputItertator>
  void insert_incident_faces(Vertices_in_constraint_iterator vcit, OutputItertator out)
  {
    Vertex_handle vh = vcit->vertex;
    Face_circulator fc = incident_faces(vh), done = fc;
    Face_circulator null ;
    if ( fc != null )
    {
      do {
        Face_handle fh = fc;
        out = fh;
        out++;
        fc++;
      }while(fc != done);
    }
  }  


void
insert_subconstraint(Vertex_handle vaa,
		     Vertex_handle vbb)
  {
    insert_subconstraint(vaa,vbb,Emptyset_iterator());
  }




template <class OutputItertator>
void
insert_subconstraint(Vertex_handle vaa,
		     Vertex_handle vbb,
		     OutputItertator out)
  // insert the subconstraint [vaa vbb] 
  // it will eventually be splitted into several subconstraints
{
  CGAL_triangulation_precondition( vaa != vbb);
  Vertex_handle vi;

  Face_handle fr;
  int i;
  if(includes_edge(vaa,vbb,vi,fr,i)) {
    mark_constraint(fr,i);
    if (vi != vbb)  {
      hierarchy.split_constraint(vaa,vbb,vi);
      insert_subconstraint(vi,vbb, out);
    }
    return;
  }
      
  List_faces intersected_faces;
  List_edges conflict_boundary_ab, conflict_boundary_ba;
     
  bool intersection  = find_intersected_faces( vaa, vbb,
			                       intersected_faces,
					       conflict_boundary_ab,
					       conflict_boundary_ba,
					       vi);
  if ( intersection) {
    if (vi != vaa && vi != vbb) {
      hierarchy.split_constraint(vaa,vbb,vi);
      insert_subconstraint(vaa,vi, out); 
      insert_subconstraint(vi,vbb, out); 
     }
    else insert_subconstraint(vaa,vbb,out);  

    
    return;
  }


  //no intersection

  List_edges edges(conflict_boundary_ab);
  std::copy(conflict_boundary_ba.begin(), conflict_boundary_ba.end(), std::back_inserter(edges));

  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba);

  get_bounded_faces(edges.begin(),
		    edges.end(),
		    out);

  if (vi != vbb) {
    hierarchy.split_constraint(vaa,vbb,vi);
    insert_subconstraint(vi,vbb, out); 
  }
  return;
}



  //to debug
public:
  void print_hierarchy() { hierarchy.print(); }

  //template member functions
public:
  template < class InputIterator >
#if defined(_MSC_VER)
  std::ptrdiff_t insert(InputIterator first, InputIterator last, int i = 0)
#else
    std::ptrdiff_t insert(InputIterator first, InputIterator last) 
#endif
  {
    size_type n = number_of_vertices();

    std::vector<Point> points (first, last);

    spatial_sort (points.begin(), points.end(), geom_traits());

    Face_handle hint;
    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
            p != end; ++p)
        hint = insert (*p, hint)->face();

    return number_of_vertices() - n;
  }

};

template <class Tr>
void
Polyline_constrained_triangulation_2<Tr>::
copy_triangulation(const Polyline_constrained_triangulation_2 &ctp)
{
  Base::copy_triangulation(ctp);
  //the following assume that the triangulation and its copy
  // iterates on their vertices on the same order 
  std::map<Vertex_handle,Vertex_handle> vmap;
  Vertex_iterator vit = ctp.vertices_begin();
  Vertex_iterator vvit = vertices_begin();
  for( ; vit != ctp.vertices_end(); ++vit, ++vvit) {
    CGAL_triangulation_assertion(vit->point() == vvit->point());
    vmap[vit] = vvit;
  }
  hierarchy.copy(ctp.hierarchy, vmap);
}

template <class Tr>
void
Polyline_constrained_triangulation_2<Tr>::
swap(Polyline_constrained_triangulation_2 &ctp)
{
  Base::swap(ctp);
  hierarchy.swap(ctp.hierarchy);
}

template <class Tr>
Polyline_constrained_triangulation_2<Tr>&
Polyline_constrained_triangulation_2<Tr>::
operator=(const Polyline_constrained_triangulation_2 &ctp)
{
  copy_triangulation(ctp);
  return *this;
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Vertex_handle
Polyline_constrained_triangulation_2<Tr>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Constraint_id
Polyline_constrained_triangulation_2<Tr>::
push_back(const Constraint &c)
{
  insert_constraint(c.first, c.second);
}
    
template < class Tr >
inline 
typename Polyline_constrained_triangulation_2<Tr>::Vertex_handle
Polyline_constrained_triangulation_2<Tr>::
insert(const Point& a, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(a, lt, li, start);
  return insert(a,lt,loc,li);
}

template < class Tr>
typename Polyline_constrained_triangulation_2<Tr>::Vertex_handle
Polyline_constrained_triangulation_2<Tr>::
insert(const Point& a, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v1, v2;
  bool insert_in_constrained_edge = false;

  if ( lt == Triangulation::EDGE && loc->is_constrained(li) ){
    insert_in_constrained_edge = true;
    v1=loc->vertex(ccw(li)); //endpoint of the constraint
    v2=loc->vertex(cw(li)); // endpoint of the constraint
  }
  Vertex_handle va = Triangulation::insert(a,lt,loc,li);
  // update the hierarchy
  if (insert_in_constrained_edge) {
    hierarchy.split_constraint(v1,v2,va);
  }
  return va;
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Constraint_id
Polyline_constrained_triangulation_2<Tr>::
insert_constraint(const Point& a, const Point& b)
  // insert endpoints first
{
  Vertex_handle va= insert(a);
  // If the segment is "short" it is a good idea to start the next insertion
  // close to point a
  // Otherwise, to start here is as good as elsewhere
  Vertex_handle vb = insert(b, va->face());
  return insert_constraint(va, vb); 
}


template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>:: Vertex_handle 
Polyline_constrained_triangulation_2<Tr>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb) 
{
  return intersect(f, i, vaa, vbb, Intersection_tag());
}

template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>:: Vertex_handle 
Polyline_constrained_triangulation_2<Tr>::

intersect(Face_handle , int , 
	  Vertex_handle ,
	  Vertex_handle ,
	  No_intersection_tag)
{
  std::cerr << " sorry, this triangulation does not deal with" 
	    <<    std::endl
	    << " intersecting constraints" << std::endl;
  CGAL_triangulation_assertion(false);
  return Vertex_handle();
}

template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>:: Vertex_handle 
Polyline_constrained_triangulation_2<Tr>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Exact_intersections_tag)
// compute the intersection of the constraint edge (f,i) 
// with the subconstraint (vaa,vbb) being inserted
// insert the intersection point
// (the  constraint edge (f,i) will be splitted in hierarchy by insert)
// and return the Vertex_handle of the new Vertex
{
  Vertex_handle  vc, vd, va, vb;
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));
  CGAL_triangulation_assertion_code( bool b1 = )
  hierarchy.enclosing_constraint(vcc,vdd,vc,vd);
  CGAL_triangulation_assertion_code( bool b2 = )
  hierarchy.enclosing_constraint(vaa,vbb,va,vb);
  CGAL_triangulation_assertion(b1);
  CGAL_triangulation_assertion(b2);

  const Point& pa = va->point();
  const Point& pb = vb->point();
  const Point& pc = vc->point();
  const Point& pd = vd->point();
  Point pi;
  Intersection_tag itag = Intersection_tag();
  CGAL_triangulation_assertion_code( bool ok = )
  intersection(geom_traits(), pa, pb, pc, pd, pi, itag );
  CGAL_triangulation_assertion(ok);

  Vertex_handle vi = insert(pi, Triangulation::EDGE, f, i);
  return vi; 
}

template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>::Vertex_handle 
Polyline_constrained_triangulation_2<Tr>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Exact_predicates_tag)
{
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));

  const Point& pa = vaa->point();
  const Point& pb = vbb->point();
  const Point& pc = vcc->point();
  const Point& pd = vdd->point();

  Point pi; //creator for point is required here
  Intersection_tag itag = Intersection_tag();
  bool ok  = intersection(geom_traits(), pa, pb, pc, pd, pi, itag );

  Vertex_handle vi;
  if ( !ok) {  //intersection detected but not computed
    int i = limit_intersection(geom_traits(), pa, pb, pc, pd, itag);
    switch(i){
    case 0 : vi = vaa; break;
    case 1 : vi = vbb; break;
    case 2 : vi = vcc; break;
    case 3 : vi = vdd; break; 
    }
    if(vi == vaa || vi == vbb) {
      Triangulation::remove_constrained_edge(f, i);
    }
  }
  else{ //computed
    Triangulation::remove_constrained_edge(f, i);
    vi = insert(pi, f);
  }

  // vi == vc or vi == vd may happen even if intersection==true
  // due to approximate construction of the intersection
  if (vi != vcc && vi != vdd) { 
    hierarchy.split_constraint(vcc,vdd,vi);
    insert_subconstraint(vcc,vi); 
    insert_subconstraint(vi, vdd);
  } 
  else {
    insert_subconstraint(vcc,vdd);
  }
  return vi; 
}

  // CONCATENATE AND SPLIT

  // concatenates two constraints
template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>::Constraint_id
Polyline_constrained_triangulation_2<Tr>::concatenate(Constraint_id first, Constraint_id second)
{
  return hierarchy.concatenate(first,second);
}

  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint 
template <class Tr>
typename Polyline_constrained_triangulation_2<Tr>::Constraint_id
Polyline_constrained_triangulation_2<Tr>::split(Constraint_id first, Vertices_in_constraint_iterator vcit)
{
  return hierarchy.split(first, vcit);
}


template <class Tr>
std::ostream &
operator<<(std::ostream& os, 
	   const Polyline_constrained_triangulation_2<Tr> &ct)
{
  ct.file_output(os);
  return os ;
}

// Constraint Hierarchy Queries

template <class Tr>
inline
typename
Polyline_constrained_triangulation_2<Tr>::Constraint_iterator
Polyline_constrained_triangulation_2<Tr>::
constraints_begin() const
{
  return hierarchy.c_begin();
}

template <class Tr>
inline
typename
Polyline_constrained_triangulation_2<Tr>::Constraint_iterator
Polyline_constrained_triangulation_2<Tr>::
constraints_end() const
{
  return hierarchy.c_end();
}

template <class Tr>
inline
typename
Polyline_constrained_triangulation_2<Tr>::Subconstraint_iterator
Polyline_constrained_triangulation_2<Tr>::
subconstraints_begin() const
{
  return hierarchy.sc_begin();
}

template <class Tr>
inline
typename
Polyline_constrained_triangulation_2<Tr>::Subconstraint_iterator
Polyline_constrained_triangulation_2<Tr>::
subconstraints_end() const
{
  return hierarchy.sc_end();
}


template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Context
Polyline_constrained_triangulation_2<Tr>::
context(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.context(va,vb);
}


template <class Tr>
inline 
typename Polyline_constrained_triangulation_2<Tr>::size_type
Polyline_constrained_triangulation_2<Tr>::
number_of_enclosing_constraints(Vertex_handle va, Vertex_handle vb)
{
 return static_cast<size_type> 
   (hierarchy.number_of_enclosing_constraints(va,vb)); 
}

template <class Tr>
inline bool
Polyline_constrained_triangulation_2<Tr>::
is_subconstraint(Vertex_handle va, Vertex_handle vb)
{
 return hierarchy.is_subconstrained_edge(va,vb); 
}


template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Context_iterator
Polyline_constrained_triangulation_2<Tr>::
contexts_begin(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.contexts_begin(va,vb);
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Context_iterator
Polyline_constrained_triangulation_2<Tr>::
contexts_end(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.contexts_end(va,vb);
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Vertices_in_constraint_iterator
Polyline_constrained_triangulation_2<Tr>::
vertices_in_constraint_begin(Constraint_id cid) const
{
  return  hierarchy.vertices_in_constraint_begin(cid);
}

template <class Tr>
inline
typename Polyline_constrained_triangulation_2<Tr>::Vertices_in_constraint_iterator
Polyline_constrained_triangulation_2<Tr>::
vertices_in_constraint_end(Constraint_id cid) const
{
  return  hierarchy.vertices_in_constraint_end(cid);
}

} //namespace CGAL

#endif //CGAL_POLYLINE_CONSTRAINED_TRIANGULATION_2_H
