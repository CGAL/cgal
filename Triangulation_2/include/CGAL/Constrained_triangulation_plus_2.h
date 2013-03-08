// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Constraint_hierarchy_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355)
//warning C4355: 'this' : used in base member initializer list
#endif

namespace CGAL {

// Comparison functor that compares two Vertex_handle.
// Used as 'Compare' functor for the constraint hierarchy.
template < class Tr >
class Ctp2_vertex_handle_less_xy {
  const Tr* tr_p;

public:
  Ctp2_vertex_handle_less_xy(const Tr* tr_p) : tr_p(tr_p) {}

  typedef typename Tr::Vertex_handle Vertex_handle;

  bool operator()(const Vertex_handle& va,
                  const Vertex_handle& vb) const
  {
    return tr_p->compare_xy(va->point(), vb->point()) == SMALLER;
  }
}; // end class template Ctp2_vertex_handle_less_xy


// Tr the base triangulation class 
// Tr has to be Constrained or Constrained_Delaunay

template < class Tr >
class Constrained_triangulation_plus_2  
  : public Tr
{
public:
  typedef Tr                                   Triangulation;
  typedef typename Tr::Intersection_tag        Intersection_tag;
  typedef Constrained_triangulation_plus_2<Tr> Self;
  typedef Tr                                   Base;

  typedef typename Triangulation::Edge             Edge;
  typedef typename Triangulation::Vertex           Vertex;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Face_handle      Face_handle;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
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

  typedef Ctp2_vertex_handle_less_xy<Self>         Vh_less_xy;
  typedef Constraint_hierarchy_2<Vertex_handle, 
                                 Vh_less_xy,
                                 bool>             Constraint_hierarchy;
  typedef Tag_true                                 Constraint_hierarchy_tag;

  // for user interface with the constraint hierarchy
  typedef typename Constraint_hierarchy::H_vertex_it    
                                            Vertices_in_constraint_iterator;
  typedef typename Constraint_hierarchy::H_context      Context;
  typedef typename Constraint_hierarchy::H_context_iterator  Context_iterator;
  typedef typename Constraint_hierarchy::H_c_iterator   Constraint_iterator;
  typedef typename Constraint_hierarchy::H_sc_iterator  Subconstraint_iterator;
 
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
  Constraint_hierarchy hierarchy;
 
public:
  Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits()) 
    : Triangulation(gt)
    , hierarchy(Vh_less_xy(this)) 
  { }

  Constrained_triangulation_plus_2(const Self& ctp)
    : Triangulation()
    , hierarchy(Vh_less_xy(this))
  { copy_triangulation(ctp);}

  virtual ~Constrained_triangulation_plus_2() {}

  Constrained_triangulation_plus_2 & operator=(const Self& ctp);

  Constrained_triangulation_plus_2(List_constraints& lc, 
				   const Geom_traits& gt=Geom_traits())
      : Triangulation(gt)
      , hierarchy(Vh_less_xy(this))
  {
    typename List_constraints::iterator lcit=lc.begin();
    for( ;lcit != lc.end(); lcit++) {
      insert_constraint( (*lcit).first, (*lcit).second);
    }
     CGAL_triangulation_postcondition( this->is_valid() );
  }

  template<class InputIterator>
  Constrained_triangulation_plus_2(InputIterator first,
				   InputIterator last,
				   const Geom_traits& gt=Geom_traits() )
     : Triangulation(gt)
    , hierarchy(Vh_less_xy(this))
  {
    while( first != last){
      insert_constraint((*first).first, (*first).second);
      ++first;
    }
      CGAL_triangulation_postcondition( this->is_valid() );
  }

    //Helping
  void clear() { Base::clear(); hierarchy.clear();}
  void copy_triangulation(const Constrained_triangulation_plus_2 &ctp);
  void swap(Constrained_triangulation_plus_2 &ctp);

  // INSERTION
  Vertex_handle insert(const Point& a, 
		       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  void insert_constraint(const Point& a, const Point& b);

  void insert_constraint(Vertex_handle va, Vertex_handle vb);
//   template < class InputIterator >
//   std::ptrdiff_t insert(InputIterator first, InputIterator last);
  Vertex_handle push_back(const Point& a);
  void          push_back(const Constraint& c);
  
  //for backward compatibility
  // not const Point&, because otherwise VC6/7 messes it up with 
  // the insert that takes an iterator range
  void insert(Point a, Point b) { insert_constraint(a, b);}
  void insert(Vertex_handle va, Vertex_handle  vb) {insert_constraint(va,vb);}

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
  void remove_constraint(Vertex_handle va, Vertex_handle vb);

  
  // Query of the constraint hierarchy
  Constraint_iterator constraints_begin() const;
  Constraint_iterator constraints_end()   const;
  Subconstraint_iterator subconstraints_begin() const;
  Subconstraint_iterator subconstraints_end() const;
  Context   context(Vertex_handle va, Vertex_handle vb); 
  size_type number_of_enclosing_constraints(Vertex_handle va, 
				      Vertex_handle vb);
  Context_iterator   contexts_begin(Vertex_handle va, 
				    Vertex_handle vb);
  Context_iterator   contexts_end(Vertex_handle va, 
				  Vertex_handle vb);
  Vertices_in_constraint_iterator 
                     vertices_in_constraint_begin(Vertex_handle va, 
						  Vertex_handle vb);
  Vertices_in_constraint_iterator 
                     vertices_in_constraint_end(Vertex_handle va, 
						Vertex_handle vb);
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
  void insert_subconstraint(Vertex_handle va,Vertex_handle vb);

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

   template <class OutputItFaces>  
     OutputItFaces  
     remove_constraint(Vertex_handle va,   
 		      Vertex_handle vb,   
 		      OutputItFaces out) {  
     CGAL_triangulation_precondition(hierarchy.is_constrained_edge(va,vb));  
     for(Vertices_in_constraint_iterator it =   
 	  hierarchy.vertices_in_constraint_begin(va,vb),  succ = it;   
 	++succ != hierarchy.vertices_in_constraint_end(va,vb);   
 	++it){  
       Face_handle fh;  
       int i;  
       CGAL_assume_code( bool b = )
 	Triangulation::is_edge(*it, *succ, fh, i);  
       CGAL_assume(b);
       // this does also flipping if necessary.  
       Triangulation::remove_constrained_edge(fh,i,out);  
 
     }   
     hierarchy.remove_constraint(va,vb);  
     return out;  
   }  

};

template <class Tr>
void
Constrained_triangulation_plus_2<Tr>::
copy_triangulation(const Constrained_triangulation_plus_2 &ctp)
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
Constrained_triangulation_plus_2<Tr>::
swap(Constrained_triangulation_plus_2 &ctp)
{
  Base::swap(ctp);
  hierarchy.swap(ctp.hierarchy);
}

template <class Tr>
Constrained_triangulation_plus_2<Tr>&
Constrained_triangulation_plus_2<Tr>::
operator=(const Constrained_triangulation_plus_2 &ctp)
{
  copy_triangulation(ctp);
  return *this;
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Tr>
inline
void
Constrained_triangulation_plus_2<Tr>::
push_back(const Constraint &c)
{
  insert_constraint(c.first, c.second);
}
    
template < class Tr >
inline 
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
insert(const Point& a, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(a, lt, li, start);
  return insert(a,lt,loc,li);
}

template < class Tr>
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
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
inline void
Constrained_triangulation_plus_2<Tr>::
insert_constraint(const Point& a, const Point& b)
  // insert endpoints first
{
  Vertex_handle va= insert(a);
  // If the segment is "short" it is a good idea to start the next insertion
  // close to point a
  // Otherwise, to start here is as good as elsewhere
  Vertex_handle vb = insert(b, va->face());
  insert_constraint(va, vb); 
}

template <class Tr>
inline void
Constrained_triangulation_plus_2<Tr>::
insert_constraint(Vertex_handle va, Vertex_handle vb)
{
  if( va != vb ){
    // protects against inserting twice the same constraint
    bool no_twice = hierarchy.insert_constraint(va, vb);
    if ( no_twice ){
      insert_subconstraint(va,vb);
    } 
  }
    return;
}


template <class Tr>
inline void
Constrained_triangulation_plus_2<Tr>::
insert_subconstraint(Vertex_handle vaa,
		     Vertex_handle vbb)
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
      insert_subconstraint(vi,vbb);
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
      insert_subconstraint(vaa,vi); 
      insert_subconstraint(vi,vbb); 
     }
    else insert_subconstraint(vaa,vbb);
    return;
  }

  //no intersection
  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba);

  if (vi != vbb) {
    hierarchy.split_constraint(vaa,vbb,vi);
    insert_subconstraint(vi,vbb); 
  }
  return;
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb) 
{
  return intersect(f, i, vaa, vbb, Intersection_tag());
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr>::

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
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr>::
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
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle 
Constrained_triangulation_plus_2<Tr>::
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


// This remove_constraint leaves all vertices in the triangulation
template <class Tr>
void
Constrained_triangulation_plus_2<Tr>::
remove_constraint(Vertex_handle va, Vertex_handle vb)
{
    remove_constraint(va, vb, Emptyset_iterator());
}

template <class Tr>
std::ostream &
operator<<(std::ostream& os, 
	   const Constrained_triangulation_plus_2<Tr> &ct)
{
  ct.file_output(os);
  return os ;
}

// Constraint Hierarchy Queries

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Constraint_iterator
Constrained_triangulation_plus_2<Tr>::
constraints_begin() const
{
  return hierarchy.c_begin();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Constraint_iterator
Constrained_triangulation_plus_2<Tr>::
constraints_end() const
{
  return hierarchy.c_end();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr>::
subconstraints_begin() const
{
  return hierarchy.sc_begin();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr>::
subconstraints_end() const
{
  return hierarchy.sc_end();
}


template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context
Constrained_triangulation_plus_2<Tr>::
context(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.context(va,vb);
}


template <class Tr>
inline 
typename Constrained_triangulation_plus_2<Tr>::size_type
Constrained_triangulation_plus_2<Tr>::
number_of_enclosing_constraints(Vertex_handle va, Vertex_handle vb)
{
 return static_cast<size_type> 
   (hierarchy.number_of_enclosing_constraints(va,vb)); 
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context_iterator
Constrained_triangulation_plus_2<Tr>::
contexts_begin(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.contexts_begin(va,vb);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context_iterator
Constrained_triangulation_plus_2<Tr>::
contexts_end(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.contexts_end(va,vb);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Vertices_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
vertices_in_constraint_begin(Vertex_handle va, Vertex_handle vb)
{
  return  hierarchy.vertices_in_constraint_begin(va,vb);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Vertices_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
vertices_in_constraint_end(Vertex_handle va, Vertex_handle vb)
{
  return  hierarchy.vertices_in_constraint_end(va,vb);
}

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif //CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
