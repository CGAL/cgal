// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Constrained_triangulation_plus_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$

// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Constraint_hierarchy_2.h>
#include <CGAL/squared_distance_2.h>

CGAL_BEGIN_NAMESPACE

struct Tag_no_intersection{};
struct Tag_exact_intersections{}; // to be used with an exact number type
struct Tag_exact_predicates{}; // to be used with filtered exact number

// Tr the base triangulation class 
// Tr has to be Constrained or Constrained_Delaunay
// I_tag to tag the optional support of constriaint intersections
template < class Tr, class I_tag = Tag_exact_predicates>
class Constrained_triangulation_plus_2  
  : public Tr
{
public:
  typedef Tr              Triangulation;
  typedef I_tag           Intersection_tag;

  typedef typename Triangulation::Edge         Edge;
  typedef typename Triangulation::Vertex        Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle   Face_handle;
  typedef typename Triangulation::Locate_type  Locate_type;
  typedef typename Triangulation::Line_face_circulator Line_face_circulator;
  typedef typename Triangulation::Geom_traits      Geom_traits;
  typedef typename Geom_traits::Point_2            Point;
  typedef typename Geom_traits::Segment_2          Segment;
  typedef typename Triangulation::Constraint       Constraint;

  typedef typename Triangulation::List_edges       List_edges;
  typedef typename Triangulation::List_faces       List_faces;
  typedef typename Triangulation::List_vertices    List_vertices;
  typedef typename Triangulation::List_constraints List_constraints;

  typedef Constraint_hierarchy_2<Vertex_handle, bool> Constraint_hierarchy;

  // for user interface with the constraint hierarchy
  typedef typename Constraint_hierarchy::H_vertex_it    Vertices_in_constraint;
  typedef typename Constraint_hierarchy::H_context      Context;
  typedef typename Constraint_hierarchy::H_context_iterator  Context_iterator;
  typedef typename Constraint_hierarchy::H_c_iterator   Constraint_iterator;
  typedef typename Constraint_hierarchy::H_sc_iterator  Subconstraint_iterator;
 

protected:
  Constraint_hierarchy hierarchy;
  Intersection_tag     itag;

public:
  Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits()) 
    : Triangulation(gt), itag() { }

  // copy constructrue et effectation operateur a revoir
//   Constrained_triangulation_plus_2(
  //  const Constrained_triangulation_plus_2& ct)
//     : Constrained_triangulation(ct) {}
  // destructeur aussi a revoir

  Constrained_triangulation_plus_2(List_constraints& lc, 
				   const Geom_traits& gt=Geom_traits())
      : Triangulation(gt), itag()
  {
    typename List_constraints::iterator lcit=lc.begin();
    for( ;lcit != lc.end(); lcit++) {
      insert( (*lcit).first, (*lcit).second);
    }
     CGAL_triangulation_postcondition( is_valid() );
  }

  template<class InputIterator>
  Constrained_triangulation_plus_2(InputIterator first,
				   InputIterator last,
				   const Geom_traits& gt=Geom_traits() )
     : Triangulation(gt)
  {
    while( first != last){
      insert((*first).first, (*first).second);
      ++first;
    }
      CGAL_triangulation_postcondition( is_valid() );
  }

  // INSERTION
  Vertex_handle insert(const Point& a, Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  void insert(Point a, Point b);
  void insert(Vertex_handle va, Vertex_handle vb);
//   template < class InputIterator >
//   int insert(InputIterator first, InputIterator last);
  Vertex_handle push_back(const Point& a);
  void          push_back(const Constraint& c);
  
 

  //SUPPRESSION
  // to be done next

  // Query of the constraint hierarchy
  Constraint_iterator constraints_begin() const;
  Constraint_iterator constraints_end()   const;
  Subconstraint_iterator subconstraints_begin() const;
  Subconstraint_iterator subconstraints_end() const;
  Context   context(Vertex_handle va, Vertex_handle vb); 
  int number_of_enclosing_constraints(Vertex_handle va, 
				      Vertex_handle vb);
  Context_iterator   contexts_begin(Vertex_handle va, 
				    Vertex_handle vb);
  Context_iterator   contexts_end(Vertex_handle va, 
				  Vertex_handle vb);
  Vertices_in_constraint vertices_in_constraint_begin(Vertex_handle va, 
						      Vertex_handle vb);
  Vertices_in_constraint vertices_in_constraint_end(Vertex_handle va, 
						    Vertex_handle vb);

protected:
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Tag_no_intersection);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Tag_exact_intersections);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Tag_exact_predicates);
  

  void  insert_subconstraint(Vertex_handle vaa,
			     Vertex_handle vbb);

  bool find_intersected_faces(Vertex_handle vaa,
			      Vertex_handle vbb,
			      List_faces & intersected_faces,
			      List_edges & list_ab, 
			      List_edges & list_ba,
			      Vertex_handle & vi);
private:
  Vertex_handle t_intersect(Vertex_handle vaa,
			    Vertex_handle vbb,
			    Vertex_handle vcc,
			    Vertex_handle vdd);

  //to debug
public:
  void print_hierarchy() { return hierarchy.print();}

  //template member functions
public:
  template < class InputIterator >
  int insert(InputIterator first, InputIterator last)
    {
      int n = number_of_vertices();
      while(first != last){
	insert(*first);
	++first;
      }
      return number_of_vertices() - n;
    }

};

template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
Constrained_triangulation_plus_2<Tr,I_tag>::
push_back(const Point &p)
{
  return insert(p);
}


template <class Tr, class I_tag >
inline
void
Constrained_triangulation_plus_2<Tr,I_tag>::
push_back(const Constraint &c)
{
  insert(c.first, c.second);
}
    
template < class Tr, class I_tag >
inline 
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
Constrained_triangulation_plus_2<Tr,I_tag>::
insert(const Point& a, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(a, lt, li, start);
  return insert(a,lt,loc,li);
}

template < class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
Constrained_triangulation_plus_2<Tr,I_tag>::
insert(const Point& a, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v1, v2;
  bool insert_in_constrained_edge = false;

  if ( lt == EDGE && loc->is_constrained(li) ){
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


template <class Tr, class I_tag >
inline void
Constrained_triangulation_plus_2<Tr,I_tag>::
insert(Point a, Point b)
  // insert endpoints first
{
  Vertex_handle va= insert(a);
  Vertex_handle vb= insert(b);
  insert(va, vb);
}

template <class Tr, class I_tag >
inline void
Constrained_triangulation_plus_2<Tr,I_tag>::
insert(Vertex_handle va, Vertex_handle vb)
{
  // protects against inserting twice the same constraint
  bool no_twice = hierarchy.insert_constraint(va, vb);
  if (va != vb && no_twice )  insert_subconstraint(va,vb); 
  return;
}


template <class Tr, class I_tag >
inline void
Constrained_triangulation_plus_2<Tr,I_tag>::
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
  List_edges new_edges;
    
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

  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba);

  if (vi != vbb) {
    hierarchy.split_constraint(vaa,vbb,vi);
    insert_subconstraint(vi,vbb); 
  }
  return;
}


    
template <class Tr, class I_tag >
bool
Constrained_triangulation_plus_2<Tr,I_tag>::
find_intersected_faces(Vertex_handle vaa,
		       Vertex_handle vbb,
		       List_faces & intersected_faces,
		       List_edges & list_ab, 
		       List_edges & list_ba,
		       Vertex_handle & vi)
  // find the facets intersected by  the subconstraint [vaa,vbb] 
  // and put in vertex-handle vi the first vertex encountered on [vaa,vbb]
  // return true if  vi is  an  intersection point 
  // when false, vi is not an intersection and : 
  // intersected_faces contains the list if faces intersected by [va,vi]
  // list_ab and list_ba represents the boundary of the union
  // of the intersected faces oriented cw
  // list_ab consists of the edges from vaa to vi (i.e. on the left of a->b)
  // list_ba    "         "        from vi to vaa (i.e. on the right of a->b)
{
  Point aa = vaa->point();
  Point bb = vbb->point();
  Line_face_circulator current_face=Line_face_circulator(vaa, this, bb);
  int ind=current_face->index(vaa);
      
  // to deal with the case where the first crossed edge
  // is constrained
  if(current_face->is_constrained(ind)) {
    vi=intersect(current_face, ind, vaa, vbb);
    return true;
  }

  Face_handle lf= current_face->neighbor(ccw(ind)); 
  Face_handle rf= current_face->neighbor(cw(ind));
  Orientation orient;
  Face_handle previous_face;
  Vertex_handle current_vertex;	

  list_ab.push_back(Edge(lf, lf->index(current_face)));
  list_ba.push_front(Edge(rf, rf->index(current_face)));
  intersected_faces.push_front(current_face);

  // initcd
  previous_face=current_face; 
  ++current_face;
  ind=current_face->index(previous_face);  
  current_vertex=current_face->vertex(ind);  

  // loop over triangles intersected by ab
  bool done = false;
  while (current_vertex != vbb && !done)  { 
    orient = orientation(aa,bb,current_vertex->point());
    int i1, i2;
    switch (orient) {
    case COLLINEAR :  
      done = true; // current_vertex is the new endpoint
      break;
    case LEFTTURN :
    case RIGHTTURN :
      if (orient == LEFTTURN) {
	i1 = ccw(ind) ; //index of second intersected edge of current_face
	i2 = cw(ind); //index of non intersected edge of current_face
      }
      else {
	i1 = cw(ind) ; //index of second intersected edge of current_face
	i2 = ccw(ind); //index of non intersected edge of current_face
      }
      if(current_face->is_constrained(i1)) {
	vi = intersect(current_face, i1, vaa,vbb);
	return true;
      }
      else {
	lf= current_face->neighbor(i2);
	intersected_faces.push_front(current_face);
	if (orient == LEFTTURN) 
	  list_ab.push_back(Edge(lf, lf->index(current_face)));
	else // orient == RIGHTTURN
	  list_ba.push_front(Edge(lf, lf->index(current_face)));
	previous_face=current_face;
	++current_face;
	ind=current_face->index(previous_face); 
	current_vertex=current_face->vertex(ind);
      }
      break;
    }
  }
    
  // last triangle 
  vi = current_vertex;
  intersected_faces.push_front(current_face);
  lf= current_face->neighbor(cw(ind));
  list_ab.push_back(Edge(lf, lf->index(current_face))); 
  rf= current_face->neighbor(ccw(ind));
  list_ba.push_front(Edge(rf, rf->index(current_face)));
  return false;
}

template <class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle 

Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb) 
{
  return intersect(f, i, vaa, vbb, Intersection_tag());
}

template <class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Tag_no_intersection)
{
  std::cerr << " sorry, this triangulation does not deal " 
	    <<    std::endl
	    << " intersecting constraints" << std::endl;
  CGAL_triangulation_assertion(false);
  return Vertex_handle();
}

template <class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Tag_exact_intersections)
// compute the intersection of the constraint edge (f,i) 
// with the subconstraint (vaa,vbb) being inserted
// insert the intersection point
// split constraint edge (f,i) in hierarchy
// and return the Vertex_handle of the new Vertex
{
  Vertex_handle  vc, vd, va, vb;
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));
  CGAL_triangulation_assertion(hierarchy.enclosing_constraint(vcc,vdd,vc,vd));
  CGAL_triangulation_assertion(hierarchy.enclosing_constraint(vaa,vbb,va,vb));
						  
  Point pi; //creator for point is required here
  Object result;
  typename Geom_traits::Intersect_2 
    compute_intersection=geom_traits().intersect_2_object();
  result = compute_intersection(Segment(vc->point(),vd->point()),
				Segment(va->point(),vb->point()));
  CGAL_triangulation_assertion(assign(pi, result));

  Vertex_handle vi = insert(pi, EDGE, f, i);
  return vi; 
}


template <class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Tag_exact_predicates)
{
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));

  Point pi; //creator for point is required here
  Object result;
  typename Geom_traits::Intersect_2 
    compute_intersection = geom_traits().intersect_2_object();
  result = compute_intersection(Segment(vcc->point(),vdd->point()),
				Segment(vaa->point(),vbb->point()));
  bool intersection = assign(pi, result);
  Vertex_handle vi;
  if ( !intersection) {  //intersection detected but not computed
    vi = t_intersect(vaa,vbb,vcc,vdd);
  }
  else{ //intersection detected but not computed
    remove_constraint(f, i);
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

template <class Tr, class I_tag >
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
t_intersect(Vertex_handle vaa,
	    Vertex_handle vbb,
	    Vertex_handle vcc,
	    Vertex_handle vdd)
{
  // intersection between edges [vaa,vbb] and [vcc,vdd]
  // has been detected by exact predicates
  // but not computed by approximate construction
  typename Geom_traits::Construct_line_2  construct_line = 
    geom_traits().construct_line_2_object();
  typename Geom_traits::Compute_squared_distance_2
    compute_squared_distance = 
    geom_traits().compute_squared_distance_2_object();
  typename Geom_traits::Line_2 l1 = construct_line(vaa->point(),
						   vbb->point());
  typename Geom_traits::Line_2 l2 = construct_line(vcc->point(),
						   vdd->point());
  Vertex_handle vi = vaa;
  typename Geom_traits::FT dd = compute_squared_distance(l2,vaa->point());
  if (compute_squared_distance(l2,vbb->point()) < dd) vi = vbb;
  if (compute_squared_distance(l1,vcc->point()) < dd) vi = vcc;
  if (compute_squared_distance(l1,vdd->point()) < dd) vi = vdd;
  return vi;
}


template <class Tr, class I_tag >
std::ostream &
operator<<(std::ostream& os, 
	   const Constrained_triangulation_plus_2<Tr,I_tag> &ct)
{
  ct.file_output(os);
  return os ;
}

// Constraint Hierarchy Queries

template <class Tr, class I_tag >
inline
typename
Constrained_triangulation_plus_2<Tr,I_tag>::Constraint_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
constraints_begin() const
{
  return hierarchy.c_begin();
}

template <class Tr, class I_tag >
inline
typename
Constrained_triangulation_plus_2<Tr,I_tag>::Constraint_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
constraints_end() const
{
  return hierarchy.c_end();
}

template <class Tr, class I_tag >
inline
typename
Constrained_triangulation_plus_2<Tr,I_tag>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
subconstraints_begin() const
{
  return hierarchy.sc_begin();
}

template <class Tr, class I_tag >
inline
typename
Constrained_triangulation_plus_2<Tr,I_tag>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
subconstraints_end() const
{
  return hierarchy.sc_end();
}


template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Context
Constrained_triangulation_plus_2<Tr,I_tag>::
context(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.context(va,vb);
}


template <class Tr, class I_tag >
inline int
Constrained_triangulation_plus_2<Tr,I_tag>::
number_of_enclosing_constraints(Vertex_handle va, Vertex_handle vb)
{
 return hierarchy.number_of_enclosing_constraints(va,vb); 
}

template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Context_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
contexts_begin(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.contexts_begin(va,vb);
}

template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Context_iterator
Constrained_triangulation_plus_2<Tr,I_tag>::
contexts_end(Vertex_handle va, Vertex_handle vb)
{
  return hierarchy.context_end(va,vb);
}

template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertices_in_constraint
Constrained_triangulation_plus_2<Tr,I_tag>::
vertices_in_constraint_begin(Vertex_handle va, Vertex_handle vb)
{
  return  hierarchy.vertices_in_constraint_begin(va,vb);
}

template <class Tr, class I_tag >
inline
typename Constrained_triangulation_plus_2<Tr,I_tag>::Vertices_in_constraint
Constrained_triangulation_plus_2<Tr,I_tag>::
vertices_in_constraint_end(Vertex_handle va, Vertex_handle vb)
{
  return  hierarchy.vertices_in_constraint_end(va,vb);
}

CGAL_END_NAMESPACE

#endif //CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H



