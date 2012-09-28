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
// Author(s)     : Mariette Yvinec, Jean-Daniel Boissonnat


#ifndef CGAL_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_2_H

#include <set>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_2.h> 
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/iterator.h>

#include <CGAL/intersections.h>
#include <CGAL/squared_distance_2.h>

namespace CGAL {

struct No_intersection_tag{};
struct Exact_intersections_tag{}; // to be used with an exact number type
struct Exact_predicates_tag{}; // to be used with filtered exact number


template < class Gt, 
           class Tds = Triangulation_data_structure_2 <
                       Triangulation_vertex_base_2<Gt>,
		       Constrained_triangulation_face_base_2<Gt> >, 
           class Itag = No_intersection_tag >
class Constrained_triangulation_2  : public Triangulation_2<Gt,Tds>
{

public:
  typedef Triangulation_2<Gt,Tds> Triangulation;
  typedef Constrained_triangulation_2<Gt,Tds,Itag>  Constrained_triangulation;
  
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::size_type size_type;
  typedef typename Triangulation::Locate_type Locate_type;
  typedef typename Triangulation::All_faces_iterator All_faces_iterator;
  typedef typename Triangulation::Face_circulator Face_circulator;
  typedef typename Triangulation::Edge_circulator Edge_circulator;
  typedef typename Triangulation::Vertex_circulator Vertex_circulator;
  typedef typename Triangulation::Line_face_circulator Line_face_circulator;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::number_of_vertices;
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::dimension;
  using Triangulation::geom_traits;
  using Triangulation::all_faces_begin;
  using Triangulation::all_faces_end;
  using Triangulation::side_of_oriented_circle;
  using Triangulation::is_infinite;
  using Triangulation::collinear_between;
  using Triangulation::incident_edges;
  using Triangulation::test_dim_down;
  using Triangulation::make_hole;
  using Triangulation::fill_hole;
  using Triangulation::delete_vertex;
  using Triangulation::delete_face;
  using Triangulation::create_face;
  using Triangulation::incident_faces;
  using Triangulation::locate;
  using Triangulation::includes_edge;
  using Triangulation::remove_first;
  using Triangulation::remove_second;
#endif

  typedef Gt                                 Geom_traits;
  typedef Itag                               Intersection_tag;
  typedef typename Geom_traits::Point_2      Point;
  typedef typename Geom_traits::Segment_2    Segment;
  typedef std::pair<Point,Point>             Constraint;
  typedef std::list<Edge>                    List_edges;
  typedef std::list<Face_handle>             List_faces;
  typedef std::list<Vertex_handle>           List_vertices;
  typedef std::list<Constraint>              List_constraints;

  // Tag to mark the presence of a hierarchy of constraints
 typedef Tag_false                           Constraint_hierarchy_tag;
   

  class Less_edge;
  typedef std::set<Edge,Less_edge> Edge_set;


  Constrained_triangulation_2(const Gt& gt = Gt()) : Triangulation(gt) { }

  Constrained_triangulation_2(const Constrained_triangulation_2& ct)
    : Triangulation(ct) {}

  Constrained_triangulation_2(std::list<Constraint>& lc, const Gt& gt=Gt())
      : Triangulation_2<Gt,Tds>(gt)
  {
    typename List_constraints::iterator lcit=lc.begin();
    for( ;lcit != lc.end(); lcit++) {
      insert( (*lcit).first, (*lcit).second);
    }
     CGAL_triangulation_postcondition( this->is_valid() );
  }

  template<class InputIterator>
  Constrained_triangulation_2(InputIterator it,
			      InputIterator last,
			      const Gt& gt=Gt() )
     : Triangulation_2<Gt,Tds>(gt)
  {
    for ( ; it != last; it++) {
      	insert_constraint((*it).first, (*it).second);
      }
      CGAL_triangulation_postcondition( this->is_valid() );
  }

  //TODO Is that destructor correct ?
  virtual ~Constrained_triangulation_2() {}

  // INSERTION
  Vertex_handle insert(const Point& p, 
			       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, 
		       int li );
  Vertex_handle push_back(const Point& a);
//   template < class InputIterator >
//   std::ptrdiff_t insert(InputIterator first, InputIterator last);
 
  void insert_constraint(const Point& a, const Point& b);
  void insert_constraint(Vertex_handle va, Vertex_handle  vb);
  void push_back(const Constraint& c);

  void remove(Vertex_handle  v);
  void remove_constrained_edge(Face_handle f, int i);
  void remove_incident_constraints(Vertex_handle  v);
  // to be used by Constrained_triangulation_plus_2
 template <class OutputItFaces>
 OutputItFaces 
 remove_constrained_edge(Face_handle f, int i, OutputItFaces out) 
 {
   remove_constrained_edge(f, i);
   return out;
 }

  //for backward compatibility
  void remove_constraint(Face_handle f, int i) {remove_constrained_edge(f,i);}
  void insert(Point a, Point b) { insert_constraint(a, b);}
  void insert(Vertex_handle va, Vertex_handle  vb) {insert_constraint(va,vb);}

  // QUERY
  bool is_constrained(Edge e) const;
  bool are_there_incident_constraints(Vertex_handle v) const;
  bool is_valid(bool verbose = false, int level = 0) const;
  // template<class OutputItEdges>
  // OutputItEdges incident_constraints(Vertex_handle v, 
  //                                     OutputItEdges out) const;
  

  class Less_edge 
    :  public std::binary_function<Edge, Edge, bool>
  {
  public:
    Less_edge() {}
    bool operator() (const Edge& e1, const Edge& e2) const
      {
	int ind1=e1.second, ind2=e2.second;
 	return( (&(*e1.first) < &(*e2.first))
 		|| ( (&(*e1.first) == &(*e2.first)) && (ind1 < ind2)));
      } 
  };


  void file_output(std::ostream& os) const;

protected:
  virtual Vertex_handle virtual_insert(const Point& a, 
				       Face_handle start = Face_handle());
  virtual Vertex_handle virtual_insert(const Point& a,
				       Locate_type lt,
				       Face_handle loc, 
				       int li );
//Vertex_handle special_insert_in_edge(const Point & a, Face_handle f, int i);
  void update_constraints_incident(Vertex_handle va, 
				   Vertex_handle c1,
				   Vertex_handle c2);
  void clear_constraints_incident(Vertex_handle va);
  void update_constraints_opposite(Vertex_handle va);
  void update_constraints(const List_edges &hole);

  void mark_constraint(Face_handle fr, int i);

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
public:
// made public for Laurent  to find out deleted faces
// when inserting a constraint with most probably 
// no intersection
  bool find_intersected_faces(Vertex_handle va, 
			      Vertex_handle vb,
			      List_faces & intersected_faces,
			      List_edges & list_ab, 
			      List_edges & list_ba,
			      Vertex_handle& vi);
protected:
  virtual void triangulate_hole(List_faces& intersected_faces,
				List_edges& conflict_boundary_ab,
				List_edges& conflict_boundary_ba);
  
  void triangulate_hole(List_faces& intersected_faces,
			List_edges& conflict_boundary_ab,
			List_edges& conflict_boundary_ba,
			List_edges& new_edges);
  
  void triangulate_half_hole(List_edges & list_edges, 
			     List_edges & new_edges);

  void remove_1D(Vertex_handle v);
  void remove_2D(Vertex_handle v);

  //templated member function
public:
  // the int parameter is a work around for VC7 to compile
  template < class InputIterator >
#if defined(_MSC_VER)
  std::ptrdiff_t insert(InputIterator first, InputIterator last, int i = 0)
#else
    std::ptrdiff_t insert(InputIterator first, InputIterator last) 
#endif
    {
      size_type n = number_of_vertices(); 

      std::vector<Point> points (first, last);
      CGAL::spatial_sort (points.begin(), points.end(), geom_traits());

      Face_handle hint;
      for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
              p != end; ++p)
          hint = insert (*p, hint)->face();

      return number_of_vertices() - n;
    }

  //deprecated
 template<class OutputIterator>
  bool are_there_incident_constraints(Vertex_handle v, 
				      OutputIterator out) const
    {
      Edge_circulator ec=incident_edges(v), done(ec);
      bool are_there = false;
      if (ec == 0) return are_there;
      do {
	if(is_constrained(*ec)) {
	  *out++ = *ec;
	  are_there = true;
	}
	ec++;
      } while (ec != done);
      return are_there;
    }

  
 template<class OutputItEdges>
 OutputItEdges  incident_constraints(Vertex_handle v, 
				      OutputItEdges out) const {
   Edge_circulator ec=incident_edges(v), done(ec);
   if (ec == 0) return  out;
   do {
     if(is_constrained(*ec))    *out++ = *ec;
     ec++;
   } while (ec != done);
   return out;
 }
 
  // the following fonctions are overloaded 
  // to take care of constraint marks 
  template<class EdgeIt>
  Vertex_handle star_hole( const Point& p, 
			   EdgeIt edge_begin,  
			   EdgeIt edge_end) {
    std::list<Face_handle> empty_list;
    return star_hole(p, 
		     edge_begin, 
		     edge_end, 
		     empty_list.begin(),
		     empty_list.end());
  }

  template<class EdgeIt, class FaceIt>
  Vertex_handle star_hole( const Point& p, 
			   EdgeIt edge_begin,
			   EdgeIt edge_end,
			   FaceIt face_begin,
			   FaceIt face_end)
{
    Vertex_handle v =  Triangulation::star_hole(p, 
						edge_begin, 
						edge_end, 
						face_begin,
						face_end);
    // restore constraint status for new faces.
    int vindex;
    Face_handle fh;
    int ih;
    Face_circulator fc = incident_faces(v), done(fc);
    do {
      vindex = fc->index(v);
      fc->set_constraint(cw(vindex), false);
      fc->set_constraint(ccw(vindex), false);
      fh = fc->neighbor(vindex);
      ih = this->mirror_index(fc,vindex);
      fc->set_constraint(vindex, fh->is_constrained(ih));
    } while (++fc != done);
    return v;
}
 
};

template < class Gt, class Tds, class Itag >
inline
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_triangulation_2<Gt,Tds,Itag>::
virtual_insert(const Point& a, Face_handle start)
// virtual version of insert
{
  return insert(a,start);
}

template < class Gt, class Tds, class Itag >
inline
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_triangulation_2<Gt,Tds,Itag>::
virtual_insert(const Point& a,
	       Locate_type lt,
	       Face_handle loc, 
	       int li )
// virtual version of insert
{
  return insert(a,lt,loc,li);
}
    
template < class Gt, class Tds, class Itag >
inline
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_triangulation_2<Gt,Tds,Itag>::
insert(const Point& a, Face_handle start)
// inserts point a 
// in addition to what is done for non constrained triangulations
// constrained edges are updated
{
  Face_handle loc;
  int li;
  Locate_type lt;
  loc = locate(a, lt, li, start);
  return Constrained_triangulation::insert(a,lt,loc,li);
}

template < class Gt, class Tds, class Itag >
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_triangulation_2<Gt,Tds,Itag>::
insert(const Point& a, Locate_type lt, Face_handle loc, int li)
// insert a point p, whose localisation is known (lt, f, i)
// in addition to what is done for non constrained triangulations
// constrained edges are updated
{
  Vertex_handle va;
  Vertex_handle v1, v2;
  bool insert_in_constrained_edge = false;

  if ( lt == Triangulation::EDGE && loc->is_constrained(li) ){
    insert_in_constrained_edge = true;
    v1=loc->vertex(ccw(li)); //endpoint of the constraint
    v2=loc->vertex(cw(li)); // endpoint of the constraint
  }

  va = Triangulation::insert(a,lt,loc,li);
  if (insert_in_constrained_edge) update_constraints_incident(va, v1,v2);
  else if(lt != Triangulation::VERTEX) clear_constraints_incident(va);
  if (dimension() == 2) update_constraints_opposite(va);
  return va;
}

// template < class Gt, class Tds, class Itag >  
// typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
// Constrained_triangulation_2<Gt, Tds, Itag>::
// special_insert_in_edge(const Point & a, Face_handle f, int i)
//   // insert  point p in edge(f,i)
//   // bypass the precondition for point a to be in edge(f,i)
//   // update constrained status
// {
//   Vertex_handle va;
//   Vertex_handle c1,c2;
//   c1 = f->vertex(cw(i));  //endpoint of edge
//   c2 = f->vertex(ccw(i)); //endpoint of edge
//   bool insert_in_constrained_edge = f->is_constrained(i);
 
//   va = this->_tds.insert_in_edge(f, i);
//   va->set_point(a);

//   if (insert_in_constrained_edge) update_constraints_incident(va, c1,c2);
//   else clear_constraints_incident(va);
//   if (dimension() == 2) update_constraints_opposite(va);
//   return va;
// }


template < class Gt, class Tds, class Itag >
inline void
Constrained_triangulation_2<Gt,Tds,Itag>::
insert_constraint(const Point& a, const Point& b)
// the algorithm first inserts a and b, 
// and then forces the constraint [va,vb]
{
  Vertex_handle va= virtual_insert(a);
  Vertex_handle vb= virtual_insert(b);
  if ( va != vb)   insert_constraint(va,vb);
}



template <class Gt, class Tds, class Itag >
inline void
Constrained_triangulation_2<Gt,Tds,Itag>::
insert_constraint(Vertex_handle  vaa, Vertex_handle vbb)
// forces the constrained [va,vb]
// [va,vb] will eventually be splitted into several edges
// if a vertex vc of t lies on segment ab
// of if ab intersect some constrained edges
{
  CGAL_triangulation_precondition( vaa != vbb);
  Vertex_handle vi;

  Face_handle fr;
  int i;
  if(includes_edge(vaa,vbb,vi,fr,i)) {
    mark_constraint(fr,i);
    if (vi != vbb)  {
      insert_constraint(vi,vbb);
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
      insert_constraint(vaa,vi); 
      insert_constraint(vi,vbb); 
     }
    else insert_constraint(vaa,vbb);
    return;
  }

  //no intersection
  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba);

  if (vi != vbb) {
    insert_constraint(vi,vbb); 
  }
  return;

}

template <class Gt, class Tds, class Itag >
bool
Constrained_triangulation_2<Gt,Tds,Itag>::
find_intersected_faces(Vertex_handle vaa,
		       Vertex_handle vbb,
		       List_faces & intersected_faces,
		       List_edges & list_ab, 
		       List_edges & list_ba,
		       Vertex_handle & vi) 
  // vi is set to the first vertex of the triangulation on [vaa,vbb].
  // Return true if an intersection with a constrained edge is
  // encountered, false otherwise
  // When false : 
  // intersected_faces contains the list if faces intersected by [va,vi]
  // list_ab and list_ba represents the boundary of the union
  // of the intersected faces oriented cw
  // list_ab consists of the edges from vaa to vi (i.e. on the left of a->b)
  // list_ba    "         "        from vi to vaa (i.e. on the right of a->b)
{
  const Point& aa = vaa->point();
  const Point& bb = vbb->point();
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
    orient = this->orientation(aa,bb,current_vertex->point());
    int i1, i2;
    switch (orient) {
    case COLLINEAR :  
      done = true; // current_vertex is the new endpoint
      break;
    case LEFT_TURN :
    case RIGHT_TURN :
      if (orient == LEFT_TURN) {
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
	if (orient == LEFT_TURN) 
	  list_ab.push_back(Edge(lf, lf->index(current_face)));
	else // orient == RIGHT_TURN
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


template <class Gt, class Tds, class Itag >
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_triangulation_2<Gt,Tds,Itag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb) 
{
  return intersect(f, i, vaa, vbb, Itag());
}

template <class Gt, class Tds, class Itag >
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_triangulation_2<Gt,Tds,Itag>::
intersect(Face_handle , int , 
	  Vertex_handle ,
	  Vertex_handle ,
	  No_intersection_tag)
{
  //SL: I added that to be able to throw while we find a better solution
  #ifdef CGAL_CT2_WANTS_TO_HAVE_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS
  CGAL_CDT2_EXTRA_ACTION_FOR_INTERSECTING_CONSTRAINTS
  #endif
  
  std::cerr << " sorry, this triangulation does not deal with" 
	    <<    std::endl
	    << " intersecting constraints" << std::endl;
  CGAL_triangulation_assertion(false);
  return Vertex_handle() ;
}

template <class Gt, class Tds, class Itag >
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_triangulation_2<Gt,Tds,Itag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Exact_intersections_tag)
// compute the intersection of the constraint edge (f,i) 
// with the subconstraint (vaa,vbb) being inserted
// insert the intersection point
// split constraint edge (f,i) 
// and return the Vertex_handle of the new Vertex
{ 
  std::cerr << "You are using an exact number types" << std::endl;
  std::cerr << "using a Constrained_triangulation_plus_2 class" << std::endl;
  std::cerr << "would avoid cascading intersection computation" << std::endl;
  std::cerr << " and be much more efficient" << std::endl;
  const Point& pa = vaa->point();
  const Point& pb = vbb->point();
  const Point& pc = f->vertex(cw(i))->point();
  const Point& pd = f->vertex(ccw(i))->point();
  Point pi;
  Itag itag = Itag();
  CGAL_triangulation_assertion_code( bool ok = )
  intersection(geom_traits(), pa, pb, pc, pd, pi, itag );
  CGAL_triangulation_assertion(ok);
  Vertex_handle vi = virtual_insert(pi, Triangulation::EDGE, f, i);
  return vi; 
}

template <class Gt, class Tds, class Itag >
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle 
Constrained_triangulation_2<Gt,Tds,Itag>::
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
  Itag itag = Itag();
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
      remove_constrained_edge(f, i);
    }
  }
  else{ //intersection computed
    remove_constrained_edge(f, i);
    vi = virtual_insert(pi, f);
  }

  // vi == vc or vi == vd may happen even if intersection==true
  // due to approximate construction of the intersection
  if (vi != vcc && vi != vdd) { 
    insert_constraint(vcc,vi); 
    insert_constraint(vi, vdd);
  } 
  else {
    insert_constraint(vcc,vdd);
  }
  return vi; 
}

template <class Gt, class Tds, class Itag >
inline
typename Constrained_triangulation_2<Gt,Tds,Itag>::Vertex_handle
Constrained_triangulation_2<Gt,Tds,Itag>::
push_back(const Point &p)
{
  return insert(p);
}


template <class Gt, class Tds, class Itag >
inline
void
Constrained_triangulation_2<Gt,Tds,Itag>::
push_back(const Constraint &c)
{
  insert(c.first, c.second);
}


template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
update_constraints_incident(Vertex_handle va, 
			    Vertex_handle c1,
			    Vertex_handle c2)
  // update status of edges incident to a 
  // after insertion in the  constrained edge c1c2
{
  if (dimension() == 0) return;
  if (dimension()== 1) {
    Edge_circulator ec=this->incident_edges(va), done(ec);
    do {
      ((*ec).first)->set_constraint(2,true);
    }while (++ec != done);
  }
  else{
    //dimension() ==2
    int cwi, ccwi, indf;
    Face_circulator fc=this->incident_faces(va), done(fc);  
    CGAL_triangulation_assertion(fc != 0);
    do {
      indf = fc->index(va);
      cwi=cw(indf);
      ccwi=ccw(indf); 
      if ((fc->vertex(cwi) == c1)||(fc->vertex(cwi) == c2)) {
	  fc->set_constraint(ccwi,true);
	  fc->set_constraint(cwi,false);
	}	
	else {
	  fc->set_constraint(ccwi,false);
	  fc->set_constraint(cwi,true);
	}
	++fc;
      } while (fc != done);
  }
}

template < class Gt, class Tds ,class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
clear_constraints_incident(Vertex_handle va)
// make the edges incident to a newly created vertex unconstrained
{
 Edge_circulator ec=this->incident_edges(va), done(ec);
 Face_handle f;
 int indf;
  if ( ec != 0){
    do {
      f = (*ec).first ;
      indf = (*ec).second;
      f->set_constraint(indf,false);
      if (dimension() == 2) {
	f->neighbor(indf)->set_constraint(this->mirror_index(f,indf),false);
      }
    } while (++ec != done);
  }
  return;
}


template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::  
update_constraints_opposite(Vertex_handle va)
  // update status of edges opposite to a
  // after insertion of a
{
  CGAL_triangulation_assertion(dimension()==2); 
  Face_handle f=va->face(), start=f;
  int indf;
  do {
    indf = f->index(va);
    if (f->neighbor(indf)->is_constrained(this->mirror_index(f,indf)) ) {
      f->set_constraint(indf,true);
    }
    else {
      f->set_constraint(indf,false);
    }
    f= f->neighbor(ccw(indf)); // turns ccw around va 
  } while (f != start);
  return;
}

template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>:: 
update_constraints( const List_edges &hole)
{
  typename List_edges::const_iterator it = hole.begin();
  Face_handle f;
  int i;
  for ( ; it != hole.end(); it ++) {
    f =(*it).first;
    i = (*it).second;
    if ( f->is_constrained(i) ) 
      (f->neighbor(i))->set_constraint(this->mirror_index(f,i),true);
    else (f->neighbor(i))->set_constraint(this->mirror_index(f,i),false);
  }
}


template < class Gt, class Tds, class Itag >
inline void
Constrained_triangulation_2<Gt,Tds,Itag>::
mark_constraint(Face_handle fr, int i)
{
  if (dimension()==1) fr->set_constraint(2, true);
  else{
    fr->set_constraint(i,true);
    fr->neighbor(i)->set_constraint(this->mirror_index(fr,i),true);
  }
  return;
}

template < class Gt, class Tds, class Itag >
inline void
Constrained_triangulation_2<Gt,Tds,Itag>::
triangulate_hole(List_faces& intersected_faces,
		 List_edges& conflict_boundary_ab,
		 List_edges& conflict_boundary_ba)
{
  List_edges new_edges;
  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba,
		   new_edges);
}



template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
triangulate_hole(List_faces& intersected_faces,
		 List_edges& conflict_boundary_ab,
		 List_edges& conflict_boundary_ba,
		 List_edges& new_edges)
  // triangulate the hole limited by conflict_boundary_ab
  // and conflict_boundary_ba
  // insert the new edges in new-edges 
  // delete the faces of intersected_faces
{
  if ( !conflict_boundary_ab.empty() ) {
    triangulate_half_hole(conflict_boundary_ab, new_edges);
    triangulate_half_hole(conflict_boundary_ba, new_edges);
	
    // the two faces that share edge ab are neighbors
    // their common edge ab is a constraint
    Face_handle fr,fl;
    fl=(*conflict_boundary_ab.begin()).first;
    fr=(*conflict_boundary_ba.begin()).first;
    fl->set_neighbor(2, fr);
    fr->set_neighbor(2, fl);
    fl->set_constraint(2, true);
    fr->set_constraint(2, true);
   
    // delete intersected faces
    while( ! intersected_faces.empty()) {
      fl = intersected_faces.front();
      intersected_faces.pop_front();
      delete_face(fl);
    }
  }
}



template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
remove(Vertex_handle  v)
  // remove a vertex and updates the constrained edges of the new faces
  // precondition : there is no incident constraints
{
  CGAL_triangulation_precondition( v != Vertex_handle() );
  CGAL_triangulation_precondition( ! is_infinite(v));
  CGAL_triangulation_precondition( ! are_there_incident_constraints(v));
    
  if  (number_of_vertices() == 1)     remove_first(v);
  else if (number_of_vertices() == 2) remove_second(v);
  else   if ( dimension() == 1) remove_1D(v);
  else  remove_2D(v);
  return;
}


template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
remove_1D(Vertex_handle  v)
{
  Edge_circulator ec = incident_edges(v), done(ec);
  do {
    (*ec).first->set_constraint(2,false);
  } while (++ec != done);
  Triangulation::remove_1D(v);
}

template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) { this->_tds.remove_dim_down(v);}
  else {
    List_edges hole;
    make_hole(v, hole);
    List_edges shell=hole; //save hole because it will be emptied by fill_hole
    fill_hole(v, hole);
    update_constraints(shell);
    delete_vertex(v);
  }
  return;       
}


template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
remove_constrained_edge(Face_handle f, int i)
{
  f->set_constraint(i, false);
  if (dimension() == 2)
    (f->neighbor(i))->set_constraint(this->mirror_index(f,i), false);
  return;
}

template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
remove_incident_constraints(Vertex_handle v)
{
   Edge_circulator ec=incident_edges(v), done(ec);
   if (ec == 0) return;
   do {
	if(is_constrained(*ec)) { remove_constrained_edge((*ec).first,
						   (*ec).second);}
	ec++;
   } while (ec != done);
   return;	
}

template < class Gt, class Tds, class Itag >
inline  bool 
Constrained_triangulation_2<Gt,Tds,Itag>::
are_there_incident_constraints(Vertex_handle v) const
{
  return are_there_incident_constraints(v, Emptyset_iterator());
}

template < class Gt, class Tds, class Itag >
inline  bool 
Constrained_triangulation_2<Gt,Tds,Itag>::
is_valid(bool verbose, int level) const
{
    bool result = Triangulation::is_valid(verbose,level);
    for( All_faces_iterator it = all_faces_begin(); 
	                    it != all_faces_end() ; it++) {
      for(int i=0; i<3; i++) {
	Face_handle n = it->neighbor(i);
	result = result && 
	  it->is_constrained(i) == n->is_constrained(n->index(it));
      }
    }
    return result;
}

template < class Gt, class Tds, class Itag >
inline  bool 
Constrained_triangulation_2<Gt,Tds,Itag>::
is_constrained(Edge e) const
{
  return (e.first)->is_constrained(e.second);
}
    
template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt,Tds,Itag>::
triangulate_half_hole(List_edges & list_edges,  List_edges & new_edges)
  // triangulates the  polygon whose boundary consists of list_edges
  // plus the edge ab joining the two endpoints of list_edges
  // the orientation of the polygon (as provided by list_edges) must
  // be cw
  // the edges of list_edges are assumed to be edges of a
  // triangulation that will be updated by the procedure
  // the edges that are created are put in list new_edges
  // takes linear time
{
  Vertex_handle va; // first vertex of list_edges 
  Face_handle newlf;
  Face_handle n1,n2,n;
  int ind1, ind2,ind;
  Orientation orient;
    
  typename List_edges::iterator current, next, tempo;
  current=list_edges.begin();

  va=((*current).first)->vertex(ccw((*current).second));
  next=current; 
  ++next;

  do
    {
      n1=(*current).first;
      ind1=(*current).second;
      // in case n1 is no longer a triangle of the new triangulation
      if ( n1->neighbor(ind1) != Face_handle() ) {
	n=n1->neighbor(ind1);
	//ind=this->mirror_index(n1,ind1); 
	// mirror_index does not work in this case
	ind = cw(n->index(n1->vertex(cw(ind1))));
	n1=n->neighbor(ind); 
	ind1= this->mirror_index(n,ind);
      }
      n2=(*next).first;
      ind2=(*next).second;
      // in case n2 is no longer a triangle of the new triangulation
      if (n2->neighbor(ind2) != Face_handle() ) {
	n=n2->neighbor(ind2); 
	// ind=this->mirror_index(n2,ind2);
	// mirror_index does not work in this case
	ind = cw(n->index(n2->vertex(cw(ind2))));
	n2=n->neighbor(ind); 
	ind2= this->mirror_index(n,ind);
      }

      Vertex_handle v0=n1->vertex(ccw(ind1));
      Vertex_handle v1=n1->vertex(cw(ind1));
      Vertex_handle v2=n2->vertex(cw(ind2));
      orient = this->orientation(v0->point(),v1->point(),v2->point());
      switch (orient) {
      case RIGHT_TURN : 	  		
	// creates the new triangle v0v1v2
	// updates the neighbors, the constraints 
	//and the list of new edges
	newlf = create_face(v0,v2,v1);
	new_edges.push_back(Edge(newlf,2));
	newlf->set_neighbor(1, n1);
	newlf->set_neighbor(0, n2);
	n1->set_neighbor(ind1, newlf);
	n2->set_neighbor(ind2, newlf);
	if (n1->is_constrained(ind1)) {
	  newlf->set_constraint(1,true);
	}
	if (n2->is_constrained(ind2)) {
	  newlf->set_constraint(0,true);
	}
	// v0, v1 or v2.face() may have been removed
	v0->set_face(newlf); 
	v1->set_face(newlf);
	v2->set_face(newlf);
	// update list_edges
	tempo=current;
	current=list_edges.insert(current, Edge(newlf,2));
	list_edges.erase(tempo);
	list_edges.erase(next);
	next=current;
	if (v0 != va) {--current;} 
	else {++next;} 
	break;
      case LEFT_TURN : 	  
	++current; ++next;
	break;
      case COLLINEAR : 
	++current; ++next;
	break;
      }
    } while (next != list_edges.end());
}

template < class Gt, class Tds, class Itag >
void
Constrained_triangulation_2<Gt, Tds, Itag>::
file_output(std::ostream& os) const
{
  Triangulation_2<Gt, Tds>::file_output(os);

  // write constrained status
  typename Tds::Face_iterator ib = this->_tds.face_iterator_base_begin();
  for( ; ib != this->_tds.face_iterator_base_end(); ++ib) {
    for(int j = 0; j < 3; ++j){
      if (ib->is_constrained(j)) { os << "C";}
      else { os << "N";}
      if(is_ascii(os)){
	if(j==2) {
	  os << "\n";
	} else {
	  os <<  ' ';
	}
      }
    }
  }
}

template < class Gt, class Tds, class Itag >
std::ostream &
operator<<(std::ostream& os, 
	   const Constrained_triangulation_2<Gt,Tds,Itag> &ct)
{
  ct.file_output(os);
  return os ;
}

template < class Gt, class Tds, class Itag >
std::istream &
operator>>(std::istream& is, 
	         Constrained_triangulation_2<Gt,Tds,Itag> &ct)
{
  typedef Constrained_triangulation_2<Gt,Tds,Itag> CDT;
  ct.clear();
  is >> static_cast<typename CDT::Triangulation&>(ct);
  for (typename CDT::All_faces_iterator fit=ct.all_faces_begin(),
                                        fit_end=ct.all_faces_end();fit_end!=fit;++fit){
    char c[3];
    is >> c[0] >>  c[1] >> c[2];
    for (int k=0;k<3;++k){
      fit->set_constraint(k,c[k]=='C');
    }
  }  
  return is;
}

//Helping functions to compute intersections of constrained edges
template<class Gt>
bool
intersection(const Gt& ,
	     const typename Gt::Point_2& , 
	     const typename Gt::Point_2& , 
	     const typename Gt::Point_2& , 
	     const typename Gt::Point_2& ,
	     typename Gt::Point_2& ,
	     No_intersection_tag)
{
  return false;
}
	     
template<class Gt>
bool
intersection(const Gt& gt,
	     const typename Gt::Point_2& pa, 
	     const typename Gt::Point_2& pb, 
	     const typename Gt::Point_2& pc, 
	     const typename Gt::Point_2& pd,
	     typename Gt::Point_2& pi,
	     Exact_intersections_tag)
{
  return compute_intersection(gt,pa,pb,pc,pd,pi);
}


template<class Gt>
inline bool
intersection(const Gt& gt,
	     const typename Gt::Point_2& pa, 
	     const typename Gt::Point_2& pb, 
	     const typename Gt::Point_2& pc, 
	     const typename Gt::Point_2& pd,
	     typename Gt::Point_2& pi,
	     Exact_predicates_tag)
{
  return compute_intersection(gt,pa,pb,pc,pd,pi);
}


template<class Gt>
bool
compute_intersection(const Gt& gt,
	     const typename Gt::Point_2& pa, 
	     const typename Gt::Point_2& pb, 
	     const typename Gt::Point_2& pc, 
	     const typename Gt::Point_2& pd,
	     typename Gt::Point_2& pi)
{
  typename Gt::Intersect_2 compute_intersec=gt.intersect_2_object();
   typename Gt::Construct_segment_2  
    construct_segment=gt.construct_segment_2_object();
  Object result = compute_intersec(construct_segment(pa,pb),
				   construct_segment(pc,pd));
  return assign(pi, result);
}


template<class Gt>
int
limit_intersection(const Gt& ,
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& ,
		   No_intersection_tag)
{
  return 0;
}

template<class Gt>
int
limit_intersection(const Gt& ,
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& , 
		   const typename Gt::Point_2& ,
		   Exact_intersections_tag)
{
  return 0;
}

template<class Gt>
int
limit_intersection(const Gt& gt,
	     const typename Gt::Point_2& pa, 
	     const typename Gt::Point_2& pb, 
	     const typename Gt::Point_2& pc, 
	     const typename Gt::Point_2& pd,
	     Exact_predicates_tag)
{
  typename Gt::Construct_line_2 line = gt.construct_line_2_object();
  typename Gt::Compute_squared_distance_2 
    distance = gt.compute_squared_distance_2_object();
  typename Gt::Line_2 l1 = line(pa,pb);
  typename Gt::Line_2 l2 = line(pc,pd);
  int i = 0;
  typename Gt::FT dx = distance(l2,pa);
  typename Gt::FT db = distance(l2,pb);
  typename Gt::FT dc = distance(l1,pc);
  typename Gt::FT dd = distance(l1,pd);
  if ( db < dx  ) { dx = db; i = 1;}
  if ( dc < dx  ) { dx = dc; i = 2;}
  if ( dd < dx  ) { i = 3;}
  return i;
}

} //namespace CGAL

#endif //CGAL_CONSTRAINED_TRIANGULATION_2_H
