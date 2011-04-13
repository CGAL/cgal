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

CGAL_BEGIN_NAMESPACE

// Tr the base triangulation class 
// Tr has to be Constrained or Constrained_Delaunay
// I_tag to tag the optional support of constriaint intersections
template < class Tr, class I_tag = Tag_false>
class Constrained_triangulation_plus_2  
  : public Tr
{
public:
  typedef Tr              Triangulation;
  typedef I_tag           Intersection_tag;

  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Locate_type Locate_type;
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

protected:
  Constraint_hierarchy hierarchy;

public:
  Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits()) 
    : Triangulation(gt) { }

  // copy constructrue et effectation operateur a revoir
//   Constrained_triangulation_plus_2(
  //  const Constrained_triangulation_plus_2& ct)
//     : Constrained_triangulation(ct) {}

  Constrained_triangulation_plus_2(List_constraints& lc, 
				   const Geom_traits& gt=Geom_traits())
      : Triangulation(gt)
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
  Vertex_handle insert(Point a);
  Vertex_handle insert(const Point& p,
		       Locate_type lt,
		       Face_handle loc, int li );
  void insert(Point a, Point b);
  void insert(Vertex_handle va, Vertex_handle vb);
  Vertex_handle push_back(const Point& a);
  void          push_back(const Constraint& c);
  
 

  //SUPPRESSION
  // to be done next

protected:
  void update_new_edges(Point a, Point b, 
			Vertex_handle vi,
			Vertex_handle vh,
			List_edges& new_edges);

  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Tag_false);
  Vertex_handle intersect(Face_handle f, int i, 
			  Vertex_handle vaa,
			  Vertex_handle vbb,
			  Tag_true);

  Vertex_handle insert_part(Vertex_handle va, 
			    Vertex_handle vb,
			    Vertex_handle vaa);

  Vertex_handle find_intersected_faces(Vertex_handle va, 
				       Vertex_handle vb,
				       Vertex_handle vaa,
				       List_faces & intersected_faces,
				       List_edges & list_ab, 
				       List_edges & list_ba,
				       List_edges & new_edges);
				       

  //to debug
public:
  void print_hierarchy() { return hierarchy.print();}

};

template <class Tr, class I_tag >
inline
Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
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
Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
Constrained_triangulation_plus_2<Tr,I_tag>::
insert(Point a)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate(a, lt, li);
  return insert(a,lt,loc,li);
}

template < class Gt, class Tds >
Constrained_triangulation_plus_2<Gt,Tds>::Vertex_handle
Constrained_triangulation_plus_2<Gt,Tds>::
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
  hierarchy.insert_constraint(va, vb);
  Vertex_handle vaa = va;
  while (vaa != vb) {
    Vertex_handle vbb = insert_part(va,vb,vaa);
    if (vbb != vb) hierarchy.split_constraint(vaa,vb,vbb);
    vaa=vbb;
  }    
  return;
}


template <class Tr, class I_tag >
inline 
Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle
Constrained_triangulation_plus_2<Tr,I_tag>::
insert_part(Vertex_handle va, 
	    Vertex_handle vb,
	    Vertex_handle vaa)
{
  Vertex_handle vbb;

  Face_handle fr;
  int i;
  if(includes_edge(vaa,vb,vbb,fr,i)) {
    mark_constraint(fr,i);
    return vbb;
  }
      
  List_faces intersected_faces;
  List_edges conflict_boundary_ab, conflict_boundary_ba;
  List_edges new_edges;
    
  vbb = find_intersected_faces( va, vb, vaa,
			       intersected_faces,
			       conflict_boundary_ab,
			       conflict_boundary_ba,
			       new_edges);
			       
  triangulate_hole(intersected_faces,
		   conflict_boundary_ab,
		   conflict_boundary_ba,
		   new_edges,
		   vbb);
  return vbb;
}


    
template <class Tr, class I_tag >
Constrained_triangulation_plus_2<Tr,I_tag>::Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
find_intersected_faces(Vertex_handle va, 
		       Vertex_handle  vb,
		       Vertex_handle vaa,
		       List_faces & intersected_faces,
		       List_edges & list_ab, 
		       List_edges & list_ba,
		       List_edges & new_edges)
  // finds all triangles intersected the current part of constraint va 
  // vb
  // vaa is the vertex at the begin of the current part
  // the procedure returns the vertex vbb at the end of the current
  // part
  // The current part ends either when the current segment 
  // encounters a vertx of the triangulation or intersect
  // a constrained edge or the vertex vb  is reached.
  // At the end of the procedure
  //  List_faces return the list of intersectef faces by the current
  //  part
  // list_ab and list_ba represents the boundary of the union
  // of the intersected faces oriented cw
  // list_ab consists of the edges from vaa to vbb (i.e. on the left of a->b)
  // list_ba    "         "        from vbb to vaa (i.e. on the right of a->b)
  // an element of the lists (an edge e) is represented as the edge of
  // the triangle incident to e that is not intersected by ab
  // new edges contains the edges created by the intersection 
  // which are not constrained (and also some of the constrained ones)
{
  Point a=va->point(), b=vb->point();
  Line_face_circulator current_face=line_walk(vaa->point(),b, vaa->face());
  int ind=current_face->index(vaa);
  Vertex_handle vh, vhh; 
  Face_handle fh;
  int ih;
    
  if(current_face->is_constrained(ind)) {
    // to deal with the case where the first crossed edge
    // is constrained
    vh = current_face->mirror_vertex(ind);
    Vertex_handle vi=intersect(current_face, ind, vaa, vb);

    // to set the edge vaa vi as constrained
    assert(is_edge(vi,vaa,fh,ih));
    fh->set_constraint(ih,true);
    (fh->neighbor(ih))->set_constraint(fh->mirror_index(ih),true);
    // new_vertices.push_back(vi);
    //no need to insert (vi,vaa) in new_edges because it is constrained....
    update_new_edges(a,b,vi,vh,new_edges);
    return vi;
  }

  Vertex_handle vbb=vb;
  Face_handle lf= current_face->neighbor(ccw(ind)); 
  Face_handle rf= current_face->neighbor(cw(ind));
  Orientation orient;
  Face_handle previous_face;
  Vertex_handle current_vertex;	

  list_ab.push_back(Edge(lf, lf->index(current_face)));
  list_ba.push_front(Edge(rf, rf->index(current_face)));
  intersected_faces.push_front(current_face);

  // init
  previous_face=current_face; 
  ++current_face;
  ind=current_face->index(previous_face);  
  current_vertex=current_face->vertex(ind);  

  // loop over triangles intersected by ab
  while (current_vertex != vbb)  { 
    orient = orientation(a,b,current_vertex->point());
    int i1, i2;
    switch (orient) {
    case COLLINEAR :  
      vbb=current_vertex; // new endpoint of the constraint
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
	vhh = current_face->vertex(i1);
	vh = current_face->mirror_vertex(i1);
	Vertex_handle vi=intersect(current_face, i1, vaa,vb);
	//new_vertices.push_back(vi);
	update_new_edges(a,b,vi,vh,new_edges);
	update_new_edges(a,b,vi,vhh,new_edges);

	current_face=line_walk(vi->point(),a,vi->face());
 	// a essayer
 	// current_face=Line_face_circulator(vi,this,a);
	current_vertex = vi;
	ind = current_face->index(current_vertex);
	vbb = vi; // new endpoint of the constraint
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
    
  // last triangle (having vbb as a vertex)
  intersected_faces.push_front(current_face);
  lf= current_face->neighbor(cw(ind));
  list_ab.push_back(Edge(lf, lf->index(current_face))); 
  rf= current_face->neighbor(ccw(ind));
  list_ba.push_front(Edge(rf, rf->index(current_face)));
  return vbb;
}


template <class Tr, class I_tag >
void
Constrained_triangulation_plus_2<Tr,I_tag>:: 
update_new_edges(Point a, Point b, 
		 Vertex_handle vi,
		 Vertex_handle vh,
		 List_edges& new_edges)
  // add edge vi vh = (fh,ih) in the list of new_edges
  // taking  care that the face fh is not intersected by ab
  // and thus will not be destroyed
{
  Face_handle fh;
  int ih;
  assert( is_edge(vi,vh,fh,ih));
  if ( orientation(a,b,vh->point()) ==
       orientation(a,b,fh->vertex(ih)->point())) {
    // face fh is not traversed by ab and will not be deleted
    new_edges.push_back(Edge(fh,ih));
  }
  else {
    new_edges.push_back(Edge(fh->neighbor(ih),fh->mirror_index(ih)));
  }
}

template <class Tr, class I_tag >
Constrained_triangulation_plus_2<Tr,I_tag>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb) 
{
  return intersect(f, i, vaa, vbb, Intersection_tag());
}

template <class Tr, class I_tag >
Constrained_triangulation_plus_2<Tr,I_tag>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Tag_false)
{
  std::cerr << " sorry, this triangulation does not deal " 
	    <<    std::endl
	    << " intersecting constraints" << std::endl;
  assert(false);
  return Vertex_handle();
}

template <class Tr, class I_tag >
Constrained_triangulation_plus_2<Tr,I_tag>:: Vertex_handle 
Constrained_triangulation_plus_2<Tr,I_tag>::
intersect(Face_handle f, int i, 
	  Vertex_handle vaa,
	  Vertex_handle vbb,
	  Tag_true)
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
  assert(hierarchy.enclosing_constraint(vcc,vdd,vc,vd));
  assert(hierarchy.enclosing_constraint(vaa,vbb,va,vb));
						  
  Point pi; //creator for point is required here
  Object result;
  result = intersection(Segment(vc->point(),vd->point()),
			Segment(va->point(),vb->point()));
  assert(assign(pi, result));

  Vertex_handle vi = special_insert_in_edge(pi,f,i);
  hierarchy.split_constraint(vcc,vdd,vi);
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


CGAL_END_NAMESPACE

#endif //CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H



