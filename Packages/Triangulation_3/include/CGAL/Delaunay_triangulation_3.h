// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Delaunay_triangulation_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================


#ifndef CGAL_DELAUNAY_TRIANGULATION_3_H
#define CGAL_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <set>
#include <list>
#include <algorithm>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Delaunay_remove_tds_3.h>
#include <CGAL/Triangulation_face_base_2.h>
CGAL_BEGIN_NAMESPACE





template < class Gt, class Tds>
class Delaunay_triangulation_3 : public Triangulation_3<Gt,Tds>
{
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream& is, Triangulation_3<Gt,Tds> &tr);

public:
  typedef Tds Triangulation_data_structure;
  typedef Gt  Geom_traits;

  typedef typename Gt::Point_3 Point;
  typedef typename Gt::Vector_3 Vector;
  typedef typename Gt::Segment_3 Segment;
  typedef typename Gt::Triangle_3 Triangle;
  typedef typename Gt::Tetrahedron_3 Tetrahedron;

  // Function objects
  typedef typename Gt::Side_of_oriented_sphere_3 Side_of_oriented_sphere;
  typedef typename Gt::Coplanar_side_of_oriented_circle_3 Coplanar_side_of_oriented_circle;
  typedef typename Gt::Cross_product Cross_product;

  typedef typename Triangulation_3<Gt,Tds>::Cell_handle Cell_handle;
  typedef typename Triangulation_3<Gt,Tds>::Vertex_handle Vertex_handle;

  typedef typename Triangulation_3<Gt,Tds>::Cell Cell;
  typedef typename Triangulation_3<Gt,Tds>::Vertex Vertex;
  typedef typename Triangulation_3<Gt,Tds>::Facet Facet;
  typedef typename Triangulation_3<Gt,Tds>::Edge Edge;

  typedef typename Triangulation_3<Gt,Tds>::Cell_circulator Cell_circulator;
  typedef typename Triangulation_3<Gt,Tds>::Cell_iterator Cell_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Facet_iterator Facet_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Edge_iterator Edge_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Vertex_iterator Vertex_iterator;

  typedef typename Triangulation_3<Gt,Tds>::Locate_type Locate_type;

protected:
  Side_of_oriented_sphere side_of_oriented_sphere;
  Coplanar_side_of_oriented_circle coplanar_side_of_oriented_circle;
  Cross_product cross_product;

public:

  Delaunay_triangulation_3()
    : Triangulation_3<Gt,Tds>() {
    init_function_objects();
  }
  
  Delaunay_triangulation_3(const Gt & gt)
    : Triangulation_3<Gt,Tds>(gt) {
    init_function_objects();
  }
  
  Delaunay_triangulation_3(const Point & p0,
			   const Point & p1,
			   const Point & p2,
			   const Point & p3)
    : Triangulation_3<Gt,Tds>(p0,p1,p2,p3){
  init_function_objects();
  } // debug

  // copy constructor duplicates vertices and cells
  Delaunay_triangulation_3(const Delaunay_triangulation_3<Gt,Tds> & tr)
    : Triangulation_3<Gt,Tds>(tr)
    { 
      init_function_objects();
      CGAL_triangulation_postcondition( is_valid() );  
    }
  

  void init_function_objects() 
    {
      side_of_oriented_sphere = geom_traits().side_of_oriented_sphere_3_object();
      coplanar_side_of_oriented_circle = geom_traits().coplanar_side_of_oriented_circle_3_object();
      cross_product = geom_traits().cross_product_object();
    }


  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }

  Vertex_handle insert(const Point & p );

  Vertex_handle insert(const Point & p, Cell_handle start);
  
  bool remove(Vertex_handle v );

  // to be made private after tests
  void make_hole_3D(Vertex_handle v,
                    std::set<Facet> & boundhole,
                    std::set<Vertex_handle> & boundvert,
                    std::list<Facet> & outside,
                    std::list<Facet> & inside); 
  
  void undo_make_hole_3D(std::list<Facet> & outside,
                         std::list<Facet> & inside);
  
  void delete_cells(std::list<Cell_handle> & cells, int dummy_for_windows=0);
  void delete_cells(std::list<Facet> & cells);

  bool fill_hole_3D(std::set<Facet> & boundhole,
	      std::set<Vertex_handle> & boundvert);


  bool remove_3D_ear(Vertex_handle v );
  bool fill_hole_3D_ear(std::list<Facet> & boundhole);
  void make_hole_3D_ear( Vertex_handle v, 
	                 std::list<Facet> & boundhole,
			 //	                 std::set<Vertex_handle> & boundvert,
	                 std::list<Cell_handle> & hole);
  void undo_make_hole_3D_ear(std::list<Facet> & boundhole,
		             std::list<Cell_handle> & hole);
  void print(Vertex_handle v) const;
  void print(Cell_handle c) const;

private:
#ifdef SYL
  typedef std::vector<void *> Conflict_set;
#else
  typedef std::set<void *> Conflict_set;
#endif

  void
  find_conflicts_3(Conflict_set & conflicts, const Point & p,
		   Cell_handle c, Cell_handle & ac, int & i);
    // 3d case
    // p is in conflict with c
    // finds the set conflicts of cells in conflict with p
    // gives a cell ac having a facet on the boundary of conflicts
    // and the index i of its facet on the boundary
  void
  find_conflicts_2(Conflict_set & conflicts, const Point & p,
		   Cell_handle c, Cell_handle & ac, int & i);
    // 2d case
    // p is in conflict with c
    // finds the set conflicts of cells in conflict with p
    // gives a cell ac having a facet on the boundary of conflicts
    // and the index i of its facet on the boundary
  bool
  violates( Vertex_handle u, 
	    Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, 
	    Facet f);
  // checks whether edge uv crosses f

  bool remove_3D(Vertex_handle v );
  //  void make_hole_3D(Vertex_handle v,
  //		    std::list<Facet> & hole) const;		 

  Oriented_side  
  side_of_oriented_circle(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & test) const;

  Bounded_side
  side_of_sphere( Vertex_handle v0, 
		  Vertex_handle v1, 
		  Vertex_handle v2, 
		  Vertex_handle v3, 
		  const Point & p) const;

  Bounded_side
  side_of_sphere_inf(const Point & p0, 
		     const Point & p1, 
		     const Point & p2, 
		     const Point & p) const;

public:

  Bounded_side
  side_of_sphere( Cell_handle c, const Point & p) const;

  Bounded_side
  side_of_circle( const Facet & f, const Point & p) const
    {
      return side_of_circle(f.first, f.second, p);
    }

  Bounded_side
  side_of_circle( Cell_handle c, int i, const Point & p) const;
    // precondition : dimension >=2
    // in dimension 3, - for a finite facet
    // returns ON_BOUNDARY if the point lies on the circle,
    // ON_UNBOUNDED_SIDE when exterior, ON_BOUNDED_SIDE
    // interior
    // for an infinite facet, considers the plane defined by the
    // adjacent finite facet of the same cell, and does the same as in 
    // dimension 2 in this plane
    // in dimension 2, for an infinite facet
    // in this case, returns ON_BOUNDARY if the point lies on the 
    // finite edge (endpoints included) 
    // ON_BOUNDED_SIDE for a point in the open half-plane
    // ON_UNBOUNDED_SIDE elsewhere
  
  bool is_valid(bool verbose = false, int level = 0) const;

  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

};

template < class Gt, class Tds >
Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
insert(const Point & p )
{
  Cell_handle start;
  if ( dimension() >= 1 ) {
    // there is at least one finite "cell" (or facet or edge)
    start = infinite_vertex()->cell()
      ->neighbor( infinite_vertex()->cell()->index( infinite_vertex()) );
  }
  else {
    start = NULL;
  }
  return insert( p, start );
}

template < class Gt, class Tds >
Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
insert(const Point & p, Cell_handle start)
{
  switch (dimension()) {
  case 3:
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, start, lt, li, lj);
      if ( lt == VERTEX ) return c->vertex(li);
      else {
//       case OUTSIDE_CONVEX_HULL:
//       case CELL:
//       case FACET:
//       case EDGE:
	  Vertex_handle v = new Vertex(p);
	  set_number_of_vertices(number_of_vertices()+1);
	  Conflict_set conflicts;
#ifdef SYL
	  conflicts.reserve(32);
#endif
	  Cell_handle aconflict;
	  int ineighbor;
	  find_conflicts_3(conflicts,p,c,aconflict,ineighbor);
	  _tds.star_region(conflicts,&(*v),&(*aconflict),ineighbor);
	  return v;
      }
      //      break;// unreachable code
    }// dim 3
  case 2:
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, start, lt, li, lj);
      switch (lt) {
      case OUTSIDE_CONVEX_HULL:
      case CELL:
      case FACET:
      case EDGE:
	{
	  Vertex_handle v = new Vertex(p);
	  set_number_of_vertices(number_of_vertices()+1);
	  Conflict_set conflicts;
#ifdef SYL
	  conflicts.reserve(16);
#endif
	  Cell_handle aconflict;
	  int ineighbor;
	  find_conflicts_2(conflicts,p,c,aconflict,ineighbor);
	  _tds.star_region(conflicts,&(*v),&(*aconflict),ineighbor);
	  return v;
	}
      case VERTEX:
	return c->vertex(li);
      case OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Delaunay, the 3d
	  // triangulation will be Delaunay
	  return
	    Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p); 
	}
      }
    }//dim 2
  default :
    // dimension <= 1
    return Triangulation_3<Gt,Tds>::insert(p);
  }
  // return Triangulation_3<Gt,Tds>::insert(p);// to avoid warning with egcs // to avoid warning with bcc
}// insert(p)

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));
  CGAL_triangulation_precondition( _tds.is_vertex(&(*v)) );

  if ( dimension() <3 ) {
    // to be implemented : removal in degenerate dimensions
    // the traingulation is now rebuilt...

    Vertex_iterator vit;
    Vertex_iterator vdone = vertices_end();
    std::list<Point> points;
    for ( vit = finite_vertices_begin(); vit != vdone ; ++vit) {
      if ( v != (*vit).handle() ) 
	{ points.push_front( vit->point() ); }
    }
    typename std::list<Point>::iterator pit;
    typename std::list<Point>::iterator pdone = points.end();
    
    clear();
    for ( pit = points.begin(); pit != pdone; ++pit) {
      insert( *pit );
    }

    return true;
  }

  return remove_3D_ear(v);
}// remove(v)

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
remove_3D(Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 3 );

  std::list<Facet> boundary; // facets on the boundary of the hole
  std::set<Vertex_handle> bdvert; // vertices on the boundary of the hole

  if ( test_dim_down(v) ) {

    // _tds.remove_dim_down(&(*v)); return; 
    // !!!!!!!!!!!! TO BE DONE !!!!!!!!!!!!!
    // the triangulation is rebuilt...

    Vertex_iterator vit;
    Vertex_iterator vdone = vertices_end();
    std::list<Point> points;
    for ( vit = finite_vertices_begin(); vit != vdone ; ++vit) {
	if ( v != (*vit).handle() ) { points.push_front( vit->point() ); }
    }
    typename std::list<Point>::iterator pit;
    typename std::list<Point>::iterator pdone = points.end();
    
    clear();
    for ( pit = points.begin(); pit != pdone; ++pit) {
      insert( *pit );
    }
    return true;
  }

  //  std::cout << "removed point " << v->point() << std::endl;

  std::list<Facet> outside, inside;
  make_hole_3D(v, boundary, bdvert, outside, inside);
  bool filled = fill_hole_3D(boundary, bdvert);
  if(filled){
    delete( &(*v) );
    delete_cells(inside);
    set_number_of_vertices(number_of_vertices()-1);
  } else {
    undo_make_hole_3D(outside, inside);
  }
  return filled;
  
}// remove_3D(v)

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
make_hole_3D( Vertex_handle v, 
              std::set<Facet> & boundhole,
              std::set<Vertex_handle> & boundvert,
                    std::list<Facet> & outside,
                    std::list<Facet> & inside)
{
  CGAL_triangulation_precondition( ! test_dim_down(v) );

  typedef std::set<Cell_handle> Hole_cells;
  Hole_cells cells;
  incident_cells( v, cells );
  int i, indv;

  typename Hole_cells::iterator cit = cells.begin();
  typename Hole_cells::iterator cdone = cells.end();

  Cell_handle opp_cit;
  Vertex_handle vi;
  do {

    indv = (*cit)->index(&(*v));
    opp_cit = (*cit)->neighbor( indv );
    
    boundhole.insert (  std::make_pair( opp_cit, opp_cit->index(*cit)) );
    outside.push_back(  std::make_pair( opp_cit, opp_cit->index(*cit)) );
    inside.push_back(   std::make_pair( Cell_handle(*cit), indv) );
    for ( i=0; i<4; i++) {
      if ( i != indv ) {
	vi = (*cit)->vertex(i);
	if ( boundvert.find( vi ) == boundvert.end() )
	  {
	    boundvert.insert( vi );
	    vi->set_cell( opp_cit );
	  }
      }
    }

    ++cit;
  } while ( cit != cdone );

  return;
}// make_hole_3D







template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
undo_make_hole_3D(std::list<Facet> & outside,
                  std::list<Facet> & inside){
  std::list<Facet>::iterator cit = inside.begin();
  for(std::list<Facet>::iterator fit = outside.begin(); 
      fit != outside.end();
      fit++) {
    Cell_handle ch = (*fit).first;
    ch->set_neighbor((*fit).second, (*cit).first);
    CGAL_triangulation_assertion( (*cit).first->neighbor((*cit).second) == ch );
    for(int i = 0; i < 4; i++) {
      ch->vertex(i)->set_cell(ch);
    }
    cit++;
  }
}


template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
delete_cells(std::list<Cell_handle> & hole, int dummy_for_windows) {
  for(std::list<Cell_handle>::iterator cit = hole.begin(); cit != hole.end(); cit++) {
    _tds.delete_cell( &*(*cit) );
  }
}


template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
delete_cells(std::list<Facet> & hole) {
    for(std::list<Facet>::iterator cit = hole.begin(); cit != hole.end(); cit++) {
     _tds.delete_cell( &*((*cit).first) );
   }
}




//debug
template < class Gt, class Tds >
void 
Delaunay_triangulation_3<Gt,Tds>::
print( Vertex_handle v ) const
{
  if ( is_infinite( v ) )
    std::cout << "inf" << "; ";
  else 
    std::cout << v->point() << "; ";
}


template < class Gt, class Tds >
void 
Delaunay_triangulation_3<Gt,Tds>::
print( Cell_handle c ) const
{
  for(int i = 0; i < 4; i++) {
    print(c->vertex(i));
  }
  std::cout << std::endl;
}



template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
fill_hole_3D( std::set<Facet> & boundhole,
	      std::set<Vertex_handle> & boundvert)
  // examines each facet of the hole in turn and finds the fourth
  // vertex to build a cell
{
  typename std::set<Facet>::iterator fit;
  typename std::set<Facet>::iterator fitset_tmp;

  std::list<Cell_handle> cells;
  Facet fit_stick;
  typename std::list<Facet>::iterator fitlist_tmp;

  typename std::set<Vertex_handle>::iterator vit;
  typename std::set<Vertex_handle>::iterator vdone = boundvert.end();

  Vertex_handle v[3]; // current facet
  Vertex_handle vf1,vf2,vf3;

  std::set<Vertex_handle> oppvert; // vertices of the hole
  // that are all copsherical with the current facet and whose sphere
  // is empty 
  typename std::set<Vertex_handle>::iterator bv;

  std::list<Facet> cosph_bound; // facets of boundhole that
  // are cospherical with the vertices of oppvert 
  std::list<Facet> not_violate;
  std::list<Facet> for_next; // facets created or not used by
  // the current vertex 

  // the facets in these sets will always be given by the cell
  // **outside** the hole and their index in this cell

  Cell_handle cnew;
  int i;

  Bounded_side sos;
  Orientation ori;
  bool opp_inf, sticked; 
  //bool created;
  int nbnew = 1; // to detect a loop in the execution due to an
  // impossible case 

  //debug
  //  int nbfacet = 0;

  while( ! boundhole.empty() ) {
    //    std::cout << boundhole.size() << " boundary faces remaining" 
    // << std::endl;
    fit = boundhole.begin();
//     for ( fit = boundhole.begin(); fit != boundhole.end(); ++fit ) {
// 	print( (*fit).first->vertex( ((*fit).second+1)&3 ) );
// 	print( (*fit).first->vertex( ((*fit).second+2)&3 ) );
// 	print( (*fit).first->vertex( ((*fit).second+3)&3 ) );
// 	std::cout << std::endl;
//     }	       

    fit = boundhole.begin();
    while( (fit != boundhole.end()) && has_vertex(*fit,infinite_vertex()) ) {
      ++fit;
    }
    if (fit == boundhole.end()) {
      // all the remaining facets have the infinite vertex 
      //      std::cout << " pb" << std::endl;
    }

    cosph_bound.clear();
    oppvert.clear();
    //    new_fac.clear();

    //debug
    //    nbfacet++;
    
    // orientation of the facet chosen so that the future cell
    // v0,v1,v2,oppvert is correctly oriented
    // (ie. orientation reversed with respect to the cell giving (*fit),
    // which is outside the hole)

    if ( ((*fit).second%2) == 0) {
      v[0] = (*fit).first->vertex( ((*fit).second+1)&3 );
      v[1] = (*fit).first->vertex( ((*fit).second+2)&3 );
      v[2] = (*fit).first->vertex( ((*fit).second+3)&3 );
    }
    else {
      v[0] = (*fit).first->vertex( ((*fit).second+2)&3 );
      v[1] = (*fit).first->vertex( ((*fit).second+1)&3 );
      v[2] = (*fit).first->vertex( ((*fit).second+3)&3 );
    }      

    //debug
//     std::cout << std::endl << nbfacet << "- facet : " << std::endl;
//     for ( i=0; i<3; i++ ) { print( v[i] ); }
//     std::cout << " index " << (*fit).second << std::endl;
//     std::cout << "   cell " ;
//     for ( i=0; i<4; i++ ) { print( (*fit).first->vertex(i) ); }
//     std::cout << std::endl;

    if ( nbnew == 0 ) {
      // the program is looping
      std::cerr << " !!!!! IMPOSSIBLE TO RETRIANGULATE !!!!! " << std::endl;
      delete_cells(cells);
      return false;
    }
    nbnew = 0;
   
    // looking for a vertex to build a cell with the current facet *fit
    vit = boundvert.begin();
    while ( (( *vit == v[0] ) || ( *vit == v[1] ) || ( *vit == v[2] )) 
	    ||
	    ( ( ! is_infinite(v[0]) ) && 
	      ( ! is_infinite(v[1]) ) &&
	      ( ! is_infinite(v[2]) ) &&
	      ( ! is_infinite(*vit) ) &&
	      orientation( v[0]->point(),v[1]->point(),
					 v[2]->point(),(*vit)->point() ) 
	      != POSITIVE ) ) 
      { ++vit; }
    oppvert.insert(*vit); // candidate to build a cell with *fit
    // it is either infinite or on the positive side of the facet
    opp_inf = is_infinite( *vit );

    // checking whether the current candidate forms an empty sphere
    ++vit;
    i=0;
    while ( vit != vdone ) {
      if ( *vit == v[0] ) {++vit; continue;};
      if ( *vit == v[1] ) {++vit; continue;};
      if ( *vit == v[2] ) {++vit; continue;};

      //debug
//       std::cout << "    vertex " ; print( *oppvert.begin() );
//       std::cout << std::endl << "            test " ; print( *vit );
//       std::cout << std::endl;

      if ( ! is_infinite(*vit) ) {
	sos = side_of_sphere(v[0],v[1],v[2],*oppvert.begin(), 
			     (*vit)->point());
	ori = orientation( v[0]->point(), v[1]->point(),
					v[2]->point(), (*vit)->point() );
	if ( ori == POSITIVE ) {
	  if ( sos == ON_BOUNDED_SIDE ) {
	    if ( is_infinite(*oppvert.begin()) ) {
	      oppvert.erase(infinite_vertex());
	      opp_inf = false;
	      //	      std::cout << "          cocy kept - inf"
	      //	      << std::endl; 
	    }
	    else {
	      for ( bv=oppvert.begin(); bv!= oppvert.end(); ++bv ) {
		if ( side_of_sphere(v[0],v[1],v[2],(*bv),(*vit)->point())
		     != ON_BOUNDARY ) break;
	      }
	      if ( bv == oppvert.end() ) {
		oppvert.erase(*oppvert.begin());
		// because all the other points are cospherical with vit
		//		std::cout << "          cocy kept" <<
		//		std::endl; 
	      }
	      else {
		oppvert.clear();
	      }
	    }
	    oppvert.insert(*vit); 
	    //	    std::cout << "          inside" << std::endl;
	  }
	  else {
	    if ( sos == ON_BOUNDARY ) {
	      oppvert.insert(*vit); 
	      //	      std::cout << "          cosph" << std::endl;
	    }; 
	  };
	}
	else {
	  if ( (ori == COPLANAR) && (opp_inf) && (sos == ON_BOUNDARY) ) {
	    oppvert.insert(*vit); 
	    //	    std::cout << "          cocycl" << std::endl;
	  }
	}
      };

      ++vit;
    }
    // now oppvert contains 
    // either the only possible vertex to build a cell with the
    // current facet *fit
    // or the set of vertices that are cospherical with the vertices
    // of the current facet *fit

    // in the case when there are several vertices in oppvert, the
    // polyhedron they form must be triangulated without violating the 
    // possibly already existing facets on the boundary of the hole.

    // creation of a list of cospherical facets already on the
    // boundary of the hole, that must not be violated
    if ( oppvert.size() > 1 ) {
      fitset_tmp = boundhole.begin();
      //debug
      std::cout << "    constraining facets " << std::endl;
      
      not_violate.clear();

      while ( fitset_tmp != boundhole.end() ) {
	vf1 = (*fitset_tmp).first->vertex(((*fitset_tmp).second+1)&3);
	vf2 = (*fitset_tmp).first->vertex(((*fitset_tmp).second+2)&3);
	vf3 = (*fitset_tmp).first->vertex(((*fitset_tmp).second+3)&3);
	if (
	    ( (vf1 == v[0]) || (vf1 == v[1]) || (vf1 == v[2]) ||
	      ( oppvert.find(vf1) != oppvert.end() ) )
	    &&
	    ( (vf2 == v[0]) || (vf2 == v[1]) || (vf2 == v[2]) ||
	      ( oppvert.find(vf2) != oppvert.end() ) )
	    &&
	    ( (vf3 == v[0]) || (vf3 == v[1]) || (vf3 == v[2]) ||
	      ( oppvert.find(vf3) != oppvert.end() ) )
	    )
	  {
	    cosph_bound.push_front( *fitset_tmp );
	    // fit will be one of the facets inserted

	    //	    print( vf1 ); print( vf2 ); print( vf3 ); 
	    //	    std::cout << " cosph" << std::endl;
	    ++fitset_tmp;
	    continue;
	  }
	for ( bv=oppvert.begin(); bv!=oppvert.end(); ++bv ) {
	  if ( ( ! is_infinite(*bv) )
	       &&
	       ( ! is_infinite( *fitset_tmp ) )
	       &&
	       ((( (vf1 == v[0]) || (vf1 == v[1]) || (vf1 == v[2]) )
		 &&
		 ( (vf2 == v[0]) || (vf2 == v[1]) || (vf2 == v[2]) )
		 &&
		 ( orientation
		   (vf1->point(),vf2->point(),
		    vf3->point(),(*bv)->point()) == COPLANAR ))
		||
		(( (vf1 == v[0]) || (vf1 == v[1]) || (vf1 == v[2]) )
		 &&
		 ( (vf3 == v[0]) || (vf3 == v[1]) || (vf3 == v[2]) )
		 &&
		 ( orientation
		   (vf1->point(),vf2->point(),
		    vf3->point(),(*bv)->point()) == COPLANAR ))
		||
		(( (vf2 == v[0]) || (vf2 == v[1]) || (vf2 == v[2]) )
		 &&
		 ( (vf3 == v[0]) || (vf3 == v[1]) || (vf3 == v[2]) )
		 &&
		 ( orientation
		   (vf1->point(),vf2->point(),
		    vf3->point(),(*bv)->point()) == COPLANAR )) ) )
	    {
	      not_violate.push_front(*fitset_tmp);
	      //	      print( vf1 ); print( vf2 ); print( vf3 ); 
	      //	      std::cout << " copl ";
	      //	      print(*bv); std::cout << std::endl;
	      ++fitset_tmp;
	      continue;
	    }
	}
	++fitset_tmp;
      }
    }
    else cosph_bound.push_front( *fit ); // *fit is the only element in
    // this case

    not_violate.insert(not_violate.end(),
		       cosph_bound.begin(),cosph_bound.end());

    // if there are cospherical vertices, the polyhedron they form is
    // triangulated by linking each vertex in its turn to its visible
    // facets, without violating the facets of cosph_bound
    while ( oppvert.size() > 0 ) {
      // look for all the facets of the polyhedron visible from
      // oppvert.begin() 
//       std::cout << std::endl << "oppvert " ; print(*oppvert.begin());
//       std::cout << std::endl;

      for_next.clear();

      while ( ! cosph_bound.empty() ) {
	fit_stick = cosph_bound.front();
	cosph_bound.pop_front(); 
	

	if ( has_vertex( fit_stick, *oppvert.begin() ) ) {
	  for_next.push_front( fit_stick );
	  continue;
	}
	     
	if ( ((fit_stick).second%2) == 0 ) {
	  v[0] = (fit_stick).first->vertex( ((fit_stick).second+1)&3 );
	  v[1] = (fit_stick).first->vertex( ((fit_stick).second+2)&3 );
	  v[2] = (fit_stick).first->vertex( ((fit_stick).second+3)&3 );
	}
	else {
	  v[0] = (fit_stick).first->vertex( ((fit_stick).second+2)&3 );
	  v[1] = (fit_stick).first->vertex( ((fit_stick).second+1)&3 );
	  v[2] = (fit_stick).first->vertex( ((fit_stick).second+3)&3 );
	}

	//	std::cout << "test facet " ;
	//	for ( i=0; i<3; i++ ) {print(v[i]);}; std::cout << std::endl;

	// testing if the current vertex of oppvert can be associated
	// with this facet fit_stick without violating any facet of
	// not_violate
	fitlist_tmp = not_violate.begin();
	while ( fitlist_tmp != not_violate.end() ) {
 	  if ( (*fitlist_tmp) == fit_stick ) {
 	    ++fitlist_tmp;
 	    continue;
 	  }
	  if ( ! (is_infinite(v[0]) ||
		  is_infinite(v[1]) || 
		  is_infinite(v[2]) ||
		  is_infinite(*oppvert.begin()) ) 
	       ) {
	    if ( violates( *oppvert.begin(),v[0],v[1],v[2],*fitlist_tmp ) ) 
	      break;
	  }
	  ++fitlist_tmp;
	}
	if ( fitlist_tmp != not_violate.end() ) {
	  // this facet cannot form a cell with oppvert.begin()
	  for_next.push_front(fit_stick);
	  continue;
	}
	
	// here fit_stick does not violate any facet
	// test whether it is visible from oppvert.begin()
	//created = false;
	if ( ( ! (is_infinite(v[0]) ||
		  is_infinite(v[1]) || 
		  is_infinite(v[2])) )
	     // already tested, ti be removed
	     && 
	     ( ( ! is_infinite((*oppvert.begin())) &&
		 ( orientation
		   (v[0]->point(),v[1]->point(),v[2]->point(),
		    ((*oppvert.begin()))->point()) 
		   == POSITIVE ) )
	       ||
	       ( is_infinite((*oppvert.begin())) &&
		 ! is_infinite( (fit_stick).first
				->vertex((fit_stick).second ) ) )
	       )
	     ) 
	  { // creation of a cell
	    cnew = new Cell (v[0],v[1],v[2],(*oppvert.begin()),
			     NULL,NULL,NULL,(fit_stick).first);
	    ++nbnew;
// 	    std::cout << "cell " ;
// 	    for (i=0;i<4;i++) { print(cnew->vertex(i)); };
// 	    std::cout << std::endl;

	    _tds.add_cell( &(*cnew) );
	    (fit_stick).first->set_neighbor((fit_stick).second,cnew);
	    //created = true;
	    ((*oppvert.begin()))->set_cell(cnew);

	    boundhole.erase( fit_stick );
	    not_violate.remove( fit_stick);
	    for_next.remove( fit_stick);

	    for ( i=0; i<3; i++ ) { 
	      v[i]->set_cell(cnew); //useless ? to be checked

	      // if the facet ( cnew, i ) is already on the boundary 
	      // of the hole, it must be sticked and removed from
	      // the hole
// 	      std::cout << "    facet " ;
// 	      for ( j=1; j<4; j++ ) {print(cnew->vertex((i+j)&3));};
// 	      std::cout << std::endl;

	      fitlist_tmp = cosph_bound.begin();
	      sticked = false;
	      while ( (fitlist_tmp != cosph_bound.end()) && (! sticked) ) {
		if ( are_equal( *fitlist_tmp, cnew, i ) ) {
		  (*fitlist_tmp).first->set_neighbor
		    ((*fitlist_tmp).second,cnew);
		  cnew->set_neighbor(i,(*fitlist_tmp).first);

		  not_violate.remove(*fitlist_tmp);
		  boundhole.erase(*fitlist_tmp);
		  cosph_bound.erase(fitlist_tmp);

		  sticked = true;
		  //		  std::cout << "      sticked" << std::endl;
		  break;
		}
		++fitlist_tmp;
	      }
	      fitlist_tmp = for_next.begin();
	      while ( (fitlist_tmp != for_next.end()) && (! sticked) ) {
		if ( are_equal( *fitlist_tmp, cnew, i ) ) {
		  (*fitlist_tmp).first->set_neighbor
		    ((*fitlist_tmp).second,cnew);
		  cnew->set_neighbor(i,(*fitlist_tmp).first);

		  not_violate.remove(*fitlist_tmp);
		  boundhole.erase(*fitlist_tmp);
		  for_next.erase(fitlist_tmp);

		  sticked = true;
		  //		  std::cout << "      sticked" << std::endl;
		  break;
		}
		++fitlist_tmp;
	      }
	      if ( ! sticked ) {
		//		std::cout << "      not sticked" << std::endl;

		for_next.push_front( std::make_pair( cnew, i ) );
		not_violate.push_front( std::make_pair( cnew, i ) );
	      }
	    }
	  }
	else 
	  {
	    // fit_stick is already in not_violate
	    for_next.push_front( fit_stick );
	  }
      }
      cosph_bound = for_next;
      oppvert.erase(oppvert.begin());
    }



//     std::cout << std::endl << "triangulated polyhedron" << std::endl;
//     std::cout << boundhole.size() << " facets in boundhole" <<
//       std::endl;
//     std::cout << cosph_bound.size() << " facets remaining in cosph_bound" 
//       << std::endl;
    // looking whether the new facets just built are already on the
    // boundary 
    // if they are, update the adjacency relations
    //
    // unefficient code, to be changed 

    while ( ! cosph_bound.empty() ) {
      fit_stick = cosph_bound.front();
      cosph_bound.pop_front();

//       std::cout << "   facet " ;
//       for ( i=1; i<4; i++ ) 
// 	{print((fit_stick).first->vertex(((fit_stick).second+i)&3));};
//       std::cout << " index " << (fit_stick).second << std::endl;
      
      sticked = false;
      fitset_tmp = boundhole.begin();
      while (fitset_tmp != boundhole.end()) {
	if ((*fitset_tmp) == fit_stick) { ++fitset_tmp; continue; };
	if ( are_equal( *fitset_tmp, fit_stick ) ) {
// 	  std::cout << "        sticked to " ;
// 	  for ( i=1; i<4; i++ ) {
// 	    print((*fitset_tmp).first->vertex(((*fitset_tmp).second+i)&3));};
// 	  std::cout << " of " ;
// 	  for ( i=0; i<4; i++ ) {
// 	    print((*fitset_tmp).first->vertex(i));}
// 	  std::cout << std::endl;

	  (*fitset_tmp).first->set_neighbor
	    ((*fitset_tmp).second,(fit_stick).first);
	  (fit_stick).first->set_neighbor
	    ((fit_stick).second,(*fitset_tmp).first);

	  boundhole.erase(*fitset_tmp);

	  sticked = true;
	  break;
	};
	++fitset_tmp;
      }
      if ( ! sticked ) {
	boundhole.insert( fit_stick );
	// 	std::cout << "        not sticked" << std::endl;
      }
    }

  }

  return true;
}// fill_hole_3D


template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
find_conflicts_3(Conflict_set & conflicts, const Point & p,
		 Cell_handle c, Cell_handle & ac, int & i)
  // 3d case
  // p is in conflict with c
  // finds the set conflicts of cells in conflict with p
  // gives a cell ac having a facet on the boundary of conflicts
  // and the index i of its facet on the boundary
{
#ifdef SYL
  conflicts.push_back( (Conflict_set::value_type) &(*c) );
  c->set_flags(1);
#else
  (void) conflicts.insert( (Conflict_set::value_type) &(*c) );
#endif

  for ( int j=0; j<4; j++ ) {
    Cell_handle test = c->neighbor(j);
#ifdef SYL
    if (test->get_flags() == 1)
#else
    if (conflicts.find( (Conflict_set::value_type) &(*test) ) != conflicts.end())
#endif
      continue; // test was already tested and found to be in conflict.
    if ( side_of_sphere( test, p ) == ON_BOUNDED_SIDE )
      find_conflicts_3(conflicts, p, test, ac, i);
    else {
      ac = c;
      i = j;
    }
  }
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
violates( Vertex_handle u, 
	  Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, 
	  Facet f)
{
    // u, v0, v1, v2 supposed to be all different
    const Point & pu = u->point();

    const Point * pf[3] = {
      &(f.first)->vertex((f.second+1)&3)->point(),
      &(f.first)->vertex((f.second+2)&3)->point(),
      &(f.first)->vertex((f.second+3)&3)->point()};
 
    if ( orientation( pu, *pf[0], *pf[1], *pf[2] ) != COPLANAR ) 
      return false;

    const Point * p[3] = {
      &v0->point(),
      &v1->point(),
      &v2->point()};

    Orientation o[3];
    o[0] = orientation(*p[0],*pf[0],*pf[1],*pf[2]);
    o[1] = orientation(*p[1],*pf[0],*pf[1],*pf[2]);
    o[2] = orientation(*p[2],*pf[0],*pf[1],*pf[2]);

    if ( ( o[0] != COPLANAR ) && ( o[1] != COPLANAR ) && ( o[2] != COPLANAR ) )
      return false;

    for (int i=0; i<3; i++ ) {
      if ( pu == *pf[i] ) {
	int j = (i+1)%3;
	int k = (i+2)%3;
	return 
	    ( (o[0] == COPLANAR) && (*p[0] != *pf[j]) && (*p[0] != *pf[k]) &&
	      ( coplanar_orientation(*pf[j],*pf[k],*pf[i],*p[0]) == NEGATIVE ) )
	    ||
	    ( (o[1] == COPLANAR) && (*p[1] != *pf[j]) && (*p[1] != *pf[k]) &&
	      ( coplanar_orientation(*pf[j],*pf[k],*pf[i],*p[1]) == NEGATIVE ) )
	    ||
	    ( (o[2] == COPLANAR) && (*p[2] != *pf[j]) && (*p[2] != *pf[k]) &&
	      ( coplanar_orientation(*pf[j],*pf[k],*pf[i],*p[2]) == NEGATIVE ) );
      }
    }

    // here pu is none of *pf[i]
    for (int i=0; i<3; i++ ) {
      if ( o[i] == COPLANAR )
	{
	  for (int l=0; l<3; l++ ) {
	    int j = (l+1)%3;
	    int k = (l+2)%3;
	    if ( *p[i] == *pf[l] ) {
	      if ( (pu != *pf[j]) && (pu != *pf[k]) &&
		   ( coplanar_orientation(*pf[j],*pf[k],*pf[l],pu) == NEGATIVE ) )
		return true;
	      else
		continue;
	    }
	    if ( coplanar_orientation(*pf[j],*pf[k],pu,*p[i]) != POSITIVE )
	      return true;
	  }
	}
      else
	continue;
    }

    return false;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
find_conflicts_2(Conflict_set & conflicts, const Point & p,
		 Cell_handle c, Cell_handle & ac, int & i)
  // 2d case
  // p is in conflict with c
  // finds the set conflicts of cells in conflict with p
  // gives a cell ac having a facet on the boundary of conflicts
  // and the index i of its facet on the boundary
{

#ifdef SYL
  conflicts.push_back( (Conflict_set::value_type) &(*c) );
  c->set_flags(1);
#else
  (void) conflicts.insert( (Conflict_set::value_type) &(*c) );
#endif

  for ( int j=0; j<3; j++ ) {
    Cell_handle test = c->neighbor(j);
#ifdef SYL
    if (test->get_flags() == 1)
#else
    if (conflicts.find( (Conflict_set::value_type) &(*test) ) != conflicts.end())
#endif
      continue;   // test was already found
    if ( side_of_circle( test, 3, p ) == ON_BOUNDED_SIDE )
      find_conflicts_2(conflicts, p, test, ac, i);
    else {
      ac = c;
      i = j;
    }
  }
}


template < class Gt, class Tds >
Oriented_side 
Delaunay_triangulation_3<Gt,Tds>::
side_of_oriented_circle(const Point & p,
			  const Point & q,
			  const Point & r,
			  const Point & test) const
{

  Vector v = cross_product(q-p, r-p);

  return coplanar_side_of_oriented_circle(p, q, r, test, v);
}



template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere(Cell_handle c, const Point & p) const
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int i3;
  if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
    return Bounded_side( side_of_oriented_sphere
			 (c->vertex(0)->point(),
			  c->vertex(1)->point(),
			  c->vertex(2)->point(),
			  c->vertex(3)->point(),p) );
//     Oriented_side
//       o = side_of_oriented_sphere(c->vertex(0)->point(),
// 						c->vertex(1)->point(),
// 						c->vertex(2)->point(),
// 						c->vertex(3)->point(),p);
//     return Bounded_side(o);
  // else infinite cell :
  int i0,i1,i2;
  if ( (i3%2) == 1 ) {
    i0 = (i3+1)&3;
    i1 = (i3+2)&3;
    i2 = (i3+3)&3;
  }
  else {
    i0 = (i3+2)&3;
    i1 = (i3+1)&3;
    i2 = (i3+3)&3;
  }
  Orientation
    o = orientation(c->vertex(i0)->point(),
				  c->vertex(i1)->point(),
				  c->vertex(i2)->point(),
				  p);
  if (o != ZERO)
    return Bounded_side(o);

  return Bounded_side( side_of_oriented_circle
		       ( c->vertex(i0)->point(), 
			 c->vertex(i1)->point(),
			 c->vertex(i2)->point(),
			 p ) );
//     {
//       Oriented_side s = 
// 	side_of_oriented_circle
// 	( c->vertex(i0)->point(),
// 	  c->vertex(i1)->point(),
// 	  c->vertex(i2)->point(),
// 	  p );
//       return ( (s == ON_NEGATIVE_SIDE) ? ON_UNBOUNDED_SIDE :
// 	       (s == ON_POSITIVE_SIDE) ? ON_BOUNDED_SIDE :
// 	       ON_BOUNDARY );
//     }
}// end side of sphere

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere(Vertex_handle v0, 
	       Vertex_handle v1, 
	       Vertex_handle v2, 
	       Vertex_handle v3, 
	       const Point & p) const
  // duplication of code with the other side_of_sphere should (?) be fixed 
{
  CGAL_triangulation_precondition( dimension() == 3 );
  Vertex_handle v[4] = { v0, v1, v2, v3 };
  int i;
  int i0,i1,i2;
  for ( i=0; i<4; i++) {
    if ( is_infinite( v[i] ) ) {
      if ( (i%2) == 1 ) {
	i0 = (i+1)&3;
	i1 = (i+2)&3;
	i2 = (i+3)&3;
      }
      else {
	i0 = (i+2)&3;
	i1 = (i+1)&3;
	i2 = (i+3)&3;
      }
      Orientation o = orientation(v[i0]->point(),
				  v[i1]->point(),
				  v[i2]->point(),
				  p);
      if (o != ZERO)
	return Bounded_side(o);

      return Bounded_side( side_of_oriented_circle
			   ( v[i0]->point(), 
			     v[i1]->point(),
			     v[i2]->point(),
			      p ) );
    }
  }
  
  // all vertices are finite

  return Bounded_side( side_of_oriented_sphere
		       ( v0->point(),   // remove the indirection through the array
			 v1->point(),
			 v2->point(),
			 v3->point(), p) );
}// end side of sphere



template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere_inf(const Point & p0, 
		   const Point & p1, 
		   const Point & p2, 
		   const Point & p) const
{
  CGAL_triangulation_precondition( dimension() == 3 );

  Orientation o = orientation(p0, p1, p2, p);
  if (o != ZERO)
    return Bounded_side(o);

  return Bounded_side( side_of_oriented_circle (p0, p1, p2, p ) );
}// end side of sphere


template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_circle(Cell_handle c, int i, const Point & p) const
  // precondition : dimension >=2
  // in dimension 3, - for a finite facet
  // returns ON_BOUNDARY if the point lies on the circle,
  // ON_UNBOUNDED_SIDE when exterior, ON_BOUNDED_SIDE
  // interior
  // for an infinite facet, considers the plane defined by the
  // adjacent finite facet of the same cell, and does the same as in 
  // dimension 2 in this plane
  // in dimension 2, for an infinite facet
  // in this case, returns ON_BOUNDARY if the point lies on the 
  // finite edge (endpoints included) 
  // ON_BOUNDED_SIDE for a point in the open half-plane
  // ON_UNBOUNDED_SIDE elsewhere
{
  CGAL_triangulation_precondition( dimension() >= 2 );
  int i3 = 5;

  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    // the triangulation is supposed to be valid, ie the facet
    // with vertices 0 1 2 in this order is positively oriented
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
      return Bounded_side( side_of_oriented_circle
			   (c->vertex(0)->point(),
			    c->vertex(1)->point(),
			    c->vertex(2)->point(),
			    p) );
//       Oriented_side
// 	o = side_of_oriented_circle(c->vertex(0)->point(),
// 						  c->vertex(1)->point(),
// 						  c->vertex(2)->point(),
// 						  p);
//       //       return ( (o == ON_NEGATIVE_SIDE) ? ON_UNBOUNDED_SIDE :
// 	       (o == ON_POSITIVE_SIDE) ? ON_BOUNDED_SIDE :
// 	       ON_BOUNDARY );
    // else infinite facet
    // v1, v2 finite vertices of the facet such that v1,v2,infinite
    // is positively oriented
    Vertex_handle 
      v1 = c->vertex( ccw(i3) ),
      v2 = c->vertex( cw(i3) );
    // does not work in dimension 3 :
    // we need a fourth point to look at orientation with respect
    // to v1v2
    Cell_handle n = c->neighbor(i3);
    Orientation o =
      coplanar_orientation( v1->point(), 
			    v2->point(), 
			    n->vertex(n->index(c))->point(),
			    p );
    if ( o != ZERO ) return Bounded_side( -o );
    // because p is in f iff
    // is does not lie on the same side of v1v2 as vn
    int i_e;
    Locate_type lt;
    // case when p collinear with v1v2
    return side_of_segment( p,
			    v1->point(), v2->point(),
			    lt, i_e );
//     switch (o) {
//     case POSITIVE:
//       // p lies on the same side of v1v2 as vn, so not in f
//       {
// 	return ON_UNBOUNDED_SIDE;
//       }
//     case NEGATIVE:
//       // p lies in f
//       { 
// 	return ON_BOUNDED_SIDE;
//       }
//     case ZERO:
//       // p collinear with v1v2
//       {
// 	int i_e;
// 	Locate_type lt;
// // 	Bounded_side side = 
// // 	  side_of_segment( p,
// // 			   v1->point(), v2->point(),
// // 			   lt, i_e );
// 	if ( side_of_segment( p,
// 			      v1->point(), v2->point(),
// 			      lt, i_e ) == ON_UNBOUNDED_SIDE )
// 	  // p lies on the line defined by the finite edge, but
// 	  // not in edge v1v2
// 	  return ON_UNBOUNDED_SIDE;
// 	// else p lies in edge v1v2 (including v1 or v2)
// 	return ON_BOUNDARY;
//       }
//     }// switch o
  }// dim 2

  // else dimension == 3
  CGAL_triangulation_precondition( (i >= 0) && (i < 4) );
  if ( ( ! c->has_vertex(infinite_vertex(),i3) ) || ( i3 != i ) ) {
    // finite facet
    // initialization of i0 i1 i2, vertices of the facet positively 
    // oriented (if the triangulation is valid)
    int i0 = (i>0) ? 0 : 1;
    int i1 = (i>1) ? 1 : 2;
    int i2 = (i>2) ? 2 : 3;
    CGAL_triangulation_precondition( orientation
				     (c->vertex(i0)->point(),
				      c->vertex(i1)->point(),
				      c->vertex(i2)->point(),
				      p) == COPLANAR );
    return Bounded_side ( side_of_oriented_circle
			  (c->vertex(i0)->point(),
			   c->vertex(i1)->point(),
			   c->vertex(i2)->point(),
			   p) );
//     Oriented_side
//       o = side_of_oriented_circle(c->vertex(i0)->point(),
// 						c->vertex(i1)->point(),
// 						c->vertex(i2)->point(),
// 						p);
//     return ( (o == ON_NEGATIVE_SIDE) ? ON_UNBOUNDED_SIDE :
// 	     (o == ON_POSITIVE_SIDE) ? ON_BOUNDED_SIDE :
// 	     ON_BOUNDARY );
  }

  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle 
    v1 = c->vertex( next_around_edge(i3,i) ),
    v2 = c->vertex( next_around_edge(i,i3) );
  Orientation o =
    coplanar_orientation( v1->point(),
			  v2->point(),
			  c->vertex(i)->point(),
			  p );
  // then the code is duplicated from 2d case
  if ( o != ZERO ) return Bounded_side( -o );
  // because p is in f iff 
  // it is not on the same side of v1v2 as c->vertex(i)
  int i_e;
  Locate_type lt;
  // case when p collinear with v1v2
  return side_of_segment( p,
			  v1->point(), v2->point(),
			  lt, i_e );
//   switch (o) {
//   case POSITIVE:
//     // p lies on the same side of v1v2 as c->vertex(i), so not in f
//     {
//       return ON_UNBOUNDED_SIDE;
//     }
//   case NEGATIVE:
//     // p lies in f
//     { 
//       return ON_BOUNDED_SIDE;
//     }
//   case ZERO:
//     // p collinear with v1v2
//     {
//       int i_e;
//       Locate_type lt;
// //       Bounded_side side = 
// // 	side_of_segment( p,
// // 			 v1->point(), v2->point(),
// // 			 lt, i_e );
//       if ( side_of_segment( p,
// 			    v1->point(), v2->point(),
// 			    lt, i_e ) == ON_UNBOUNDED_SIDE )
// 	// p lies on the line defined by the finite edge, but
// 	// not in edge v1v2
// 	return ON_UNBOUNDED_SIDE;
//       // else p lies in edge v1v2 (including v1 or v2)
//       return ON_BOUNDARY;
//     }
//   }// switch o
//   // end infinite facet
//   return ON_BOUNDARY; // to avoid warning with egcs
}// side_of_circle

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if ( ! tds().is_valid(verbose,level) ) {
    if (verbose) { std::cerr << "invalid data structure" << std::endl; }
    CGAL_triangulation_assertion(false); return false;
  }
    
  if ( &(*infinite_vertex()) == NULL ) {
    if (verbose) { std::cerr << "no infinite vertex" << std::endl; }
    CGAL_triangulation_assertion(false); return false;
  }

  int i;

  switch ( dimension() ) {
  case 3:
    {
      Cell_iterator it;
      for ( it = finite_cells_begin(); it != cells_end(); ++it ) {
	is_valid_finite((*it).handle());
	for ( i=0; i<4; i++ ) {
	  if ( side_of_sphere
	       ( (*it).handle(), 
		 it->vertex( (it->neighbor(i))->index((*it).handle() ) )
		 ->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty sphere " << std::endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	}
      }
      break;
    }
  case 2:
    {
      Facet_iterator it;
      for ( it = finite_facets_begin(); it != facets_end(); ++it ) {
	is_valid_finite((*it).first);
	for ( i=0; i<2; i++ ) {
	  if ( side_of_circle
	       ( (*it).first, 3,
		 (*it).first
		 ->vertex( (((*it).first)->neighbor(i))
			   ->index((*it).first) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty circle " << std::endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	}
      }
      break;
    }
  case 1:
    {
      Edge_iterator it;
      for ( it = finite_edges_begin(); it != edges_end(); ++it ) {
	is_valid_finite((*it).first);
      }
      break;
    }
  }
  if (verbose) { std::cerr << "Delaunay valid triangulation" << std::endl;}
  return true;
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_valid(Cell_handle c, bool verbose, int level) const
{
  int i;

  if ( ! (&(*c))->is_valid(dimension(),verbose,level) ) {
    if (verbose) { 
      std::cerr << "combinatorically invalid cell" ;
      for ( i=0; i <= dimension(); i++ ) {
	std::cerr << c->vertex(i)->point() << ", " ;
      }
      std::cerr << std::endl;
    }
    CGAL_triangulation_assertion(false); return false;
  }
  switch ( dimension() ) {
  case 3:
    {
      if ( ! is_infinite(c) ) {
	is_valid_finite(c,verbose,level);
	for ( i=0; i<4; i++ ) {
	  if ( side_of_sphere
	       ( c, c->vertex( (c->neighbor(i))->index(c) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty sphere " << std::endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	}
      }
      break;
    }
  case 2:
    {
      if ( ! is_infinite(c,3) ) {
	for ( i=0; i<2; i++ ) {
	  if ( side_of_circle
	       ( c , 3, c->vertex(c->neighbor(i)->index(c))->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty circle " << std::endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	}
      }
      break;
    }
  }
  if (verbose) { std::cerr << "Delaunay valid cell" << std::endl;}
  return true;
}



template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
remove_3D_ear(Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 3 );

  std::list<Facet> boundary; // facets on the boundary of the hole
  //  std::set<Vertex_handle> bdvert; // vertices on the boundary of the hole

  if ( test_dim_down(v) ) {

    // _tds.remove_dim_down(&(*v)); return; 
    // !!!!!!!!!!!! TO BE DONE !!!!!!!!!!!!!
    // the triangulation is rebuilt...

    Vertex_iterator vit;
    Vertex_iterator vdone = vertices_end();
    std::list<Point> points;
    for ( vit = finite_vertices_begin(); vit != vdone ; ++vit) {
	if ( v != (*vit).handle() ) { points.push_front( vit->point() ); }
    }
    typename std::list<Point>::iterator pit;
    typename std::list<Point>::iterator pdone = points.end();
    
    clear();
    for ( pit = points.begin(); pit != pdone; ++pit) {
      insert( *pit );
    }
    return true;
  }

  std::list<Cell_handle> hole;
  make_hole_3D_ear(v, boundary,
		   //bdvert, 
		   hole);


  bool filled = fill_hole_3D_ear(boundary);
  if(filled){
    delete( &(*v) );
    delete_cells(hole);
    set_number_of_vertices(number_of_vertices()-1);
  } else {
    undo_make_hole_3D_ear(boundary, hole);
  }

  return filled;

}// remove_3D_ear(v)



template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
make_hole_3D_ear( Vertex_handle v, 
	      std::list<Facet> & boundhole,
		  //std::set<Vertex_handle> & boundvert,
	      std::list<Cell_handle> & hole)
{
  CGAL_triangulation_precondition( ! test_dim_down(v) );

  typedef std::set<Cell_handle> Hole_cells;
  Hole_cells cells;
  incident_cells( v, cells );
  int i, indv;

  typename Hole_cells::iterator cit = cells.begin();
  typename Hole_cells::iterator cdone = cells.end();

  Cell_handle opp_cit;
  Vertex_handle vi;
  do {

    indv = (*cit)->index(&(*v));
    opp_cit = (*cit)->neighbor( indv );
    hole.push_back(*cit);    
    boundhole.push_back( std::make_pair( opp_cit, opp_cit->index(*cit)) );

    for ( i=0; i<4; i++) {
      if ( i != indv ) {
	vi = (*cit)->vertex(i);
	//if ( boundvert.find( vi ) == boundvert.end() )
	//  {
	//    boundvert.insert( vi );
	    vi->set_cell( opp_cit );
	    //  }
      }
    }

    ++cit;
  } while ( cit != cdone );

  return;
}// make_hole_3D_ear


template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
undo_make_hole_3D_ear(std::list<Facet> & boundhole,
		  std::list<Cell_handle> & hole){
  std::list<Cell_handle>::iterator cit = hole.begin();
  for(std::list<Facet>::iterator fit = boundhole.begin();	
      fit != boundhole.end();
      fit++) {
    Cell_handle ch = (*fit).first;
    ch->set_neighbor((*fit).second, *cit);
    for(int i = 0; i < 4; i++) {
      ch->vertex(i)->set_cell(ch);
    }
    cit++;
  }
}// undo_make_hole_3D_ear





template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
fill_hole_3D_ear( std::list<Facet> & boundhole)
{
  typedef Delaunay_remove_tds_3_2<Delaunay_triangulation_3> Surface;
  typedef Surface::Face_3_2 Face_3_2;
  typedef Surface::Vertex_3_2 Vertex_3_2;
  typedef Surface::Vertex_circulator Vertex_circulator_3_2;

  std::list<Cell_handle> cells;

  Surface surface(boundhole);

  int size = 5* boundhole.size();
  int opcount = 0;

  Face_3_2 *f = &(* surface.faces_begin());
  Face_3_2 *last_op;
  int k = -1;

  // This is a loop over the halfedges of the surface of the hole
  // we have a current face f, and look at its incident edges (k = 0,1,2)
  for(;;){
    k++;
    if(k == 3) {
      f = (Face_3_2*)f->n();
      if(f == last_op) {
	std::cerr << "\nUnable to find an ear\n" <<  surface << std::endl;
	delete_cells(cells);
	return false;
      }
      k = 0;
    }
    if(f->edge(k)) {
      bool violation = false;
      Vertex_3_2 *w0, *w1, *w2, *w3;
      Vertex *v0, *v1, *v2, *v3;
      int i = ccw(k);
      int j = cw(k);
      Face_3_2 *n = f->neighbor(k);
      int fi = n->index(f);

      w0 = f->vertex(k);
      w1 = f->vertex(i);
      w2 = f->vertex(j);
      w3 = n->vertex(fi);

      v0 = w0->info();
      v1 = w1->info();
      v2 = w2->info();
      v3 = w3->info();

      const Point & p0 = v0->point();
      const Point & p1 = v1->point();
      const Point & p2 = v2->point();
      const Point & p3 = v3->point();


      bool inf_0 = is_infinite(Vertex_handle(v0));
      bool inf_1 = is_infinite(Vertex_handle(v1));
      bool inf_2 = is_infinite(Vertex_handle(v2));
      bool inf_3 = is_infinite(Vertex_handle(v3));

      if( inf_1 || inf_2 ){
	// there will be another ear, so let's ignore this one,
	// because it is complicated to treat
      } else if( inf_0 || inf_3 ||
		 (orientation(p0, p1, p2, p3) == POSITIVE) ) {
	// the two faces form a concavity, in which we might plug a tetrahedron

	int cospheric = 0;
	std::set<Vertex_3_2*> cospheric_vertices;
	bool on_unbounded_side = false;
	// we now look at all vertices that are on the boundary of the hole
	for(Surface::Vertex_iterator vit = surface.vertices_begin();
	    (! violation ) && (vit != surface.vertices_end());
	    vit++) {
	  Vertex *v = (*vit).info();
	  if( (! is_infinite(Vertex_handle(v)))
	      && (v != v0) && (v != v1) && (v != v2) && (v != v3)) {
	    const Point & p = v->point();

	    Bounded_side bs; // = side_of_sphere(v0, v1, v2, v3, p);

	    if(inf_0) {
	      bs = side_of_sphere_inf(p2, p1, p3, p);
	    } else if(inf_3) {
	      bs = side_of_sphere_inf(p0, p1, p2, p);
	    } else {
	      bs = Bounded_side( side_of_oriented_sphere(p0, p1, p2, p3, p) );
	    }
	    //CGAL_triangulation_assertion(bs == side_of_sphere(v0, v1, v2, v3, p));

	    if((bs == ON_BOUNDARY)) {
	      cospheric++;
	      cospheric_vertices.insert(&(*vit));
	    }
	    violation = (bs == ON_BOUNDED_SIDE);
	    on_unbounded_side |= (bs == on_unbounded_side);
	  }
	}

	// if there are cospheric points we have to test a little bit more
	if( (! violation) && (cospheric > 0) ) {
	  if(inf_0 || inf_3) {
	    // the cospheric points are on the boundary of the convex hull we don't care
	    if(! on_unbounded_side) {
	      //std::cout << "all points are coplanar" << std::endl;
	    }
	  }  else  {
	    // for all edges that are incident to w2, check if the other vertex is cospheric
	    // and if the edge  is coplanar to plane(v0, v2, v3) 
	    Vertex_circulator_3_2 vc = w2->incident_vertices();
	    Vertex_circulator_3_2 done = vc;
	    std::set<Vertex_3_2*>::iterator co_it, not_found = cospheric_vertices.end();
	    do {
	      if( ! ((co_it = cospheric_vertices.find(&(*vc))) == not_found) ) {
		const Point & pc = (*co_it)->info()->point();
		violation = orientation(p0,p2,p3,pc) == COPLANAR;
	      }
	    } while( (! violation) && (++vc != done));

	    if( ! violation ) {
	      // for all edges that are incident to w1, check if the other vertex is cospheric
	      // and if the edge  is coplanar to plane(v0, v1, v3) 
	      Vertex_circulator_3_2 vc = w1->incident_vertices();
	      Vertex_circulator_3_2 done = vc;
	      std::set<Vertex_3_2*>::iterator co_it, not_found = cospheric_vertices.end();
	      do {
		if( ! ((co_it = cospheric_vertices.find(&(*vc))) == not_found) ) {
		  const Point & pc = (*co_it)->info()->point();
		  violation = (orientation(p0,p1,p3,pc) == COPLANAR);
		}
	      } while( (! violation) && (++vc != done));
	    }
	  }
	  if( (! violation) && (! inf_0) && (! inf_1) && (! inf_2) && (! inf_3) ) {
	    for(std::set<Vertex_3_2*>::iterator it = cospheric_vertices.begin();
		it != cospheric_vertices.end();
		it++) {
	      const Point & pit = (*it)->info()->point();
	      Vertex_circulator_3_2 vc = (*it)->incident_vertices();
	      Vertex_circulator_3_2 done = vc;
	      std::set<Vertex_3_2*>::iterator co_it, not_found = cospheric_vertices.end();
	      do {
		if( ! ((co_it = cospheric_vertices.find(&(*vc))) == not_found) ) {
		  const Point & pc = (*co_it)->info()->point();
		  violation = (orientation(p0,p3,pc, pit) == COPLANAR) 
		    && (coplanar_orientation(p0,p3,pc, pit) == NEGATIVE);
		}
	      } while( (! violation) && (++vc != done));
	    }
	  }
	} // test a little bit more

	Face_3_2 *m_i = f->neighbor(i);
	Face_3_2 *m_j = f->neighbor(j); 
	bool neighbor_i = m_i == n->neighbor(cw(fi));
	bool neighbor_j = m_j == n->neighbor(ccw(fi));

	if((! violation) && (! (((! neighbor_i) && (! neighbor_j)) 
				&& surface.is_edge(f->vertex(k), n->vertex(fi))))) {
	  // none of the vertices violates the Delaunay property
	  // and if we are in the flip case, the edge that would get introduced is not on the surface
	  // We are ready to plug the tetrahedron
	  // It may touch 2 triangles, 

	  Cell_handle ch = create_cell(Vertex_handle(v0), Vertex_handle(v1), Vertex_handle(v2), Vertex_handle(v3),
				       NULL, NULL, NULL, NULL);
	  cells.push_back(ch);
	  Cell* c = handle2pointer(ch);
	  v0->set_cell(c);
	  v1->set_cell(c);
	  v2->set_cell(c);
	  v3->set_cell(c);
	 
	  //print(ch);
	  /*
Removal of point 265 : -7 -18 -15

-9 -18 -14; -7 -21 -10; -8 -22 -16; -3 -21 -11;
flip
-1 -20 -15; -8 -22 -16; -3 -21 -11; -9 -18 -14;
flip
-9 -18 -14; -1 -20 -15; -8 -22 -16; -7 -19 -24;
remove
-3 -11 -20; -7 -19 -24; -1 -20 -15; -9 -18 -14;
flip
-9 -18 -14; -3 -11 -20; -7 -19 -24; -11 -14 -13;
remove
-9 -18 -14; -3 -11 -20; -11 -14 -13; -7 -9 -11;
flip
-11 -14 -13; -7 -9 -11; -5 -12 -9; -3 -11 -20;
flip
-5 -12 -9; -9 -18 -14; -11 -14 -13; -7 -9 -11;
flip
-5 -12 -9; -7 -21 -10; -9 -18 -14; -3 -21 -11;
remove
-5 -12 -9; -3 -21 -11; -9 -18 -14; -1 -20 -15;
remove
-3 -11 -20; -1 -20 -15; -5 -12 -9; -9 -18 -14;
remove
-3 -11 -20; -9 -18 -14; -5 -12 -9; -7 -9 -11;
remove

Unable to find an ear
4 4 2
-3 -11 -20; -5 -12 -9; -11 -14 -13; -7 -9 -11
	   */

	  Facet fac = n->info();
	  c->set_neighbor(0, fac.first);
	  fac.first->set_neighbor(fac.second, c);
	  fac = f->info();
	  c->set_neighbor(3, fac.first);
	  fac.first->set_neighbor(fac.second, c);
	  CGAL_triangulation_assertion(c->neighbor(0) != c->neighbor(3));
	  // 3, or even 4 if it is the last tetrahedron

	  if(neighbor_i) {
	    fac = m_i->info();
	    c->set_neighbor(1, fac.first);
	    CGAL_triangulation_assertion(c->neighbor(1) != c->neighbor(3));
	    CGAL_triangulation_assertion(c->neighbor(1) != c->neighbor(0));
	    fac.first->set_neighbor(fac.second, c);
	  }
	  if(neighbor_j) {
	    fac = m_j->info();
	    c->set_neighbor(2, fac.first);
	    CGAL_triangulation_assertion(c->neighbor(2) != c->neighbor(3));
	    CGAL_triangulation_assertion(c->neighbor(2) != c->neighbor(0));
	    if(neighbor_i) {
	    CGAL_triangulation_assertion(c->neighbor(2) != c->neighbor(1));
	    }
	    fac.first->set_neighbor(fac.second, c);
	  }

	  if((! neighbor_i) && (! neighbor_j)) {
	    surface.flip(f,k);
	    //std::cout << "flip" << std::endl;
	    int ni = f->index(n);
	    f->set_edge(ni, false);
	    f->set_edge(cw(ni));
	    f->set_edge(ccw(ni));
	    last_op = f; k = -1;

	    int fi = n->index(f);
	    n->set_edge(fi, false);
	    n->set_edge(cw(fi), f);
	    n->set_edge(ccw(fi), f);

	    f->set_info(std::make_pair(Cell_handle(c),2));
	    n->set_info(std::make_pair(Cell_handle(c),1));
	  } else if (neighbor_i && (! neighbor_j)) {
	    f->unlink(j);
	    surface.remove_degree_3(f->vertex(j), f);
	    //std::cout << "remove" << std::endl;
	    f->set_edge();
	    last_op = f; k = -1;
	    f->set_info(std::make_pair(Cell_handle(c),2));
	  } else if ((! neighbor_i) && neighbor_j)  {
	    f->unlink(i);
	    surface.remove_degree_3(f->vertex(i), f);
	    //std::cout << "remove" << std::endl;
	    f->set_edge();
	    last_op = f; k = -1;
	    f->set_info(std::make_pair(Cell_handle(c),1));
	  } else {
	    if(surface.number_of_vertices() != 4) {
	      delete_cells(cells);
	      return false;
	    } else {
	      // when we leave the function the vertices and faces of the surface
	      // are deleted by the destructor
	      return true;
	    }
	  }
	}// if(! violation)
      }// if(geom_traits()..
    }// if(f->edge(k))
  } // for(;;)
}



CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_3_H
