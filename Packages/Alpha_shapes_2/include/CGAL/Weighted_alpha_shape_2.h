// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/Weighted_alpha_shape_2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_2_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_2_H

#ifdef CGAL_STL_SGI_CC
#define STL_HASH_TABLES
#endif

#include <CGAL/basic.h>

// SGI MIPSpro C++ 7.x
#ifdef STL_HASH_TABLES
#include <hash_set>
#include <hash_map>
#include <map>
#else
#include <set>
#include <map>
#endif

#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/less_partial.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Rt >
class Weighted_alpha_shape_2 : public Alpha_shape_2<Rt>
{

  //------------------------- TYPES ------------------------------------

public:
  typedef typename Alpha_shape_2<Rt>::Gt Gt;
  typedef typename Alpha_shape_2<Rt>::Triangulation_data_structure Tds;

  typedef typename Gt::Coord_type Coord_type;
  typedef typename Gt::Point Point;

  typedef typename Gt::Distance Distance;
  typedef typename Gt::Ray Ray;
  typedef typename Gt::Line Line;
  typedef typename Gt::Direction Direction;


  typedef typename Alpha_shape_2<Rt>::Face_handle Face_handle;
  typedef typename Alpha_shape_2<Rt>::Vertex_handle Vertex_handle;
  typedef typename Alpha_shape_2<Rt>::Edge Edge;

  typedef typename Alpha_shape_2<Rt>::Face_circulator Face_circulator;
  typedef typename Alpha_shape_2<Rt>::Edge_circulator Edge_circulator;
  typedef typename Alpha_shape_2<Rt>::Vertex_circulator Vertex_circulator;

  typedef typename Alpha_shape_2<Rt>::Face_iterator Face_iterator;
  typedef typename Alpha_shape_2<Rt>::Edge_iterator Edge_iterator;
  typedef typename Alpha_shape_2<Rt>::Vertex_iterator Vertex_iterator;

  typedef typename Alpha_shape_2<Rt>::Locate_type Locate_type;
  
  typedef typename Alpha_shape_2<Rt>::Mode Mode;
  
public:

  //------------------------- CONSTRUCTORS ------------------------------
 
  // Introduces an empty alpha-shape `A' for a positive
  // alpha-value `alpha'. Precondition: `alpha' >= 0.
  Weighted_alpha_shape_2(Coord_type alpha = 0, 
			 Mode m = GENERAL)
    : Alpha_shape_2<Rt>(alpha, m)
    {}
 
  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES

  template <class InputIterator>
  Weighted_alpha_shape_2( InputIterator first,  
			  InputIterator last,  
			  const Coord_type& alpha = 0,
			  Mode = GENERAL)
    : Alpha_shape_2<Rt>(first, last, alpha, m) 
    {}

#else //CGAL_CFG_NO_MEMBER_TEMPLATES
#if defined(LIST) || defined(__SGI_STL_LIST)

  Weighted_alpha_shape_2(typename std::list<Point>::const_iterator first,
			 typename std::list<Point>::const_iterator last,
			 const Coord_type& alpha = 0,
			 Mode m = GENERAL) 
    : Alpha_shape_2<Rt>(first, last, alpha, m) 
    {}

#endif //LIST

#if defined(VECTOR) || defined(__SGI_STL_VECTOR)

  Weighted_alpha_shape_2(typename std::vector<Point>::const_iterator first,
			 typename std::vector<Point>::const_iterator last,
			 const Coord_type& alpha = 0,
			 Mode m = GENERAL) 
    : Alpha_shape_2<Rt>(first, last, alpha, m) 
    {}

#endif // VECTOR

#endif //CGAL_CFG_NO_MEMBER_TEMPLATES 


  //---------------heuristic initialization of weights----------------

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES

  template <class Iterator>
  void
  initialize_weights_to_the_nearest_voronoi_vertex
  (Iterator  first, Iterator  last, const Coord_type &k);

  //---------------------------------------------------------------------

  template <class Iterator>
  void
  initialize_weights_to_the_nearest_voronoi_edge
  (Iterator  first, Iterator  last, const Coord_type &k);

  //---------------------------------------------------------------------

  template <class Iterator>
  void
  initialize_weights_to_the_nearest_vertex
  (Iterator  first, Iterator  last, const Coord_type &k);

  
#else //CGAL_CFG_NO_MEMBER_TEMPLATES
#if defined(LIST) || defined(__SGI_STL_LIST)

  void
  initialize_weights_to_the_nearest_voronoi_vertex
  (typename std::list<Point>::iterator  first,
   typename std::list<Point>::iterator  last, const Coord_type &k);

  //---------------------------------------------------------------------

  void
  initialize_weights_to_the_nearest_voronoi_edge
  (typename std::list<Point>::iterator  first,
   typename std::list<Point>::iterator  last, const Coord_type &k);


  //---------------------------------------------------------------------

  void
  initialize_weights_to_the_nearest_vertex
  (typename std::list<Point>::iterator  first,
   typename std::list<Point>::iterator  last, const Coord_type &k);
  

#endif //LIST

#if defined(VECTOR) || defined(__SGI_STL_VECTOR)

  void
  initialize_weights_to_the_nearest_voronoi_vertex
  (typename std::vector<Point>::iterator  first,
   typename std::vector<Point>::iterator  last, const Coord_type &k);

  //---------------------------------------------------------------------
  
  void
  initialize_weights_to_the_nearest_voronoi_edge
  (typename std::vector<Point>::iterator  first,
   typename std::vector<Point>::iterator  last, const Coord_type &k);

  //---------------------------------------------------------------------

  void
  initialize_weights_to_the_nearest_vertex
  (typename std::vector<Point>::iterator  first,
   typename std::vector<Point>::iterator  last, const Coord_type &k);
  
  //---------------------------------------------------------------------
   

#endif // VECTOR
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
};

//---------------------- MEMBER FUNCTIONS -----------------------------

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES

template <class Rt> template <class Iterator>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_vertex
(Iterator  first, Iterator  last, const Coord_type &k)
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  Iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      Face_circulator face_circ=D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;
      double d=DBL_MAX;
      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  Point p = D.dual(f);
		  double dd = squared_distance(p, *point_it);
		  d = std::min(dd, d);
		}
	    }
	  while(++face_circ != done);
	}
      (*point_it) = Point((*point_it).point(), k*k*d);
    }
}


//---------------------------------------------------------------------

template <class Rt> template <class Iterator>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_edge
(Iterator  first, Iterator  last, const Coord_type &k)
{ 
  typedef  Delaunay_triangulation_2<Gt, Tds> Dt_int;
  Dt_int D;	
  Iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      typename Dt_int::Face_circulator face_circ=
	D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;

      double d = DBL_MAX;
      double dd = DBL_MAX;

      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      typename Dt_int::Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  for ( int i=0; i!=3; ++i)		    
		    {
		      typename Dt_int::Edge e(f,i);
		      
		      if ((!D.is_infinite(e.first))&&
			  (!D.is_infinite(e.first->neighbor(e.second))))
			{
			  typename Gt::Segment seg(D.dual(e.first),
						   D.dual(e.first
							  ->neighbor(e.second)));
			  dd = squared_distance(seg, *point_it);
			}
		      d = std::min(dd, d);
		    }
		}
	    }
	  while(++face_circ != done);
	}
      if (d != DBL_MAX)	
	(*point_it) = Point((*point_it).point(), k*k*d); 
      else	
	(*point_it) = Point((*point_it).point(), Coord_type(0));
    }
}


//---------------------------------------------------------------------

template <class Rt> template <class Iterator>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_vertex
(Iterator  first, Iterator  last, const Coord_type &k)
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  Iterator point_it;
    
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it) 
    { 

      D.remove(D.nearest_vertex(*point_it));
      
      Point neighbor = D.nearest_vertex(*point_it)->point();

      (*point_it) = Point((*point_it).point(),k*k* 
			  squared_distance(neighbor, *point_it));
			
      D.insert(*point_it);
    }
}
 
#else // CGAL_CFG_NO_MEMBER_TEMPLATES

//-------------------------------------------------------------------

#if defined(VECTOR) || defined(__SGI_STL_VECTOR)

template<class Rt> 
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_vertex
(typename std::vector<Point>::iterator  first,
 typename std::vector<Point>::iterator  last, const Coord_type &k)
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  typename std::vector<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      Face_circulator face_circ=D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;
      double d=DBL_MAX;
      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  Point p = D.dual(f);
		  double dd = squared_distance(p, *point_it);
		  d = std::min(dd, d);
		}
	    }
	  while(++face_circ != done);
	}
      (*point_it) = Point((*point_it).point(), k*k*d);
    }
}

//-----------------------------------------------------------------------
  
template<class Rt>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_edge
(typename std::vector<Point>::iterator  first,
 typename std::vector<Point>::iterator  last, const Coord_type &k)
{ 
  typedef  Delaunay_triangulation_2<Gt, Tds> Dt_int;
  Dt_int D;	
  typename std::vector<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      typename Dt_int::Face_circulator face_circ=
	D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;

      double d = DBL_MAX;
      double dd = DBL_MAX;

      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      typename Dt_int::Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  for ( int i=0; i!=3; ++i)		    
		    {
		      typename Dt_int::Edge e(f,i);
		      
		      if ((!D.is_infinite(e.first))&&
			  (!D.is_infinite(e.first->neighbor(e.second))))
			{
			  typename Gt::Segment seg(D.dual(e.first),
						   D.dual(e.first
							  ->neighbor(e.second)));
			  dd = squared_distance(seg, *point_it);
			}
		      d = std::min(dd, d);
		    }
		}
	    }
	  while(++face_circ != done);
	}
      if (d != DBL_MAX)	
	(*point_it) = Point((*point_it).point(), k*k*d); 
      else	
	(*point_it) = Point((*point_it).point(), Coord_type(0));
    }
}

//-----------------------------------------------------------------------

template<class Rt>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_vertex
(typename std::vector<Point>::iterator  first,
 typename std::vector<Point>::iterator  last, const Coord_type &k) 
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  typename std::vector<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it) 
    { 

      D.remove(D.nearest_vertex(*point_it));
      
      Point neighbor = D.nearest_vertex(*point_it)->point();

      (*point_it) = Point((*point_it).point(),k*k* 
			  squared_distance(neighbor, *point_it));
      D.insert(*point_it);
    }
}

#endif //VECTOR

//-----------------------------------------------------------------------

#if defined(LIST) || defined(__SGI_STL_LIST)

template<class Rt> 
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_vertex
(typename std::list<Point>::iterator  first,
 typename std::list<Point>::iterator  last, const Coord_type &k)
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  typename std::list<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      Face_circulator face_circ=D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;
      double d=DBL_MAX;
      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  Point p = D.dual(f);
		  double dd = squared_distance(p, *point_it);
		  d = std::min(dd, d);
		}
	    }
	  while(++face_circ != done);
	}
      (*point_it) = Point((*point_it).point(), k*k*d);
    }
}

//-----------------------------------------------------------------------
  
template<class Rt>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_voronoi_edge
(typename std::list<Point>::iterator  first,
 typename std::list<Point>::iterator  last, const Coord_type &k)
{ 
  typedef  Delaunay_triangulation_2<Gt, Tds> Dt_int;
  Dt_int D;	
  typename std::list<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it)    
    { 
      typename Dt_int::Face_circulator face_circ=
	D.incident_faces(D.nearest_vertex(*point_it)),
	done = face_circ;

      double d = DBL_MAX;
      double dd = DBL_MAX;

      if (!face_circ.is_empty())	
	{
	  do	    
	    {
	      typename Dt_int::Face_handle f = face_circ;
	      if (!D.is_infinite(f))		
		{
		  for ( int i=0; i!=3; ++i)		    
		    {
		      typename Dt_int::Edge e(f,i);
		      
		      if ((!D.is_infinite(e.first))&&
			  (!D.is_infinite(e.first->neighbor(e.second))))
			{
			  typename Gt::Segment seg(D.dual(e.first),
						   D.dual(e.first
							  ->neighbor(e.second)));
			  dd = squared_distance(seg, *point_it);
			}
		      d = std::min(dd, d);
		    }
		}
	    }
	  while(++face_circ != done);
	}
      if (d != DBL_MAX)	
	(*point_it) = Point((*point_it).point(), k*k*d); 
      else	
	(*point_it) = Point((*point_it).point(), Coord_type(0));
    }
}

//-----------------------------------------------------------------------

template<class Rt>
void
Weighted_alpha_shape_2<Rt>::initialize_weights_to_the_nearest_vertex
(typename std::list<Point>::iterator  first,
 typename std::list<Point>::iterator  last, const Coord_type &k) 
{ 
  Delaunay_triangulation_2<Gt, Tds> D;
  typename std::list<Point>::iterator point_it;
  
  D.insert(first, last);

  for( point_it = first; 
       point_it != last; 
       ++point_it) 
    { 

      D.remove(D.nearest_vertex(*point_it));
      
      Point neighbor = D.nearest_vertex(*point_it)->point();

      (*point_it) = Point((*point_it).point(),k*k* 
			  squared_distance(neighbor, *point_it));
      D.insert(*point_it);
    }
}

//-------------------------------------------------------------------

#endif // LIST
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_WEIGHTED_ALPHA_2_H
