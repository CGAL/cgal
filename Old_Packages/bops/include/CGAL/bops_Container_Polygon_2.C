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
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_Container_Polygon_2.C
// package       : bops (2.2)
// source        : include/CGAL/bops_Container_Polygon_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_CONTAINER_POLYGON_2_C
#define CGAL_BOPS_CONTAINER_POLYGON_2_C


#include <CGAL/bops_simple_polygons_2.h>
#include <CGAL/bops_Convex_Polygon_2.h>
#include <CGAL/bops_Container_Polygon_2.h>
#include <CGAL/bops_assertions.h>

CGAL_BEGIN_NAMESPACE

/*
   transforms a sequence of special objects
   into a sequence of type Object
*/
template <class ForwardIterator, class OutputIterator>
OutputIterator _transform_to_object(
         ForwardIterator first, ForwardIterator last, OutputIterator result) {
    while (first != last) *result++ = make_object(*first++);
    return result;
}




/***********************************/ 
/* For Polygons with traits class: */
/***********************************/ 

template < class ForwardIterator, class Traits > 
bool do_intersect(ForwardIterator Afirst, ForwardIterator Alast,
		       ForwardIterator Bfirst, ForwardIterator Blast,
		       Traits &)
{
  typename Traits::Input_polygon A(Afirst,Alast);
  typename Traits::Input_polygon B(Bfirst,Blast);
  
  Bops_Simple_Polygons_2_Intersection<Traits> intersection(A,B);
  return intersection.do_intersect();
}


template < class ForwardIterator, class OutputIterator, class Traits > 
OutputIterator intersection(
	ForwardIterator Afirst,    // first "polygon", defined by a sequence of points
	ForwardIterator Alast,
	ForwardIterator Bfirst,    // second "polygon"
	ForwardIterator Blast,
	Traits &,                     // traits class
	OutputIterator result    // the result
)
{
  if( Afirst == Alast || Bfirst == Blast )  // if one polygon is empty
    return result;                          // then the intersection is empty

  typename Traits::Input_polygon A(Afirst,Alast);
  typename Traits::Input_polygon B(Bfirst,Blast);
  

  /*** NOT IMPLEMENTED YET
  if( A.size() == 3 && B.size() == 3 ) 
     return intersection(
	Triangle_2<Traits::R>(*Afirst++,*Afirst++,*Afirst),
        Triangle_2<Traits::R>(*Bfirst++,*Bfirst++,*Bfirst),
        result);
  */
  
	Bops_Polygons_2<Traits> *bops;
  
	if( A.is_convex() && B.is_convex() )
		bops= new Bops_Convex_Polygons_2_Intersection<Traits>(A,B);
	else
		bops= new Bops_Simple_Polygons_2_Intersection<Traits>(A,B);
		
  	bops->operation();
	std::copy(bops->begin(), bops->end(), result);
	
	return result;
}



template < class ForwardIterator, class OutputIterator, class Traits > 
OutputIterator Union(
	ForwardIterator Afirst,    // first "polygon", defined by a sequence of points
	ForwardIterator Alast,
	ForwardIterator Bfirst,    // second "polygon"
	ForwardIterator Blast,
	Traits &,                     // traits class
	OutputIterator result    // the result
)
{
  if( Afirst == Alast ) // if polygon A is empty then return polygon B
    return _transform_to_object(Bfirst, Blast, result);

  if( Bfirst == Blast ) // if polygon B is empty then return polygon A
    return _transform_to_object(Afirst, Alast, result);

  typename Traits::Input_polygon A(Afirst,Alast);
  typename Traits::Input_polygon B(Bfirst,Blast);
  Bops_Simple_Polygons_2_Union<Traits> bops(A,B);
  bops.operation();
  return std::copy(bops.begin(), bops.end(), result);
}

template < class ForwardIterator, class OutputIterator, class Traits > 
OutputIterator difference(
	ForwardIterator Afirst,    // first "polygon", defined by a sequence of points
	ForwardIterator Alast,
	ForwardIterator Bfirst,    // second "polygon"
	ForwardIterator Blast,
	Traits &,                     // traits class
	OutputIterator result    // the result
)
{
  if( Afirst == Alast ) // if polygon A is empty then the difference is empty
    return result;      

  if( Bfirst == Blast ) // if polygon B is empty then return polygon A
    return _transform_to_object(Afirst, Alast, result);

  typename Traits::Input_polygon A(Afirst,Alast);
  typename Traits::Input_polygon B(Bfirst,Blast);
  Bops_Simple_Polygons_2_Difference<Traits> bops(A,B);
  bops.operation();
  return std::copy(bops.begin(), bops.end(), result);
}

CGAL_END_NAMESPACE

#endif // CGAL_BOPS_CONTAINER_POLYGON_2_C
