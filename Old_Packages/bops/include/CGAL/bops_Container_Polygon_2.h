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
// file          : include/CGAL/bops_Container_Polygon_2.h
// package       : bops (2.2)
// source        : include/CGAL/bops_Container_Polygon_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_CONTAINER_POLYGON_2_H
#define CGAL_BOPS_CONTAINER_POLYGON_2_H

CGAL_BEGIN_NAMESPACE

/***********************************/ 
/* For Polygons with traits class: */
/***********************************/ 

template < class ForwardIterator, class Traits >
bool do_intersect( ForwardIterator Afirst, ForwardIterator Alast,
		        ForwardIterator Bfirst, ForwardIterator Blast,
		        Traits &);


template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
intersection( ForwardIterator Afirst, ForwardIterator Alast,
     	           ForwardIterator Bfirst, ForwardIterator Blast,
                   Traits &, OutputIterator list_of_objects_it);


template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
Union(        ForwardIterator Afirst, ForwardIterator Alast,
                   ForwardIterator Bfirst, ForwardIterator Blast,
		   Traits &, OutputIterator list_of_objects_it);


template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
difference(   ForwardIterator Afirst, ForwardIterator Alast,
		   ForwardIterator Bfirst, ForwardIterator Blast,
		   Traits &, OutputIterator list_of_objects_it);

CGAL_END_NAMESPACE 
 
#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_Container_Polygon_2.C>
#endif

#endif // CGAL_BOPS_CONTAINER_POLYGON_2_H
