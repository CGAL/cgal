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
// file          : include/CGAL/nsquare_intersecting.h
// package       : bops (2.2)
// source        : include/CGAL/nsquare_intersecting.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Carl Van Geem <Carl.Van.Geem@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

/*
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   An n^2 algorithm for the computation of the intersections 
    of segments belonging to two differents lists. It's the 
    first step in the computation of Boolean Operations on 
    simple polygons in 2D.

   Last changes: 26.Mar.1997
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
#ifndef CGAL_NSQUARE_INTERSECTING_H
#define CGAL_NSQUARE_INTERSECTING_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#include <list>

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif
#ifndef CGAL_SEGMENT_2_H
#include <CGAL/Segment_2.h>
#endif
#ifndef CGAL_OBJECT_H
#include <CGAL/Object.h>
#endif

#ifndef CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H
#include <CGAL/Segment_2_Segment_2_intersection.h>
#endif
#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif

#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class R >
class Intersectionresult
{
  typedef Point_2<R> Point;
  typedef Segment_2<R> Segment;
  typedef std::list<Segment> Segment_list;
  //typedef Polygon_2<Polygon_traits_2<R>, Container > Polygon;
  //typedef typename Polygon::Edge_const_iterator edge_iterator;

 private:
  Object                    _intersection_object;
  Segment_list      _segments_poly1;
  Segment_list      _segments_poly2;
  bool                           _is_vertex_of_poly1;
  bool                           _is_vertex_of_poly2;
  bool                           _is_edge_of_poly1;
  bool                           _is_edge_of_poly2;
  bool                           _is_vertex_vertex_intersection;
  bool                           _is_vertex_edge_intersection;
  bool                           _is_edge_edge_intersection;
 public:
  Intersectionresult() ;
  Intersectionresult (const Intersectionresult<R> &ires ) ;
  Intersectionresult (const Point &iobj,
                           const Segment &iseg1,
                           const Segment &iseg2,
                           const bool &i_is_vertex_of_poly1,
                           const bool &i_is_vertex_of_poly2,
                           const bool &i_is_edge_of_poly1,
                           const bool &i_is_edge_of_poly2,
                           const bool &i_is_vertex_vertex_intersection,
                           const bool &i_is_vertex_edge_intersection,
                           const bool &i_is_edge_edge_intersection) ;
  Intersectionresult (const Segment &iobj,
                           const Segment &iseg1,
                           const Segment &iseg2,
                           const bool &i_is_vertex_of_poly1,
                           const bool &i_is_vertex_of_poly2,
                           const bool &i_is_edge_of_poly1,
                           const bool &i_is_edge_of_poly2,
                           const bool &i_is_vertex_vertex_intersection,
                           const bool &i_is_vertex_edge_intersection,
                           const bool &i_is_edge_edge_intersection) ;
  Intersectionresult<R>  &operator=(const Intersectionresult<R> 
                                                                   &ires) ; 
  ~Intersectionresult() ;


  Object     intersection_object() const { return _intersection_object;}
  Object&    intersection_object() { return _intersection_object;}
  Segment_list  segments_poly1() const { return _segments_poly1;}
  Segment_list& segments_poly1() { return _segments_poly1;}
  Segment_list  segments_poly2() const { return _segments_poly2;}
  Segment_list& segments_poly2() { return _segments_poly2;}
  bool   is_vertex_of_poly1() const { return _is_vertex_of_poly1;}
  bool&  is_vertex_of_poly1() { return _is_vertex_of_poly1;}
  bool   is_vertex_of_poly2() const { return _is_vertex_of_poly2;}
  bool&  is_vertex_of_poly2() { return _is_vertex_of_poly2;}
  bool   is_edge_of_poly1()   const { return _is_edge_of_poly1;}
  bool&  is_edge_of_poly1()   { return _is_edge_of_poly1;}
  bool   is_edge_of_poly2()   const { return _is_edge_of_poly2;}
  bool&  is_edge_of_poly2()   { return _is_edge_of_poly2;}
  bool   is_vertex_vertex_intersection() const
           { return _is_vertex_vertex_intersection;}
  bool&  is_vertex_vertex_intersection() 
           { return _is_vertex_vertex_intersection;}
  bool   is_vertex_edge_intersection()   const  
           { return _is_vertex_edge_intersection;}
  bool&  is_vertex_edge_intersection()   
           { return _is_vertex_edge_intersection;}
  bool   is_edge_edge_intersection()     const  
           { return _is_edge_edge_intersection;}
  bool&  is_edge_edge_intersection()     
           { return _is_edge_edge_intersection;}


};






template < class R > 
class List_of_intersections
{
 public:
  typedef Point_2<R> Point;
  typedef Segment_2<R> Segment;
  typedef std::list<Segment> Segment_list;
  //typedef Polygon_2<Polygon_traits_2<R>, Container > Polygon;
  //typedef typename Polygon::Edge_const_iterator edge_iterator;

 private:
   std::list<Intersectionresult<R> >    _list_of_intersections;
 public:
   List_of_intersections() ;
   List_of_intersections(const List_of_intersections<R> &ilist) ;
   List_of_intersections<R> &operator=(
                              const List_of_intersections<R> &ilist) ;
   ~List_of_intersections() ;


   std::list<Intersectionresult<R> >  list_of_intersections() const {
                                   return _list_of_intersections; }
   std::list<Intersectionresult<R> >&  list_of_intersections() {
                                   return _list_of_intersections; }
   void add( const Point &iobj,
             const Segment &iseg1,
             const Segment &iseg2,
             const bool &i_is_vertex_of_poly1,
             const bool &i_is_vertex_of_poly2,
             const bool &i_is_edge_of_poly1,
             const bool &i_is_edge_of_poly2,
             const bool &i_is_vertex_vertex_intersection,
             const bool &i_is_vertex_edge_intersection,
             const bool &i_is_edge_edge_intersection) 
      {
       Intersectionresult<R> ires( 
            iobj,
            iseg1,
            iseg2,
            i_is_vertex_of_poly1,
            i_is_vertex_of_poly2,
            i_is_edge_of_poly1,
            i_is_edge_of_poly2,
            i_is_vertex_vertex_intersection,
            i_is_vertex_edge_intersection,
            i_is_edge_edge_intersection);

       _list_of_intersections.push_back(ires); 
      }


   void add( const Segment &iobj,
             const Segment &iseg1,
             const Segment &iseg2,
             const bool &i_is_vertex_of_poly1,
             const bool &i_is_vertex_of_poly2,
             const bool &i_is_edge_of_poly1,
             const bool &i_is_edge_of_poly2,
             const bool &i_is_vertex_vertex_intersection,
             const bool &i_is_vertex_edge_intersection,
             const bool &i_is_edge_edge_intersection) 
      {
       Intersectionresult<R> ires( 
            iobj,
            iseg1,
            iseg2,
            i_is_vertex_of_poly1,
            i_is_vertex_of_poly2,
            i_is_edge_of_poly1,
            i_is_edge_of_poly2,
            i_is_vertex_vertex_intersection,
            i_is_vertex_edge_intersection,
            i_is_edge_edge_intersection);

       _list_of_intersections.push_back(ires); 
      }

 };



template < class R, class Container >
class nsquareintersection {
 public:
  typedef Point_2<R> Point;
  typedef Segment_2<R> Segment;
  typedef std::list<Segment> Segment_list;
  typedef Polygon_2<Polygon_traits_2<R>, Container > Polygon;
  typedef typename Polygon::Edge_const_iterator edge_iterator;
  typedef typename Segment_list::iterator segment_iterator;

  nsquareintersection() {}
  std::list<Intersectionresult<R> > operator()(
                          segment_iterator beginpoly1,
                          segment_iterator endpoly1,
                          segment_iterator beginpoly2,
                          segment_iterator endpoly2) 
 {
  List_of_intersections<R> tobereturned;
  std::list<Intersectionresult<R> >::iterator it;
  segment_iterator it1;
  segment_iterator it2;
  Object    result;
  Point ipoint;

  Segment isegment;
  Point acgalpoint;
  bool notfound;
  bool i_is_vertex_of_poly1;
  bool i_is_vertex_of_poly2;
  bool i_is_edge_of_poly1;
  bool i_is_edge_of_poly2;
  bool i_is_vertex_vertex_intersection;
  bool i_is_vertex_edge_intersection;
  bool i_is_edge_edge_intersection;

  for(it1 = beginpoly1; it1 != endpoly1; it1++ )
   {
    for(it2 = beginpoly2; it2 != endpoly2; it2++ )
     {
      if ( do_intersect(*it1, *it2) )
       {
        result = intersection(*it1, *it2);
        if ( assign(ipoint, result) )
         {/* report intersection */
          it = tobereturned.list_of_intersections().begin();
          notfound = true;
          while ((it != tobereturned.list_of_intersections().end())
                  && (notfound) )
            {
             if ( ( assign(acgalpoint, (*it).intersection_object()) )
              &&  ( acgalpoint == ipoint ) )
               {
                notfound = false;
                /* adapt the information in *it */
                (*it).segments_poly1().push_back(*it1);
                (*it).segments_poly2().push_back(*it2);
/*              (*it).is_edge_of_poly1() = false;
                (*it).is_edge_of_poly2() = false;

                (*it).is_edge_edge_intersection() = false; */

                if ((ipoint == (*it1).min())||(ipoint == (*it1).max()))
                  {  (*it).is_vertex_of_poly1() = true; }
                if ((ipoint == (*it2).min())||(ipoint == (*it2).max()))
                  {  (*it).is_vertex_of_poly2() = true; }

                if (  ((*it).is_vertex_of_poly1())
                    &&((*it).is_vertex_of_poly2()))
                  {  (*it).is_vertex_vertex_intersection() = true; }
                 else
                  {
                   (*it).is_vertex_vertex_intersection() = false;
                   if ( ((*it).is_vertex_of_poly1())
                     || ((*it).is_vertex_of_poly2()))
                     { (*it).is_vertex_edge_intersection() = true; 
		       if ( ( (*it).is_vertex_of_poly1() ) &&
			    ( (*it).segments_poly2().size() > 1 ) )
			 (*it).segments_poly2().pop_back();
		       if ( ( (*it).is_vertex_of_poly2() ) &&
			    ( (*it).segments_poly1().size() > 1 ) )
			 (*it).segments_poly1().pop_back();
		     }
                    else
                     { (*it).is_vertex_edge_intersection() = false;}
                  }
               }
              it++;
            }
          if (notfound)

           {/* add this intersection point to the list */
            i_is_vertex_of_poly1 =
               ((ipoint == (*it1).min())||(ipoint == (*it1).max()));
            i_is_vertex_of_poly2 =
               ((ipoint == (*it2).min())||(ipoint == (*it2).max()));
            i_is_vertex_vertex_intersection =
               ( i_is_vertex_of_poly1 && i_is_vertex_of_poly2 );
            i_is_vertex_edge_intersection =
               ( (( i_is_vertex_of_poly1)&&(!i_is_vertex_of_poly2)) ||
                 ((!i_is_vertex_of_poly1)&&( i_is_vertex_of_poly2)) ) ;
            i_is_edge_edge_intersection = !i_is_vertex_edge_intersection;

            tobereturned.add(ipoint,*it1,*it2, i_is_vertex_of_poly1,
                         i_is_vertex_of_poly2,false,false,
                         i_is_vertex_vertex_intersection,
                         i_is_vertex_edge_intersection,
                         i_is_edge_edge_intersection) ;
           }
         }
        else if ( assign(isegment, result) )
         {/* report intersection */
          i_is_edge_of_poly1 = (isegment == *it1) ;
          i_is_edge_of_poly2 = (isegment == *it2) ;

          tobereturned.add(isegment,*it1,*it2,false,false,
                                         i_is_edge_of_poly1,
                                         i_is_edge_of_poly2,
                                         false,false,true );
         }
       }
     }
   }
  return tobereturned.list_of_intersections();
 }

  std::list<Intersectionresult<R> > operator()(
         edge_iterator beginpoly1,
         edge_iterator endpoly1,
         edge_iterator beginpoly2,
         edge_iterator endpoly2)
 {
  List_of_intersections<R> tobereturned;
  std::list<Intersectionresult<R> >::iterator it;
  edge_iterator it1;
  edge_iterator it2;
  Object    result;
  Point ipoint;

  Segment isegment;
  Point acgalpoint;
  bool notfound;
  bool i_is_vertex_of_poly1;
  bool i_is_vertex_of_poly2;
  bool i_is_edge_of_poly1;
  bool i_is_edge_of_poly2;
  bool i_is_vertex_vertex_intersection;
  bool i_is_vertex_edge_intersection;
  bool i_is_edge_edge_intersection;
 
  for(it1 = beginpoly1; it1 != endpoly1; ++it1 )
   {
    for(it2 = beginpoly2; it2 != endpoly2; ++it2 )
     {
      if ( do_intersect(*it1, *it2) )
       {
        result = intersection(*it1, *it2);
        if ( assign(ipoint, result) )
         {/* report intersection */
          it = tobereturned.list_of_intersections().begin();
          notfound = true;
          while ((it != tobereturned.list_of_intersections().end())
                  && (notfound) )
            {
             if ( ( assign(acgalpoint, (*it).intersection_object()) )
              &&  ( acgalpoint == ipoint ) )
               {
                notfound = false;
                /* adapt the information in *it */
                (*it).segments_poly1().push_back(*it1);
                (*it).segments_poly2().push_back(*it2);
/*              (*it).is_edge_of_poly1() = false;
                (*it).is_edge_of_poly2() = false;

                (*it).is_edge_edge_intersection() = false; */

                if ((ipoint == (*it1).min())||(ipoint == (*it1).max()))
                  {  (*it).is_vertex_of_poly1() = true; }
                if ((ipoint == (*it2).min())||(ipoint == (*it2).max()))
                  {  (*it).is_vertex_of_poly2() = true; }

                if (  ((*it).is_vertex_of_poly1())
                    &&((*it).is_vertex_of_poly2()))
                  {  (*it).is_vertex_vertex_intersection() = true; }
                 else
                  {
                   (*it).is_vertex_vertex_intersection() = false;
                   if ( ((*it).is_vertex_of_poly1())
                     || ((*it).is_vertex_of_poly2()))
                     { (*it).is_vertex_edge_intersection() = true; 
		       if ( ( (*it).is_vertex_of_poly1() ) &&
			    ( (*it).segments_poly2().size() > 1 ) )
			 (*it).segments_poly2().pop_back();
		       if ( ( (*it).is_vertex_of_poly2() ) &&
			    ( (*it).segments_poly1().size() > 1 ) )
			 (*it).segments_poly1().pop_back();
		     }
                    else
                     { (*it).is_vertex_edge_intersection() = false;}
                  }
               }
              it++;
            }
          if (notfound)

           {/* add this intersection point to the list */
            i_is_vertex_of_poly1 =
               ((ipoint == (*it1).min())||(ipoint == (*it1).max()));
            i_is_vertex_of_poly2 =
               ((ipoint == (*it2).min())||(ipoint == (*it2).max()));
            i_is_vertex_vertex_intersection =
               ( i_is_vertex_of_poly1 && i_is_vertex_of_poly2 );
            i_is_vertex_edge_intersection =
               ( (( i_is_vertex_of_poly1)&&(!i_is_vertex_of_poly2)) ||
                 ((!i_is_vertex_of_poly1)&&( i_is_vertex_of_poly2)) ) ;
            i_is_edge_edge_intersection = !i_is_vertex_edge_intersection;

            tobereturned.add(ipoint,*it1,*it2, i_is_vertex_of_poly1,
                         i_is_vertex_of_poly2,false,false,
                         i_is_vertex_vertex_intersection,
                         i_is_vertex_edge_intersection,
                         i_is_edge_edge_intersection) ;
           }
         }
        else if ( assign(isegment, result) )
         {/* report intersection */
          i_is_edge_of_poly1 = (isegment == *it1) ;
          i_is_edge_of_poly2 = (isegment == *it2) ;

          tobereturned.add(isegment,*it1,*it2,false,false,
                                         i_is_edge_of_poly1,
                                         i_is_edge_of_poly2,
                                         false,false,true );
         }
       }
     }
   }
  return tobereturned.list_of_intersections();
 }


};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/nsquare_intersecting.C>
#endif

#endif //  CGAL_NSQUARE_INTERSECTING_H
