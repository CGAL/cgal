// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H
#define CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H 1

#include <CGAL/Iterators_for_face_graph_wrapper.h>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

////////////////////////////////////////////////////////////////////////////////
/** This file contains the following functors for Face_graph_wrapper:
 * Is_free<typename HEG, unsigned int i>::
         operator() (const HEG& heg, Dart_const_handle dh)
 * Get_beta<typename HEG, unsigned int i>::
         operator() (const HEG& heg, Dart_const_handle dh)
 // The following functors are not needed: we can use directly the ones in CMap
 // * Belong_to_same_cell<unsigned int i>::    // TODO add second parameter unsigned int d
 //         operator() (const HEG& heg, Dart_const_handle dh1, Dart_const_handle dh2)
 // * Mark_cell<unsigned int i>::              // TODO add second parameter unsigned int d
 //         operator() (const HEG& heg, Dart_const_handle dh)
 // * Unmark_cell<unsigned int i>::            // TODO add second parameter unsigned int d
 //         operator() (const HEG& heg, Dart_const_handle dh)
*/
////////////////////////////////////////////////////////////////////////////////
/// Is_free
template<typename HEG, unsigned int i>
class Is_free
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

  static bool value(const HEG& /*heg*/, Dart_const_handle /*dh*/) 
  { CGAL_static_assertion(i==0 || i==1); return false; }
};
template<typename HEG>
class Is_free<HEG, 2>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static bool value(const HEG& heg, Dart_const_handle dh) 
  { return CGAL::is_border(CGAL::opposite(dh, heg), heg); }
};
////////////////////////////////////////////////////////////////////////////////
/// Get_beta
template<typename HEG, unsigned int i>
class Get_beta
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

  static Dart_const_handle value(const HEG& /*heg*/, Dart_const_handle /*dh*/) 
  { /* CGAL_static_assertion(false);*/
    std::cout<<"ERROR Get_beta<HEG, "<<i<<">"<<std::endl;
    CGAL_assertion(false);
    return nullptr;
  }
};
template<typename HEG>
class Get_beta<HEG, 0>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  { return CGAL::prev(dh, heg); }
};
template<typename HEG>
class Get_beta<HEG, 1>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  { return CGAL::next(dh, heg); }
};
template<typename HEG>
class Get_beta<HEG, 2>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  {
    if (Is_free<HEG, 2>::value(heg, dh)) return nullptr;
    return CGAL::opposite(dh, heg);
  }
};
// ////////////////////////////////////////////////////////////////////////////////
// /// Belong_to_same_cell
// template<typename HEG, unsigned int i>
// class Belong_to_same_cell
// {
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   static bool value(const HEG& /*heg*/,
//                     Dart_const_handle /*dh1*/, Dart_const_handle /*dh2*/) 
//   {  /* CGAL_static_assertion(false);*/
//     std::cout<<"ERROR Belong_to_same_cell<HEG, "<<i<<">"<<std::endl;
//     CGAL_assertion(false);
//     return false;
//   }
// };
// template<typename HEG>
// class Belong_to_same_cell<HEG, 0>
// {
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   static bool value(const HEG& heg,
//                           Dart_const_handle dh1, Dart_const_handle dh2) 
//   {
//     bool onborder=false;
//     Dart_const_handle curd=dh1;
//     do
//     {
//       if (curd==dh2) { return true; }
//       if (Is_free<HEG, 2>::value(heg, curd)) { onborder=true; }
//       else { curd=Get_beta<HEG, 1>::value(heg, Get_beta<HEG, 2>(curd)); }
//     }
//     while(curd!=dh1 && !onborder);

//     if (onborder)
//     { // We need to iterate in the other way
//       Dart_const_handle curd=dh1;
//       do
//       {
//         if (curd==dh2) { return true; }
//         if (Is_free<HEG, 2>::value(heg, curd)) { onborder=true; }
//         else { curd=Get_beta<HEG, 2>::value(heg, Get_beta<HEG, 0>(curd)); }
//         CGAL_assertion(curd!=dh1);
//       }
//       while(!onborder);
//     }
    
//     return false;
//   }
// };
// template<typename HEG>
// class Belong_to_same_cell<HEG, 1>
// {
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   static bool value(const HEG& heg,
//                     Dart_const_handle dh1, Dart_const_handle dh2) 
//   {
//     return dh1==dh2 ||
//       (!Is_free<HEG, 2>::value(dh1) && Get_beta<HEG, 2>::value(dh1)==dh2);
//   }
// };
// template<typename HEG>
// class Belong_to_same_cell<HEG, 2>
// {
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   static bool value(const HEG& heg,
//                     Dart_const_handle dh1, Dart_const_handle dh2) 
//   {
//     Dart_const_handle curd=dh1;
//     do
//     {
//       if (curd==dh2) { return true; }
//       curd=Get_beta<HEG, 1>::value(heg, curd);
//     }
//     while(curd!=dh1);
//     return false;
//   }
// };
////////////////////////////////////////////////////////////////////////////////

#endif // CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H //
// EOF //
