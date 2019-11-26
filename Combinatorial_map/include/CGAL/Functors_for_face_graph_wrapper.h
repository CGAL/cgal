// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H
#define CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H 1

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Iterators_for_face_graph_wrapper.h>

////////////////////////////////////////////////////////////////////////////////
/** This file contains the following functors for Face_graph_wrapper:
 * Is_free<typename HEG, unsigned int i>::
         operator() (const HEG& heg, Dart_const_handle dh)
 * Get_beta<typename HEG, unsigned int i>::
         operator() (const HEG& heg, Dart_const_handle dh)
*/
////////////////////////////////////////////////////////////////////////////////
namespace CGAL 
{
/// Is_free
//template<typename HEG, unsigned int i>
//struct Is_free
//{
//  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//  static bool value(const HEG& /*heg*/, Dart_const_handle /*dh*/)
//  { CGAL_static_assertion(i==0 || i==1); return false; }
//};
//template<typename HEG>
//struct Is_free<HEG, 2>
//{
//  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
//  static bool value(const HEG& heg, Dart_const_handle dh)
//  { return is_border(opposite(dh, heg), heg); }
//};
////////////////////////////////////////////////////////////////////////////////
/// Get_beta
template<typename HEG, unsigned int i>
struct Get_beta
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
struct Get_beta<HEG, 0>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  { return prev(dh, heg); }
};
template<typename HEG>
struct Get_beta<HEG, 1>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  { return next(dh, heg); }
};
template<typename HEG>
struct Get_beta<HEG, 2>
{
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  static Dart_const_handle value(const HEG& heg, Dart_const_handle dh) 
  {
    //if (Is_free<HEG, 2>::value(heg, dh)) return Dart_const_handle();
    return opposite(dh, heg);
  }
};
////////////////////////////////////////////////////////////////////////////////
} // namespace CGAL
  
#endif // CGAL_FUNCTORS_FOR_FACE_GRAPH_WRAPPER_H //
// EOF //
