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
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_euclidean_traits_wrapper_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H
#define CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H

CGAL_BEGIN_NAMESPACE


template<class Gt_base>
class Apollonius_graph_gt_wrapper : public Gt_base
{
public:
  struct Triangle_2 {};

  Apollonius_graph_gt_wrapper() {}
  Apollonius_graph_gt_wrapper(const Gt_base& gtb) : Gt_base(gtb) {}

};




CGAL_END_NAMESPACE


#endif // CGAL_APOLLONIUS_GRAPH_EUCLIDEAN_TRAITS_WRAPPER_2_H
