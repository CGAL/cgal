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
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/IO/Arr_file_writer.h
// package       : Arrangement (1.81)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_IO_ARR_FILE_WRITER_H
#define CGAL_IO_ARR_FILE_WRITER_H 1

#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H

#ifndef CGAL_IO_PM_FILE_WRITER_H
#include <CGAL/IO/Pm_file_writer.h>
#endif // CGAL_IO_PM_FILE_HEADER_H

#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

CGAL_BEGIN_NAMESPACE


template <class Arrangement>
class Arr_file_writer : public  Pm_file_writer<Arrangement> {

public:
  //typedef Arrangement_                                           Arrangement;
  typedef typename Arrangement::Curve_iterator            Curve_iterator;
  typedef typename Arrangement::Subcurve_iterator         Subcurve_iterator;
  typedef typename Arrangement::Edge_iterator             Edge_iterator;
  typedef typename Arrangement::Curve_const_iterator      Curve_const_iterator;
  typedef typename Arrangement::Subcurve_const_iterator   
                                                    Subcurve_const_iterator;
  typedef typename Arrangement::Edge_const_iterator      Edge_const_iterator;

  Arr_file_writer(std::ostream& o, 
                  const Arrangement& arr, 
                  bool verbose = false) : 
    Pm_file_writer<Arrangement>(o, arr, verbose) {}

  Arr_file_writer(std::ostream& o, 
                  const File_header& h) : 
    Pm_file_writer<Arrangement>(o, h) {}
  

  void write_curve (Curve_iterator cv){
    out () << cv->curve() << std::endl;
  }

  void write_curve (Curve_const_iterator cv){
    out () << cv->curve() << std::endl;
  }

  void write_subcurve (Subcurve_iterator scv){
    out () << scv->curve() << std::endl;
  }

  void write_subcurve (Subcurve_const_iterator scv){
    out () << scv->curve() << std::endl;
  }
  
  void write_edge(Edge_iterator edge){
    out () << edge->curve() << std::endl;
  }
  
  void write_edge(Edge_const_iterator edge){
    out () << edge->curve() << std::endl;
  }

  //void write_edge_nodes_end() {
  //  out() << std::endl;
  // }

  /*void write_footer() {
    if (m_header.comments())
    out() << "#------------------- End of Arrangement #";
    out() << std::endl;
    } */
};
CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_WRITER_ARR_H //
// EOF //









