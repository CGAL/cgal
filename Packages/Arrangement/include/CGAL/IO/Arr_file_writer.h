// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

#ifndef CGAL_IO_ARR_FILE_WRITER_H
#define CGAL_IO_ARR_FILE_WRITER_H 1

#include <CGAL/IO/binary_file_io.h>
#include <CGAL/IO/Pm_file_writer.h>
#include <iostream>
#include <cstddef>

CGAL_BEGIN_NAMESPACE

template <class Arrangement>
class Arr_file_writer : public  Pm_file_writer<Arrangement> {

public:
  typedef Pm_file_writer<Arrangement>                   Base;
  typedef typename Arrangement::Curve_iterator          Curve_iterator;
  typedef typename Arrangement::Subcurve_iterator       Subcurve_iterator;
  typedef typename Arrangement::Edge_iterator           Edge_iterator;
  typedef typename Arrangement::Curve_const_iterator    Curve_const_iterator;
  typedef typename Arrangement::Subcurve_const_iterator   
                                                        Subcurve_const_iterator;
  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Base::out;
#endif
  Arr_file_writer(std::ostream & o,
                  const Arrangement & arr,
                  bool verbose = false) : 
    Pm_file_writer<Arrangement>(o, arr, verbose) {}

  Arr_file_writer(std::ostream & o, const File_header& h) : 
    Pm_file_writer<Arrangement>(o, h) {}
  
  void write_curve (Curve_iterator cv){
    out () << cv->curve() << std::endl;
  }

  void write_curve (Curve_const_iterator cv){
    out () << cv->curve() << std::endl;
  }

  void write_subcurve (Subcurve_iterator scv){
    out () << scv->x_curve() << std::endl;
  }

  void write_subcurve (Subcurve_const_iterator scv){
    out () << scv->x_curve() << std::endl;
  }
  
  void write_edge(Edge_iterator edge){
    out () << edge->x_curve() << std::endl;
  }
  
  void write_edge(Edge_const_iterator edge){
    out () << edge->x_curve() << std::endl;
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
