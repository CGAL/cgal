// Copyright (c) 2001  Tel-Aviv University (Israel).
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

#ifndef CGAL_IO_ARR_FILE_SCANNER_H
#define CGAL_IO_ARR_FILE_SCANNER_H    1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#include <CGAL/known_bit_size_integers.h>
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H

//#ifndef CGAL_IO_FILE_HEADER_PM_H
//#include <CGAL/IO/File_header_pm.h>
//#endif // CGAL_IO_FILE_HEADER_PM_H

#ifndef CGAL_IO_PM_FILE_SCANNER_H
#include <CGAL/IO/Pm_file_scanner.h>
#endif // CGAL_IO_PM_FILE_SCANNER_H

#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

CGAL_BEGIN_NAMESPACE

template <class Arrangement>
class Arr_file_scanner : public  Pm_file_scanner<Arrangement> {
public:
  typedef Pm_file_scanner<Arrangement> Base;
  typedef typename Arrangement::Curve_node                Curve_node;
  typedef typename Arrangement::Subcurve_node             Subcurve_node;
  typedef typename Arrangement::Edge_node                 Edge_node;
  
  typedef typename Arrangement::Traits                   Traits;
  typedef typename Traits::Point                         Point;
  typedef typename Traits::X_curve                       X_curve;
  typedef typename Traits::Curve                         Curve;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Base::skip_comment;
  using Base::in;
#endif


  Arr_file_scanner(std::istream& in) : Pm_file_scanner<Arrangement>(in) {}

  Arr_file_scanner(std::istream& in, const File_header& h) : 
    Pm_file_scanner<Arrangement>(in, h) {}

  void scan_Curve_node(Curve_node* cn){

    skip_comment();

    // providing default reading function.
    Curve curve;
    in() >> curve;

    cn->set_curve(curve);
  }
  
  void scan_Subcurve_node(Subcurve_node* scn){
    
    skip_comment();

    // providing default reading function.
    X_curve curve;
    in() >> curve;
    
    scn->set_x_monotone_curve(curve);
  }
    
  void scan_Edge_node(Edge_node* en){
    
    skip_comment();

    // providing default reading function.
    X_curve curve;
    in() >> curve;
    
    en->set_x_monotone_curve(curve);
  }

};


CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_SCANNER_PM_H //
// EOF //









