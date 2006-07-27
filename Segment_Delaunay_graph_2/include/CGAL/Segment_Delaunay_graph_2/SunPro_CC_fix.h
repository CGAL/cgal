// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_SUNPRO_CC_FIX_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_SUNPRO_CC_FIX_H 1

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian_converter.h>
//#include <CGAL/number_utils_classes.h>


CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

#if defined(__sun) && defined(__SUNPRO_CC)
// workaround for the Sun CC-5.30 compiler; it does not like default
// template parameters that are themselves templates and have
// templated classes as parameters, which have then nested types as
// arguments... oooof!!!
//
// In case you did understand what I just described you are most
// probably crazy... If you did not, look below to see what kind of
// code CC-5.30 did not like.
namespace Internal {

  template<class CK, class FK>
  struct SUNPRO_CC_Interval_converter
    : public Cartesian_converter<CK, FK,
                                 To_interval< typename CK::RT > >
  {
  };

}
#endif

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_SUNPRO_CC_FIX_H
