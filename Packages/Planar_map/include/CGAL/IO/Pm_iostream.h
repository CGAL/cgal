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
// file          : include/CGAL/IO/Pm_iostream.h
// package       : pm (5.45)
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

#ifndef CGAL_PM_IOSTREAM_H
#define CGAL_PM_IOSTREAM_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_H
#include <CGAL/Planar_map_2.h>
#endif

//#ifndef CGAL_INVERSE_INDEX_H
//#include <CGAL/Inverse_index.h>
//#endif

#ifndef CGAL_IO_PM_FILE_WRITER_H
#include <CGAL/IO/Pm_file_writer.h>
#endif // CGAL_IO_PM_FILE_WRITER_H

#ifndef CGAL_IO_WRITE_PM_H
#include <CGAL/IO/write_pm.h>
#endif // CGAL_IO_WRITE_PM_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Dcel, class Traits> inline
::std::ostream& operator << (::std::ostream& o, const Planar_map_2<Dcel,Traits>& pm) 
{

  Pm_file_writer< Planar_map_2<Dcel,Traits> >  writer(o, pm);
  
  write_pm(pm, writer, o);
  
  return o;
}
 
template <class Dcel, class Traits> inline
::std::istream& operator >> ( std::istream& in, Planar_map_2<Dcel,Traits>& pm) {
  
  pm.read(std::cin);

  return in;
}

CGAL_END_NAMESPACE


#endif











