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
::std::ostream & operator << (::std::ostream & o,
                              const Planar_map_2<Dcel,Traits> & pm) 
{

  Pm_file_writer< Planar_map_2<Dcel,Traits> >  writer(o, pm);
  
  write_pm(pm, writer, o);
  
  return o;
}
 
template <class Dcel, class Traits> inline
::std::istream & operator >> (std::istream& in, Planar_map_2<Dcel,Traits> & pm)
{
  
  pm.read(std::cin);

  return in;
}

CGAL_END_NAMESPACE


#endif











