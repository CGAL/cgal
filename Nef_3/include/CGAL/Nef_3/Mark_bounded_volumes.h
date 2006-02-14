// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef MARK_BOUNDED_VOLUMES_H
#define MARK_BOUNDED_VOLUMES_H

#include <CGAL/Nef_polyhedron_3.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_3>
class Mark_bounded_volumes : public Modifier_base<typename Nef_3::SNC_structure>
{  typedef typename Nef_3::SNC_structure      SNC_structure;
   typedef typename SNC_structure::Infi_box   Infi_box;
   typedef typename Nef_3::Volume_iterator    Volume_iterator;

   bool flag;

public:
   Mark_bounded_volumes (bool b=true) : flag(b) {}

   void operator()(SNC_structure &snc)
   {  // mark bounded volumes
      Volume_iterator vol_it = snc.volumes_begin();
      CGAL_assertion ( vol_it != snc.volumes_end() );
      if ( Infi_box::extended_kernel() ) ++vol_it; // skip Infi_box
      CGAL_assertion ( vol_it != snc.volumes_end() );
      ++vol_it; // skip unbounded volume
      for (; vol_it != snc.volumes_end(); ++vol_it)
      {  vol_it->mark() = flag; // mark
      }
      // simplify(); // in Nef_3.delegate()
   }
};

CGAL_END_NAMESPACE
#endif // MARK_BOUNDED_VOLUMES_H
