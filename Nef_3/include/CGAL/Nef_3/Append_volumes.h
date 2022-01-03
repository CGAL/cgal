// Copyright (c) 2022  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//

#ifndef CGAL_NEF3_APPEND_VOLUMES_H
#define CGAL_NEF3_APPEND_VOLUMES_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/shell_to_nef_3.h>

namespace CGAL {

template<typename Nef_3>
class Append_volumes : public Modifier_base<typename Nef_3::SNC_structure>
{
  typedef typename Nef_3::SNC_structure SNC_structure;
  typedef typename Nef_3::Shell_entry_const_iterator Shell_entry_const_iterator;
public:
   Append_volumes (const Nef_3& n) : nef(n) {}

   void operator()(SNC_structure& snc)
   {
     Shell_entry_const_iterator si;
     CGAL_forall_shells_of(si,nef.volumes_begin())
       CGAL::shell_to_nef_3(nef,si,snc,true);
   }
private:
  const Nef_3& nef;
};

} //namespace CGAL
#endif // CGAL_NEF3_APPEND_VOLUMES_H
