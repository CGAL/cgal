// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel        <seel@mpi-sb.mpg.de>
//                 Miguel Granados     <granados@mpi-sb.mpg.de>
//                 Susan Hert          <hert@mpi-sb.mpg.de>
//                 Lutz Kettner        <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NEF_VOLUME_H
#define CGAL_NEF_VOLUME_H

#include <CGAL/license/Nef_3.h>


#include <string>
#include <sstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename Refs>
class Volume_base  {

  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Volume_handle  Volume_handle;
  typedef typename Refs::Volume_const_handle  Volume_const_handle;
  typedef typename Refs::Object_list   Object_list;
  typedef typename Refs::Shell_entry_iterator
    Shell_entry_iterator;
  typedef typename Refs::Shell_entry_const_iterator
    Shell_entry_const_iterator;

  Mark         mark_;
  Object_list shell_entry_objects_; // SFaces

 public:

  Volume_base() {}

  Volume_base(Mark m) : mark_(m) {}

    ~Volume_base() {
      CGAL_NEF_TRACEN("  destroying Volume_base item "<<&*this);
    }

    Volume_base(const Volume_base<Refs>& v)
      { mark_ = v.mark_;
	shell_entry_objects_ = v.shell_entry_objects_;
      }

    Volume_base<Refs>& operator=(const Volume_base<Refs>& v)
      { if (this == &v) return *this;
	mark_ = v.mark_;
	shell_entry_objects_ = v.shell_entry_objects_;
	return *this;
      }

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Object_list& shell_entry_objects() { return shell_entry_objects_; }
    const Object_list& shell_entry_objects() const { 
      return shell_entry_objects_; 
    }

    Shell_entry_iterator shells_begin()
    { return shell_entry_objects_.begin(); }
    Shell_entry_iterator shells_end()
    { return shell_entry_objects_.end(); }
    Shell_entry_const_iterator shells_begin() const
    { return shell_entry_objects_.begin(); }
    Shell_entry_const_iterator shells_end() const
    { return shell_entry_objects_.end(); }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::Volume_base::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid = (!shell_entry_objects_.empty());

      verr << "end of CGAL::SNC_items<...>::Volume_base::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

}; // Volume_base

} //namespace CGAL
#endif //CGAL_NEF_VOLUME_H
