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
// Author(s)     : Eti Ezra          <estere@post.tau.ac.il>
#ifndef CGAL_MAP_OVERLAY_BASE_H
#define CGAL_MAP_OVERLAY_BASE_H

CGAL_BEGIN_NAMESPACE

template <class Arrangement_, class Map_overlay_ChangeNotification_>
class Map_overlay_base
{
public:
  typedef Arrangement_                     Arrangement;
  typedef Map_overlay_ChangeNotification_  Map_overlay_ChangeNotification;

  Map_overlay_base() {}
  
  virtual void map_overlay(const Arrangement &a1, 
                           const Arrangement &a2, 
                           Map_overlay_ChangeNotification *pmwx_change_notf, 
                           Arrangement &result) = 0;

  virtual ~Map_overlay_base() {}
};

CGAL_END_NAMESPACE

#endif
