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
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Map_overlay_base.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
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
