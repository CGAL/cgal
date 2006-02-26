// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source: /CVSROOT/CGAL/Packages/Envelope_3/include/CGAL/Envelope_base.h,v $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_BASE_3_H
#define CGAL_ENVELOPE_BASE_3_H

CGAL_BEGIN_NAMESPACE

enum Envelope_type
{
  LOWER   = 1,
  UPPER
};

// the type of intersection curve between 2 xy-monotone surfaces
enum Intersection_type
{
  TRANSVERSAL = 1, // the 2 surfaces vertically change place at the intersection curve
  TANGENT,         // the 2 surfaces keep their relative vertical position
  UNKNOWN          // the type is not known
};

enum Dac_decision
{
  FIRST = -1,
  BOTH,
  SECOND,
  NOT_SET
};

//template <class _Traits>
//class Envelope_base_3
//{
//public:
//  typedef _Traits     Traits;
//
//  // virtual destructor.
//  virtual ~Envelope_base_3(){}
//
//};

CGAL_END_NAMESPACE

#endif // CGAL_ENVELOPE_BASE_3_H
