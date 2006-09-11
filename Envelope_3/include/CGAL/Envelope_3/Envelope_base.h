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
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_BASE_H
#define CGAL_ENVELOPE_BASE_H

CGAL_BEGIN_NAMESPACE

// Envelope types:
typedef unsigned int   Envelope_type;

const Envelope_type      LOWER = 1;
const Envelope_type      UPPER = 2;

// Types of intersection curve between 2 xy-monotone surfaces:
typedef unsigned int   Intersection_type;

const Intersection_type  UNKNOWN = 0;
const Intersection_type  TRANSVERSAL = 1;
const Intersection_type  TANGENT = 2;

// Decision mark for DCEL features:
enum Dac_decision
{
  FIRST = -1,
  BOTH,
  SECOND,
  NOT_SET
};

CGAL_END_NAMESPACE

#endif // CGAL_ENVELOPE_BASE_3_H
