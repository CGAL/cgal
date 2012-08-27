// Copyright (c) 2005  Tel-Aviv University (Israel).
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
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>
//                 Efi Fogel              <efif@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_BASE_H
#define CGAL_ENVELOPE_BASE_H

namespace CGAL {

// Envelope types:
enum Envelope_type {
  ENVELOPE_LOWER = 1,
  ENVELOPE_UPPER = 2
};

// Types of intersection curve between 2 xy-monotone surfaces:
enum Intersection_type {
  INTERSECTION_UNKNOWN = 0,
  INTERSECTION_TRANSVERSAL = 1,
  INTERSECTION_TANGENT = 2
};

// Decision mark for DCEL features:
enum Dac_decision
{
  DAC_DECISION_FIRST = -1,
  DAC_DECISION_BOTH,
  DAC_DECISION_SECOND,
  DAC_DECISION_NOT_SET
};

} //namespace CGAL

#endif // CGAL_ENVELOPE_BASE_3_H
