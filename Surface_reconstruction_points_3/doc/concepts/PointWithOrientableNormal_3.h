// Copyright (c) 2007-2008  INRIA (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez


/// A PointWithOrientableNormal_3 concept represents a 3D point with:
/// - a position
/// - a normal (oriented or not).
///
/// It is strictly identical to PointWithNormal_3 concept, except that its
/// normal type is a model of OrientableNormal_3.
///
/// @heading Has Models:
/// - Point_with_normal_3<GeomTraits, Normal_3> if Normal_3 is a model of OrientableNormal_3.

class PointWithOrientableNormal_3 
  : public PointWithNormal_3,
    public DefaultConstructible, public CopyConstructible, public Assignable, public EqualityComparable
{
public:
    typedef xxx Normal; ///< Model of OrientableNormal_3 concept.

    // Everything else is identical to PointWithNormal_3 concept.
};

