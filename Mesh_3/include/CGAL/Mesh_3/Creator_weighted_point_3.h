// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// This file contains the definition and implementation of class
// Creator_weighted_point_3. See class description for further information.
//******************************************************************************

#ifndef CGAL_MESH_3_CREATOR_WEIGHTED_POINT_3_H
#define CGAL_MESH_3_CREATOR_WEIGHTED_POINT_3_H

namespace CGAL {


namespace Mesh_3 {


/**
 * \class Creator_weighted_point_3
 *
 * This class is designed to handle the creation of Weighted_point from its 3
 * coordinates only.
 *
 * It is a model of concept Creator, used by many Geometric Object Generators
 */
template<class Arg_, class Weighted_point>
class Creator_weighted_point_3
{
public:
  typedef Arg_ argument_type;

  /// Constructor
  Creator_weighted_point_3() {}
  /// Destructor
  ~Creator_weighted_point_3() {}

  /// Constructs a weighted point with default weight from 3 \c Arg_
  /// We do not pass \c const \c Arg_& because \c Arg_ is a ring number type
  Weighted_point operator()(Arg_ a1, Arg_ a2, Arg_ a3) const
  {
    typedef typename Weighted_point::Point Bare_point;
    return Weighted_point(Bare_point(a1,a2,a3));
  }

private:
  // Disabled copy constructor & assignment operator
  typedef Creator_weighted_point_3<Arg_, Weighted_point> Self;
  Creator_weighted_point_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Creator_weighted_point_3


}  // end namespace Mesh_3


}  // end namespace CGAL

#endif // CGAL_MESH_3_CREATOR_WEIGHTED_POINT_3_H
