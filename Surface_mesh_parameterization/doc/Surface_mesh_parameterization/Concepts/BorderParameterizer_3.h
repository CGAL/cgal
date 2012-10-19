// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


/// \ingroup PkgSurfaceParameterizationConcepts
/// \cgalconcept
/// BorderParameterizer_3 is a concept of class that parameterizes a given type of mesh,
/// <i>Adaptor</i>, which is a model of the ParameterizationMesh_3 concept.
///
/// Models of BorderParameterizer_3 know only the
/// ParameterizationMesh_3 class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///


class BorderParameterizer_3
{
// Public types
public:
    /// Export ParameterizationMesh_3 template parameter
    typedef Hidden_type Adaptor;

    /// The various errors detected by this package
  typedef Hidden_type Error_code;

// Public operations
public:
    // Construction and destruction are undefined

    /// Assign to mesh's border vertices a 2D position (i.e.\ a `(u,v)` pair)
    /// on border's shape. Mark them as <i>parameterized</i>.
    Error_code parameterize_border (Adaptor& mesh);

    /// Indicate if border's shape is convex.
    bool  is_border_convex ();
};

