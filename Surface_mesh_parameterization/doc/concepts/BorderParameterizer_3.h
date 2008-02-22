// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


/// BorderParameterizer_3 is a concept of class that parameterizes a given type of mesh,
/// 'Adaptor', which is a model of the ParameterizationMesh_3 concept.
///
/// Implementation note:
/// To simplify the implementation, BorderParameterizer_3 models know only the
/// ParameterizationMesh_3 class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// @heading Design Pattern:
/// BorderParameterizer_3 models are Strategies [GHJV95]: they implement
/// a strategy of border parameterization for models of ParameterizationMesh_3.

class BorderParameterizer_3
{
// Public types
public:
    /// Export ParameterizationMesh_3 template parameter
    typedef xxx Adaptor;

    /// The various errors detected by this package
    typedef xxx Error_code;

// Public operations
public:
    // Construction and destruction are undefined

    /// Assign to mesh's border vertices a 2D position (i.e. a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
    /// Return false on error.
    Error_code parameterize_border (Adaptor& mesh);

    /// Indicate if border's shape is convex.
    bool  is_border_convex ();
};

