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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


/// Concept of parameterization objects for a given type of mesh,
/// 'Adaptor', which is a model of the MeshAdaptor_3 concept.
class ParametizerTraits_3
{
// Public types
public:
    /// The various errors detected by this package
    enum ErrorCode {
    OK,
    ERROR_EMPTY_MESH,               ///< input mesh is empty
    ERROR_NON_TRIANGULAR_MESH,      ///< input mesh is not triangular
    ERROR_NO_SURFACE_MESH,          ///< input mesh is not a surface
    ERROR_INVALID_BOUNDARY,         ///< parameterization requires a convex border
    ERROR_BAD_MATRIX_CONDITIONING,  ///< result is mathematically unstable
    ERROR_CANNOT_SOLVE_LINEAR_SYSTEM,///< cannot solve linear system
    ERROR_NO_1_TO_1_MAPPING,        ///< parameterization does not ensure 1 to 1 mapping
    ERROR_NOT_ENOUGH_MEMORY,        ///< it's time to buy some RAM :-)
    ERROR_WRONG_PARAMETER           ///< a method received an unexpected parameter
    };

    typedef xxx Adaptor;
    typedef xxx NT;

// Public operations
public:
    /// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh'
    /// to a piece of the 2D space.
    /// The mapping is linear by pieces (linear in each triangle).
    /// The result is the (u,v) pair image of each vertex of the 3D surface.
    ///
    /// Preconditions:
    /// - 'mesh' must be a surface with 1 connected component and no hole
    /// - 'mesh' must be a triangular mesh
    ErrorCode  parameterize (Adaptor* mesh);
};

