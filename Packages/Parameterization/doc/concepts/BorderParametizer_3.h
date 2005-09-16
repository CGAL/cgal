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


/// Concept BorderParametizer_3 of a class that parametrizes a given type of mesh,
/// 'Adaptor', which is a model of the MeshAdaptor_3 concept.
///
/// Implementation note:
/// To simplify the implementation, BorderParametizer_3 models know only the
/// MeshAdaptor_3 class. They don't know the parameterization algorithm
/// requirements nor the kind of sparse linear system used.
///
/// Design pattern:
/// BorderParametizer_3 models are Strategies (see [GOF95]): they implement
/// a strategy of boundary parameterization for models of MeshAdaptor_3.

class BorderParametizer_3
{
// Public types
public:
    typedef xxx Adaptor;
    typedef typename Parametizer_traits_3<Adaptor>::ErrorCode
                                            ErrorCode;

// Public operations
public:
    /// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
    /// Return false on error.
    ErrorCode parameterize_border (Adaptor* mesh);

    /// Indicate if border's shape is convex
    bool  is_border_convex ();
};

