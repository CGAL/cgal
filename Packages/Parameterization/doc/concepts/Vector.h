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


/// Concept SparseLinearAlgebraTraits_d::Vector
/// is a concept of a vector that can be multiplied by a sparse matrix. 
///
/// Sub-concept: This is a sub-concept of LinearAlgebraTraits_d::Vector.
///
/// Models:
/// - Taucs_vector
/// - OpenNL::FullVector

class Vector
{
// Public types
public:

    typedef xxx NT;

// Public operations
public:

    /// Create a vector initialized with zeros
    Vector (int rows);

    /// Copy constructor
    Vector(const Vector& toCopy);

    /// operator =()
    Vector& operator=(const Vector& toCopy);

    /// Return the vector's number of coefficients
    int  dimension () const;

    /// Read/write access to 1 vector coefficient.
    ///
    /// Precondition: 0 <= row < dimension().
    NT  operator[] (int row) const;
    NT&  operator[] (int row);
};

