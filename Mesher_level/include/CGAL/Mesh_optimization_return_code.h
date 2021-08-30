// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_OPTIMIZATION_RETURN_CODE_H
#define CGAL_MESH_OPTIMIZATION_RETURN_CODE_H

namespace CGAL {

enum Mesh_optimization_return_code
{
  MESH_OPTIMIZATION_UNKNOWN_ERROR=-1,
  BOUND_REACHED=0,
  TIME_LIMIT_REACHED,
  CANT_IMPROVE_ANYMORE,
  CONVERGENCE_REACHED,
  MAX_ITERATION_NUMBER_REACHED,
  ALL_VERTICES_FROZEN
};


} //namespace CGAL

#endif // CGAL_MESH_OPTIMIZATION_RETURN_CODE_H
