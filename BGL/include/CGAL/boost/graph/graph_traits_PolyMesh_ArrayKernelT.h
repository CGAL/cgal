// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri, Philipp Moeller

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H

// http://openmesh.org/Documentation/OpenMesh-Doc-Latest/classOpenMesh_1_1Concepts_1_1KernelT.html
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#define OPEN_MESH_CLASS OpenMesh::PolyMesh_ArrayKernelT<K>
#include <CGAL/boost/graph/graph_traits_OpenMesh.h>

#endif // CGAL_BOOST_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H
