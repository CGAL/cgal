// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj

#ifndef CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
#define CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/tags.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_triangulation_3`
is a class template which provides a valid triangulation type
that can be used as the 3D triangulation input for
the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Concurrency_tag enables sequential versus parallel implementation of the
triangulation data structure.
Possible values are `Sequential_tag` (the default), `Parallel_tag`,
and `Parallel_if_available_tag`.

\tparam Vb is a vertex base class from which `Remeshing_vertex_base_3` derives.
It must be a model of the `TriangulationVertexBase_3` concept.
It has the default value `Triangulation_vertex_base_3<Gt>`.

\tparam Cb is a cell base class from which `Remeshing_cell_base_3` derives.
It must be a model of the `TriangulationCellBase_3` concept.
It has the default value `Triangulation_cell_base_3<Gt>`.

*/
template<typename Gt,
         typename Concurrency_tag = CGAL::Sequential_tag,
         typename Vb = CGAL::Triangulation_vertex_base_3<Gt>,
         typename Cb = CGAL::Triangulation_cell_base_3<Gt>
>
class Remeshing_triangulation_3
  : public CGAL::Triangulation_3<Gt,
      CGAL::Triangulation_data_structure_3<
        Remeshing_vertex_base_3<Gt, Vb>,
        Remeshing_cell_base_3<Gt, Cb>
      >
    >
{
public:
  typedef Remeshing_vertex_base_3<Gt, Vb> Remeshing_Vb;
  typedef Remeshing_cell_base_3<Gt, Cb>   Remeshing_Cb;

  typedef CGAL::Triangulation_data_structure_3<
            Remeshing_Vb, Remeshing_Cb, Concurrency_tag>  Tds;
  typedef CGAL::Triangulation_3<Gt, Tds>                  Base;

  using Base::Base;
};

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
