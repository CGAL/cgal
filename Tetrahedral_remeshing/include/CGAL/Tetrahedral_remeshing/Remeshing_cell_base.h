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

#ifndef CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_H
#define CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_H

#include <CGAL/Mesh_cell_base_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
struct Fake_MD_C
{
  typedef int Subdomain_index;
  typedef int Surface_patch_index;
  typedef int Index;
};
}

/*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_cell_base` is a model of the concept `MeshCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Cb is a cell base class from which `Remeshing_cell_base` derives.
It must be a model of the `TriangulationCellBase_3` concept.
It has the default value `Triangulation_cell_base_3<Gt>`.

\cgalModels `MeshCellBase_3`

*/
template<typename Gt,
         typename Cb = CGAL::Triangulation_cell_base_3<Gt> >
class Remeshing_cell_base
#ifndef DOXYGEN_RUNNING
  : public CGAL::Mesh_cell_base_3<Gt, internal::Fake_MD_C, Cb>
#endif
{
  typedef CGAL::Mesh_cell_base_3<Gt, internal::Fake_MD_C, Cb> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle   Cell_handle;

public:
  // To get correct cell type in TDS
  template < class TDS2 >
  struct Rebind_TDS
  {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
    typedef Remeshing_cell_base<Gt, Cb2> Other;
  };

  using Base::Base;

#ifndef DOXYGEN_RUNNING
  /// TODO : remove this function from here
  /// Returns `true` if facet lies on a surface patch
  bool is_facet_on_surface(const int facet) const
  {
    CGAL_precondition(facet >= 0 && facet<4);
    return this->subdomain_index() != this->neighbor(facet)->subdomain_index();
  }
#endif
};



template < class Gt, class Cb >
std::istream&
operator>>(std::istream &is, Remeshing_cell_base<Gt, Cb> &c)
{
  typename Remeshing_cell_base<Gt, Cb>::Subdomain_index index;
  if (is_ascii(is))
    is >> index;
  else
    read(is, index);
  if (is) {
    c.set_subdomain_index(index);
//      for (int i = 0; i < 4; ++i)
//      {
//        typename Compact_mesh_cell_base_3<GT, MT, Cb>::Surface_patch_index i2;
//        if (is_ascii(is))
//          is >> iformat(i2);
//        else
//        {
//          read(is, i2);
//        }
//        c.set_surface_patch_index(i, i2);
//      }
  }
  return is;
}

template < class Gt, class Cb >
std::ostream&
operator<<(std::ostream &os, const Remeshing_cell_base<Gt, Cb> &c)
{
  if (is_ascii(os))
    os << c.subdomain_index();
  else
    write(os, c.subdomain_index());
  //for (int i = 0; i < 4; ++i)
  //{
  //  if (is_ascii(os))
  //    os << ' ' << oformat(c.surface_patch_index(i));
  //  else
  //    write(os, c.surface_patch_index(i));
  //}
  return os;
}


}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_H
