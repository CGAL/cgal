// Copyright (c) 2019 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_H
#define CGAL_TET_ADAPTIVE_REMESHING_CELL_BASE_H


#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include <utility>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  /*!
\ingroup PkgTetrahedralRemeshingClasses

The class `Remeshing_cell_base` is a model of the concept `RemeshingCellBase_3`.
It is designed to serve as cell base class for the 3D triangulation
used in the tetrahedral remeshing process.

\tparam Gt is the geometric traits class.
It has to be a model of the concept `RemeshingTriangulationTraits_3`.

\tparam Info is the information the user would like to add to a cell.
It has to be `DefaultConstructible` and `Assignable`.

\tparam Cb is a cell base class from which `Remeshing_cell_base` derives.
It must be a model of the `TriangulationCellBase_3` concept.
It has the default value `Triangulation_cell_base_3<Gt>`.

\cgalModels `RemeshingCellBase_3`
\cgalRefines `Triangulation_cell_base_with_info_3`

*/
  template<typename Gt,
           typename Info,
           typename Cb = CGAL::Triangulation_cell_base_3<Gt> >
  class Remeshing_cell_base
    : public CGAL::Triangulation_cell_base_with_info_3<Info, Gt, Cb>

  {
    typedef CGAL::Triangulation_cell_base_with_info_3<Info, Gt, Cb> Base;
    typedef typename Base::Vertex_handle Vertex_handle;
    typedef typename Base::Cell_handle   Cell_handle;

  public:
    typedef int Subdomain_index;
    typedef int Surface_patch_index;

  private:
    Subdomain_index subdomain_index_;
        //  0 is undefined
        // -1 for infinite cells
        //  1 to n for subdomains
        // n + 1 for imaginary cells
    std::size_t time_stamp_;

  public:
    // To get correct cell type in TDS
    template < class TDS2 >
    struct Rebind_TDS
    {
      typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
      typedef Remeshing_cell_base<Gt, Info, Cb2> Other;
    };

    Remeshing_cell_base()
      : subdomain_index_(0)
      , time_stamp_(-1)
    {}

    Remeshing_cell_base(Vertex_handle v0,
                        Vertex_handle v1,
                        Vertex_handle v2,
                        Vertex_handle v3)
      : Base(v0, v1, v2, v3)
      , subdomain_index_(0)
      , time_stamp_(-1)
    {}

    Remeshing_cell_base(Vertex_handle v0,
                        Vertex_handle v1,
                        Vertex_handle v2,
                        Vertex_handle v3,
                        Cell_handle n0,
                        Cell_handle n1,
                        Cell_handle n2,
                        Cell_handle n3)
      : Base(v0, v1, v2, v3, n0, n1, n2, n3)
      , subdomain_index_(0)
      , time_stamp_(-1)
    {}

    const Subdomain_index& subdomain_index() const
    {
      return subdomain_index_;
    }
    void set_subdomain_index(const Subdomain_index& si)
    {
      subdomain_index_ = si;
    }

    void set_surface_patch_index(const int, const Surface_patch_index&)
    {/*nothing to do because we use incident subdomain indices*/ }

    const Surface_patch_index surface_patch_index(const int& i) const
    {
      CGAL_precondition(i >= 0 && i < 4);
      if(is_facet_on_surface(i))
        return 1;
      else
        return 0;
    }

    /// Returns `true` if facet lies on a surface patch
    bool is_facet_on_surface(const int& facet) const
    {
      CGAL_precondition(facet >= 0 && facet<4);
      return this->subdomain_index() != this->neighbor(facet)->subdomain_index();
    }

    typedef Tag_true Has_timestamp;
    std::size_t time_stamp() const {
      return time_stamp_;
    }
    void set_time_stamp(const std::size_t& ts) {
      time_stamp_ = ts;
    }

 };



  template < class Gt, class Info, class Cb >
  std::istream&
    operator>>(std::istream &is, Remeshing_cell_base<Gt, Info, Cb> &c)
  {
    typename Remeshing_cell_base<Gt, Info, Cb>::Subdomain_index index;
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

  template < class Gt, class Info, class Cb >
  std::ostream&
    operator<<(std::ostream &os, const Remeshing_cell_base<Gt, Info, Cb> &c)
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
