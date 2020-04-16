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

#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base.h>

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
Possible values are `Sequential_tag` (the default) and `Parallel_tag`.

\tparam Cb is a cell base class from which `Remeshing_cell_base` derives.
It must be a model of the `TriangulationCellBase_3` concept.
It has the default value `Triangulation_cell_base_3<Gt>`.

\tparam Vb is a vertex base class from which `Remeshing_vertex_base` derives.
It must be a model of the `TriangulationVertexBase_3` concept.
It has the default value `Triangulation_vertex_base_3<Gt>`.

\cgalRefines `Triangulation_3`

*/
template<typename Gt,
         typename Concurrency_tag = CGAL::Sequential_tag,
         typename Cb = CGAL::Triangulation_cell_base_3<Gt>,
         typename Vb = CGAL::Triangulation_vertex_base_3<Gt>
>
class Remeshing_triangulation_3
  : public CGAL::Triangulation_3<Gt,
      CGAL::Triangulation_data_structure_3<
        Remeshing_vertex_base<Gt, Vb>,
        Remeshing_cell_base<Gt, Cb>
      >
    >
{
public:
  typedef Remeshing_vertex_base<Gt, Vb> Remeshing_Vb;
  typedef Remeshing_cell_base<Gt, Cb>   Remeshing_Cb;

  typedef CGAL::Triangulation_data_structure_3<
            Remeshing_Vb, Remeshing_Cb, Concurrency_tag>  Tds;
  typedef CGAL::Triangulation_3<Gt, Tds>                  Base;

  using Base::Base;
};

namespace internal
{
template<typename TDS_src, typename TDS_tgt>
struct Vertex_converter
{
  //This operator is used to create the vertex from v_src.
  typename TDS_tgt::Vertex operator()(const typename TDS_src::Vertex& v_src) const
  {
    typedef typename CGAL::Kernel_traits<
    typename TDS_src::Vertex::Point>::Kernel GT_src;
    typedef typename CGAL::Kernel_traits<
    typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
    CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

    typedef typename TDS_tgt::Vertex::Point Tgt_point;

    typename TDS_tgt::Vertex v_tgt;
    v_tgt.set_point(Tgt_point(conv(point(v_src.point()))));
    v_tgt.set_time_stamp(-1);
    v_tgt.set_dimension(3);//-1 if unset, 0,1,2, or 3 if set
    return v_tgt;
  }
  //This operator is meant to be used in case heavy data should transferred to v_tgt.
  void operator()(const typename TDS_src::Vertex& v_src,
                  typename TDS_tgt::Vertex& v_tgt) const
  {
    typedef typename CGAL::Kernel_traits<
    typename TDS_src::Vertex::Point>::Kernel GT_src;
    typedef typename CGAL::Kernel_traits<
    typename TDS_tgt::Vertex::Point>::Kernel GT_tgt;
    CGAL::Cartesian_converter<GT_src, GT_tgt> conv;

    typedef typename TDS_tgt::Vertex::Point Tgt_point;

    v_tgt.set_point(Tgt_point(conv(point(v_src.point()))));
    v_tgt.set_dimension(3);//v_src.info());
  }
};

template<typename TDS_src, typename TDS_tgt>
struct Cell_converter
{
  //This operator is used to create the cell from c_src.
  typename TDS_tgt::Cell operator()(const typename TDS_src::Cell& c_src) const
  {
    typename TDS_tgt::Cell c_tgt;
    c_tgt.set_subdomain_index(1);//c_src.subdomain_index());
//        c_tgt.info() = c_src.info();
    c_tgt.set_time_stamp(-1);
    return c_tgt;
  }
  //This operator is meant to be used in case heavy data should transferred to c_tgt.
  void operator()(const typename TDS_src::Cell& c_src,
                  typename TDS_tgt::Cell& c_tgt) const
  {
//        c_tgt.set_subdomain_index(c_src.subdomain_index());
    //        c_tgt.info() = c_src.info();
  }
};

}


template<typename T3, typename Gt>
void build_remeshing_triangulation(const T3& tr,
                                   Remeshing_triangulation_3<Gt>& remeshing_tr)
{
  typedef typename T3::Triangulation_data_structure Tds;
  typedef typename Remeshing_triangulation_3<Gt>::Tds   RTds;

  remeshing_tr.clear();

  remeshing_tr.set_infinite_vertex(
    remeshing_tr.tds().copy_tds(
      tr.tds(),
      tr.infinite_vertex(),
      internal::Vertex_converter<Tds, RTds>(),
      internal::Cell_converter<Tds, RTds>()));
}

template<typename T3, typename Gt>
void build_from_remeshing_triangulation(
  const Remeshing_triangulation_3<Gt>& remeshing_tr,
  T3& tr)
{
  typedef typename T3::Triangulation_data_structure Tds;
  typedef typename Remeshing_triangulation_3<Gt>::Tds   RTds;

  tr.clear();

  tr.set_infinite_vertex(
    tr.tds().copy_tds(
      remeshing_tr.tds(),
      remeshing_tr.infinite_vertex(),
      internal::Vertex_converter<RTds, Tds>(),
      internal::Cell_converter<RTds, Tds>()));
}

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
