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


#ifndef CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H
#define CGAL_TETRAHEDRAL_REMESHING_TRIANGULATION_H

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  template<typename K,
           typename Info = void,
           typename Cb = CGAL::Triangulation_cell_base_3<K> >
  class Remeshing_triangulation_3
    : public CGAL::Triangulation_3<K,
        CGAL::Triangulation_data_structure_3<
          Remeshing_vertex_base<K>,
          Remeshing_cell_base<K, Info, Cb>
        >
      >
  {
    typedef Remeshing_vertex_base<K>                       RVb;
    typedef Remeshing_cell_base<K, Info, Cb>               RCb;

  public:
    typedef CGAL::Triangulation_data_structure_3<RVb, RCb> Tds;
    typedef CGAL::Triangulation_3<K, Tds>                  Self;
    typedef Self                                           type;
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
        c_tgt.set_subdomain_index(c_src.subdomain_index());
//        c_tgt.info() = c_src.info();
        c_tgt.set_time_stamp(-1);
        return c_tgt;
      }
      //This operator is meant to be used in case heavy data should transferred to c_tgt.
      void operator()(const typename TDS_src::Cell& c_src,
                      typename TDS_tgt::Cell& c_tgt) const
      {
        c_tgt.set_subdomain_index(c_src.subdomain_index());
        //        c_tgt.info() = c_src.info();
      }
    };

  }

  template<typename T3, typename K, typename Info>
  void build_remeshing_triangulation(const T3& tr,
                                     Remeshing_triangulation_3<K, Info>& remeshing_tr)
  {
    typedef typename T3::Triangulation_data_structure Tds;
    typedef Remeshing_triangulation_3<K, Info>::Tds   RTds;

    remeshing_tr.clear();

    remeshing_tr.set_infinite_vertex(
      remeshing_tr.tds().copy_tds(
        tr.tds(),
        tr.infinite_vertex(),
        internal::Vertex_converter<Tds, RTds>(),
        internal::Cell_converter<Tds, RTds>()));
  }

  template<typename T3, typename K, typename Info>
  void build_from_remeshing_triangulation(
    const Remeshing_triangulation_3<K, Info>& remeshing_tr,
    T3& tr)
  {
    typedef typename T3::Triangulation_data_structure Tds;
    typedef Remeshing_triangulation_3<K, Info>::Tds   RTds;

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
