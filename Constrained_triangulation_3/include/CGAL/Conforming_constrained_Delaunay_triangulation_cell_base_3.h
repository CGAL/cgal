// Copyright (c) 2019-2024  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_data_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/SMDS_3/io_signature.h>

namespace CGAL {

/**
 * @ingroup PkgCT_3Classes
 * @brief Cell base class for the 3D conforming constrained Delaunay triangulation.
 *
 * This class is derived from the `Triangulation_cell_base_3` class and provides additional functionality
 * required by `make_conforming_constrained_Delaunay_triangulation_3()`.
 *
 * @tparam Traits The geometric traits class, which must be a model of `ConformingConstrainedDelaunayTriangulationTraits_3`.
 *         It should be the same as the geometric traits class of the triangulation.
 * @tparam Cell_base The base class for the cell, which must be a model of `TriangulationCellBase_3`.
 *
 * @cgalModels{ConformingConstrainedDelaunayTriangulationCellBase_3, SimplicialMeshCellBase_3, RemeshingCellBase_3}
 *
 * \note This cell base class also models the `SimplicialMeshCellBase_3` and `RemeshingCellBase_3` concepts, allowing the use of functionality from \ref Chapter_Tetrahedral_Remeshing "Tetrahedral Remeshing" and \ref Chapter_3D_Simplicial_Mesh_Data_Structure "3D Simplicial Mesh Data Structures", if the corresponding vertex base also models the right concepts.
 * \todo After discussion with Jane. Maybe there should be a second pair of Vb/Cb, designed to model the concepts of simplicial mesh and remeshing.
 *
 * \sa `CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3`
 */
template <typename Traits, typename Cell_base = Triangulation_cell_base_3<Traits> >
class Conforming_constrained_Delaunay_triangulation_cell_base_3
  : public Base_with_time_stamp<Cell_base>
{
  using Base = Base_with_time_stamp<Cell_base>;
  Conforming_constrained_Delaunay_triangulation_cell_data_3 ccdt_3_data_;

  mutable bool sliver_cache_validity_ = false;
  CDT_3_face_index subdomain_index_ = -1;
  double sliver_value_ = 0.;
public:
  // To get correct cell type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Cell_base::template Rebind_TDS<TDS3>::Other Cb3;
    typedef Conforming_constrained_Delaunay_triangulation_cell_base_3 <Traits, Cb3> Other;
  };

  // Constructors inherited from the base class
  using Base::Base;

  Conforming_constrained_Delaunay_triangulation_cell_data_3& ccdt_3_data() {
    return ccdt_3_data_;
  }

  const Conforming_constrained_Delaunay_triangulation_cell_data_3& ccdt_3_data() const {
    return ccdt_3_data_;
  }

  // model of SimplicialMeshCellBase_3
  using Surface_patch_index = CDT_3_face_index;
  using Subdomain_index = CDT_3_face_index;
  bool is_facet_on_surface(int i) const { return ccdt_3_data().is_facet_constrained(i); }
  Surface_patch_index surface_patch_index(int i) const { return ccdt_3_data().face_constraint_index(i) + 1; }
  void set_surface_patch_index(int i, Surface_patch_index index)
  {
    ccdt_3_data().set_face_constraint_index(i, index - 1);
  }
  Subdomain_index subdomain_index() const { return subdomain_index_; }
  void set_subdomain_index(Subdomain_index i) { subdomain_index_ = i; }

  // model of RemeshingCellBase_3
  void set_sliver_value(double value) {
    sliver_cache_validity_ = true;
    sliver_value_ = value;
  }
  double sliver_value() const {
    CGAL_assertion(is_cache_valid());
    return sliver_value_;
  }
  bool is_cache_valid() const { return sliver_cache_validity_; }
  void reset_cache_validity() const { sliver_cache_validity_ = false; }

  static std::string io_signature() {
    static_assert(
        std::is_same_v<
            decltype(std::declval<Conforming_constrained_Delaunay_triangulation_cell_data_3>().face_constraint_index(0)), int>);

    return Get_io_signature<Base>()() + "+(" + Get_io_signature<int>()() + ")[4]";
  }

  friend std::ostream&
  operator<<(std::ostream& os,
             const Conforming_constrained_Delaunay_triangulation_cell_base_3& c)
  {
    os << static_cast<const Base&>(c);
    for( unsigned li = 0; li < 4; ++li ) {
      if(IO::is_ascii(os)) {
        os << " " << c.ccdt_3_data().face_constraint_index(li);
      } else {
        CGAL::write(os, c.ccdt_3_data().face_constraint_index(li));
      }
    }
    return os;
  }
  friend std::istream&
  operator>>(std::istream& is,
             Conforming_constrained_Delaunay_triangulation_cell_base_3& c)
  {
    is >> static_cast<Base&>(c);
    if(!is) return is;
    for( int li = 0; li < 4; ++li ) {
      int i;
      if(IO::is_ascii(is)) {
        is >> i;
      } else {
        CGAL::read(is, i);
      }
      if(!is) return is;
      c.face_id[li] = i;
    }
    return is;
  }
};

} // namespace CGAL

#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_CELL_BASE_3_H
