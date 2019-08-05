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


#ifndef CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_H
#define CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_H

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  namespace internal
  {
    class Fake_MD
    {
    public:
      typedef int Index;
      typedef int Surface_patch_index;
      typedef int Subdomain_index;
    };
  }

  template<typename GT,
           typename Vb = CGAL::Triangulation_vertex_base_3<GT> >
  class Remeshing_vertex_base
    : public CGAL::Mesh_vertex_base_3<GT, internal::Fake_MD, Vb>
  {
    typedef CGAL::Mesh_vertex_base_3<GT, internal::Fake_MD, Vb> Base;

  private:
    short dimension_;
    std::size_t time_stamp_;

  public:
    Remeshing_vertex_base() : dimension_(-1)
					   // time_stamp_ // do not initialize
    {}
    typedef int                   Index;

    // To get correct vertex type in TDS
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
      typedef Remeshing_vertex_base<GT, Vb3> Other;
    };

    // Returns the dimension of the lowest dimensional face of the input 3D
    // complex that contains the vertex
    int in_dimension() const {
      if (dimension_ < -1) return -2 - dimension_;
      else return dimension_;
    }

    // Sets the dimension of the lowest dimensional face of the input 3D complex
    // that contains the vertex
    void set_dimension(const int dimension) {
      CGAL_assertion(dimension < 4);
      dimension_ = short(dimension);
    }

    /// For the determinism of Compact_container iterators
    ///@{
    typedef Tag_true Has_timestamp;
    std::size_t time_stamp() const {
      return time_stamp_;
    }
    void set_time_stamp(const std::size_t& ts) {
      time_stamp_ = ts;
    }
    ///@}

  };

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TET_ADAPTIVE_REMESHING_VERTEX_BASE_H
