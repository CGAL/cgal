// Copyright (c) 1997-2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_P2T2_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
#define CGAL_P2T2_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_face_base_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Delaunay_triangulation_on_lattice_2.h>

#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Triangulation_face_base_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Delaunay_triangulation_on_square_flat_torus_2.h>

#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>

#include <CGAL/Default.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace CGAL {

// Distinguish based on whether the domain is an iso_rectangle or a lattice (generic)

// Using CGAL::Default because the default is not the same for both specializations
template <class Gt_,
          class Tds_ = CGAL::Default,
          class Domain_ = typename Gt_::Domain>
class Periodic_2_Delaunay_triangulation_2
  : public Delaunay_triangulation_on_lattice_2<
             Gt_, typename Default::Get<Tds_,
                    Triangulation_data_structure_2<
                      Periodic_2_triangulation_vertex_base_2<Gt_>,
                      Triangulation_face_base_on_lattice_2<Gt_> > >::type>
{
private:
  typedef Delaunay_triangulation_on_lattice_2<
            Gt_, typename Default::Get<Tds_,
                   Triangulation_data_structure_2<
                     Periodic_2_triangulation_vertex_base_2<Gt_>,
                     Triangulation_face_base_on_lattice_2<Gt_> > >::type> Base;

public:
  template <typename ...Args>
  Periodic_2_Delaunay_triangulation_2(const Args& ...args) : Base(args...) { }

  template <typename ...Args>
  Periodic_2_Delaunay_triangulation_2(const Args&& ...args) : Base(std::forward<Args>(args)...) { }
};

template <class Gt_, class Tds_, class K_>
class Periodic_2_Delaunay_triangulation_2<Gt_, Tds_, CGAL::Iso_rectangle_2<K_> >
  : public Delaunay_triangulation_on_square_flat_torus_2<
             Gt_, typename Default::Get<Tds_,
                    Triangulation_data_structure_2<
                      Periodic_2_triangulation_vertex_base_2<Gt_>,
                      Triangulation_face_base_on_square_flat_torus_2<Gt_> > >::type>
{
private:
  typedef Delaunay_triangulation_on_square_flat_torus_2<
            Gt_, typename Default::Get<Tds_,
                   Triangulation_data_structure_2<
                     Periodic_2_triangulation_vertex_base_2<Gt_>,
                     Triangulation_face_base_on_square_flat_torus_2<Gt_> > >::type> Base;

public:
  template <typename ...Args>
  Periodic_2_Delaunay_triangulation_2(const Args& ...args) : Base(args...) { }

  template <typename ...Args>
  Periodic_2_Delaunay_triangulation_2(const Args&& ...args) : Base(std::forward<Args>(args)...) { }
};

} // namespace CGAL

#endif // CGAL_P2T2_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
