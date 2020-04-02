// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_SURFACE_MESHER_TYPES_GENERATORS_H
#define CGAL_SURFACE_MESHER_TYPES_GENERATORS_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>

namespace CGAL {
  namespace Surface_mesher {
    namespace details {

      template <typename Base>
      class Triangulation_generator {
        typedef typename Base::Complex_2_in_triangulation_3 C2T3;
        typedef typename C2T3::Triangulation Triangulation;
      public:
        typedef Triangulation Type;
        typedef Type type;
      }; // end Triangulation_generator<Base>

      template <typename Base>
      class Facet_generator {
        typedef typename Triangulation_generator<Base>::type Tr;
      public:
        typedef typename Tr::Facet Type;
        typedef Type type;
      }; // end Facet_generator<Base>

      template <typename Base>
      class Edge_generator {
        typedef typename Triangulation_generator<Base>::type Tr;
      public:
        typedef typename Tr::Edge Type;
        typedef Type type;
      }; // end Edge_generator<Base>

      template <typename Base, typename Self, typename Element,
                typename PreviousLevel = Null_mesher_level>
      class Mesher_level_generator {
        typedef typename Base::Complex_2_in_triangulation_3 C2T3;
        typedef typename C2T3::Triangulation Triangulation;
        typedef Triangulation_mesher_level_traits_3<Triangulation> Tr_m_l_traits_3;
      public:
        typedef Mesher_level <Triangulation,
                              Self,
                              Element,
                              PreviousLevel,
                              Tr_m_l_traits_3> Type;
        typedef Type type;
      }; // end class Mesher_level_generator<Base, Self, Element>

    } // end namespace details
  } // end namespace Surface_mesher
} // end namespace CGAL

#endif // CGAL_SURFACE_MESHER_TYPES_GENERATORS_H
