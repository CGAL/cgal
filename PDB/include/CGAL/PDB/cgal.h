/* Copyright 2004
   Stanford University

   This file is part of the DSR PDB Library.

   The DSR PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. */

#ifndef CGAL_DSR_PDB_CGAL_H
#define CGAL_DSR_PDB_CGAL_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/iterator.h>
#include <vector>
#include <cassert>

CGAL_PDB_BEGIN_NAMESPACE

 /*! Output a CGAL::Weighted_point_3 for each atom. It takes an
    instance of CGAL::Regular_triangulation_traits_3 as an argument to
    get the types from it.

    The weight is the squared radius.
  */
  template <class It, class K, class Voit>
  inline void all_weighted_points(It b, It e, K, Voit out){
    typedef typename K::Weighted_point WP;
    typedef typename K::Bare_point BP;
    for (; b != e; ++b){
      *out= WP(BP(b->second.cartesian_coords().x(), b->second.cartesian_coords().y(), b->second.cartesian_coords().z()), 
	       b->second.radius()*b->second.radius());
      ++out;
    }
  }

  /*! Output a CGAL::Weighted_point_3 for each atom. It takes an
    instance of CGAL::Regular_triangulation_traits_3 as an argument to
    get the types from it. The points for HETATM records are also
    returned.

    The weight is the squared radius.
  */
  template <class K, class Voit>
  inline void all_weighted_points(const Model &m, K k, Voit out){
    for (typename Model::Const_chains_iterator it2= m.chains_begin();
	 it2 != m.chains_end(); ++it2){
      all_weighted_points(it2->atoms_begin(), it2->atoms_end(), k, out);
    }
    all_weighted_points(m.hetatoms_begin(), m.hetatoms_end(), k, out);
    
  }



 /*! Output a CGAL::Weighted_point_3 for each atom. It takes an
    instance of CGAL::Regular_triangulation_traits_3 as an argument to
    get the types from it.

    The weight is the squared radius.
  */
  template < class K, class It,  class Voit>
  inline void all_spheres(It b, It e, Voit out){
    typedef typename K::Sphere_3 S;
    typedef typename K::Point_3 P;
    for (; b != e; ++b){
      *out= S(P(b->second.cartesian_coords().x(),
		b->second.cartesian_coords().y(), 
		b->second.cartesian_coords().z()), 
	       b->second.radius()*b->second.radius());
      ++out;
    }
  }

  /*! Output a CGAL::Sphere_3 for each atom. The points for HETATM
    records are also returned.

    The weight is the squared radius.
  */
  template <class K, class Voit>
  inline void all_spheres(const Model &m,Voit out){
    for (typename Model::Const_chains_iterator it2= m.chains_begin();
	 it2 != m.chains_end(); ++it2){
      all_spheres<K>(it2->atoms_begin(), it2->atoms_end(), out);
    }
    all_spheres<K>(m.hetatoms_begin(), m.hetatoms_end(), out);
  }
  
CGAL_PDB_END_NAMESPACE
#endif
