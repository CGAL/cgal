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

#ifndef CGAL_DSR_PDB_GEOMETRY_INTERNAL_H
#define CGAL_DSR_PDB_GEOMETRY_INTERNAL_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Protein.h>
CGAL_PDB_BEGIN_INTERNAL_NAMESPACE


template <class Filt, class Ait, class Voit>
void filtered_coordinates(Ait ab, Ait ae, Filt f, Voit out){
  for (; ab != ae; ++ab){
    if (f(*ab)) {
      *out = ab->second.cartesian_coords();
      ++out;
    }
  }
}

template <class Filt, class Voit, class Boit>
void filtered_coordinates_and_bonds(const CGAL_PDB_NS::Protein &p, Filt f, Voit out, Boit bout){
  int na= p.number_of_atoms();
  std::vector<int> map(na, -1);
  int seen=0;
  for (typename CGAL_PDB_NS::Protein::Const_atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
    unsigned int index= it->second.index().to_index();
    if (index >= map.size()){
      map.resize(index+1, -1);
    }
    if (f(*it)) {
      map[index]= seen;
      ++seen;
      *out = it->second.cartesian_coords();
      ++out;
    }
  }
      
  for (typename CGAL_PDB_NS::Protein::Bonds_iterator it= p.bonds_begin(); it != p.bonds_end(); ++it){
    assert(it->first != CGAL_PDB_NS::Atom::Index()
	   && it->second != CGAL_PDB_NS::Atom::Index());
    unsigned int a= it->first.to_index();
    unsigned int b= it->second.to_index();
    assert(a < map.size() && b < map.size());
    if (map[a] !=-1 && map[b]!=-1) { 
      *bout = std::pair<int,int>(map[a], map[b]);
      ++bout;
    }
  }
}
CGAL_PDB_END_INTERNAL_NAMESPACE
#endif
