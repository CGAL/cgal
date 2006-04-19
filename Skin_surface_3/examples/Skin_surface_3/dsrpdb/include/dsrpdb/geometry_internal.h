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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_PDB_GEOMETRY_INTERNAL_H
#define DSR_PDB_GEOMETRY_INTERNAL_H

#include <dsrpdb/Protein.h>

namespace dsrpdb {
  namespace internal {

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
    void filtered_coordinates_and_bonds(const Protein &p, Filt f, Voit out, Boit bout){
      std::vector<int> map(p.number_of_atoms(), -1);
      int seen=0;
      for (typename Protein::Const_atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
	unsigned int index= it->second.index();
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
      
      for (typename Protein::Bonds_iterator it= p.bonds_begin(); it != p.bonds_end(); ++it){
	assert(it->first && it->second);
	unsigned int a= it->first;
	unsigned int b= it->second;
	assert(a < map.size() && b < map.size());
	if (map[a] !=-1 && map[b]!=-1) { 
	  *bout = std::pair<int,int>(map[a], map[b]);
	  ++bout;
	}
      }
    }
  } 
}

#endif
