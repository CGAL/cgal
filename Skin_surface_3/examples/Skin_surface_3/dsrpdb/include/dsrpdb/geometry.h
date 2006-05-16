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

#ifndef DSR_PDB_GEOMETRY_H
#define DSR_PDB_GEOMETRY_H

#include <dsrpdb/Protein.h>
#include <dsrpdb/iterator.h>
#include <vector>
#include <cassert>
#include <dsrpdb/geometry_internal.h>

namespace dsrpdb {
  //! Insert all the atom coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    dsrpdb::Residue::Atoms_iterator or dsrpdb::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void coordinates(Ait ab, Ait ae, Voit out){
    internal::filtered_coordinates(ab, ae, Yes(), out);
  }
  
  //! Insert all the backbone atom coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    dsrpdb::Residue::Atoms_iterator or dsrpdb::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void backbone_coordinates(Ait ab, Ait ae, Voit out){
    internal::filtered_coordinates(ab, ae, Is_backbone(), out);
  }
  
  //! Insert all the c_alpha coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    dsrpdb::Residue::Atoms_iterator or dsrpdb::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void ca_coordinates(Ait ab, Ait ae, Voit out){
    internal::filtered_coordinates(ab, ae, Is_CA(), out);
  }
  
  


  //! Insert all the atom coordinates from the protein in to the output iterator and an edge corresponding to each bond in the second output iterator
  /*!  The bonds differ from those returned by the bonds returned by
    the protein in that they are renumbered to correspond to the
    sequence positions of the atom coordinates.
  */
  template <class Voit, class Boit>
  inline void coordinates_and_bonds(const Protein &p, Voit out, Boit bout){
    internal::filtered_coordinates_and_bonds(p, Yes(), out, bout);
  }

//! Insert the backbone atom coordinates from the protein in to the output iterator and an edge corresponding to each bond in the second output iterator.
  /*!  The bonds differ from those returned by the bonds returned by
    the protein in that they are renumbered to correspond to the
    sequence positions of the atom coordinates.
  */
  template <class Voit, class Boit>
  inline void backbone_coordinates_and_bonds(const Protein &p, Voit out, Boit bout){
    internal::filtered_coordinates_and_bonds(p, Is_backbone(), out, bout);
  }


//! Insert the ca atom coordinates from the protein in to the output iterator and an edge corresponding to each bond in the second output iterator.
  /*!  The bonds differ from those returned by the bonds returned by
    the protein in that they are renumbered to correspond to the
    sequence positions of the atom coordinates.
  */
  template <class Voit, class Boit>
  inline void ca_coordinates_and_bonds(const Protein &p,Voit out, Boit bout){
    internal::filtered_coordinates_and_bonds(p, Is_CA(), out, bout);
  }


  //! Inserts the backbone atoms and the sidechain centroids and appropriate bonds into the iterators
  /*!  Note that this now returns all the backbone atoms and puts all
    the sidechains after allthe backbone atoms. The number of backbone
    atoms is returned. 
  */
  template <class Voit, class Boit>
  inline int simplified_coordinates_and_bonds(const Protein &p,Voit out, Boit bout){
    int index=-1;
    std::vector<int> ca_indices;
    //const Residue::Atom_label als[]={ Residue::AL_CA};
    for (Protein::Const_residues_iterator rit= p.residues_begin(); 
	 rit != p.residues_end(); ++rit){
      //for (unsigned int i=0; i< 3; ++i){
      if (rit->has_atom(Residue::AL_N)){
	*out= rit->atom(Residue::AL_N).cartesian_coords();
	++out;
	++index;
	if (index != 0) {
	  *bout = std::pair<int,int>(index-1, index);
	  ++bout;
	}
      }
      if (rit->has_atom(Residue::AL_CA)){
	*out= rit->atom(Residue::AL_CA).cartesian_coords();
	++out;
	
	++index;
	if (index != 0) {
	  *bout = std::pair<int,int>(index-1, index);
	  ++bout;
	}
	ca_indices.push_back(index);
      }
      if (rit->has_atom(Residue::AL_C)){
	*out= rit->atom(Residue::AL_C).cartesian_coords();
	++out;
	++index;
	if (index != 0) {
	  *bout = std::pair<int,int>(index-1, index);
	  ++bout;
	}
      }
    }
    int nbackbone=index+1;
    int rindex=0;
    for (Protein::Const_residues_iterator rit= p.residues_begin(); 
	 rit != p.residues_end(); ++rit){
      //for (unsigned int i=0; i< 3; ++i){
      if (rit->has_atom(Residue::AL_CA)){
	
	*out = rit->sidechain_point();
	++out;
	
	++index;

	*bout = std::pair<int,int>(ca_indices[rindex], index);
	++bout;
	
	++rindex;
      }
    }
    return nbackbone;
  }

 /*! Output a CGAL::Weighted_point_3 for each atom. It takes an
   instance of CGAL::Regular_triangulation_traits_3 as an argument to
   get the types from it.

   The weight is the squared radius.
 */
 template <class It, class K, class Voit>
 inline void all_weighted_points(It b, It e, K, Voit out){
   typedef typename K::FT                          Weight;
   typedef typename K::Point_3                     BP;
   typedef CGAL::Weighted_point<BP,Weight>         WP;
   for (; b != e; ++b){
     *out= WP(BP(b->second.cartesian_coords()),
	      b->second.radius()*b->second.radius());
     ++out;
   }
 }

/*! Output a CGAL::Weighted_point_3 for each atom. It takes an
   instance of CGAL::Regular_triangulation_traits_3 as an argument to
   get the types from it.

   The weight is the squared radius.
 */
  template <class K, class Voit>
  inline void all_weighted_points(const PDB &pdb, K k, Voit out){
    for (typename PDB::Const_models_iterator it= pdb.models_begin(); it != pdb.models_end(); ++it){
      for (typename Model::Const_chains_iterator it2= it->chains_begin();
	   it2 != it->chains_end(); ++it2){
	all_weighted_points(it2->atoms_begin(), it2->atoms_end(), k, out);
      }
    }
  }

}

#endif
