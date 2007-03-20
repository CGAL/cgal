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

#ifndef CGAL_DSR_PDB_GEOMETRY_H
#define CGAL_DSR_PDB_GEOMETRY_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/iterator.h>
#include <vector>
#include <cassert>
#include <CGAL/PDB/geometry_internal.h>
CGAL_PDB_BEGIN_NAMESPACE
  //! Insert all the atom coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    CGAL::PDB::Residue::Atoms_iterator or CGAL::PDB::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void coordinates(Ait ab, Ait ae, Voit out){
    CGAL_PDB_INTERNAL_NS::filtered_coordinates(ab, ae, Yes(), out);
  }
  
  //! Insert all the backbone atom coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    CGAL::PDB::Residue::Atoms_iterator or CGAL::PDB::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void backbone_coordinates(Ait ab, Ait ae, Voit out){
    CGAL_PDB_INTERNAL_NS::filtered_coordinates(ab, ae, Is_backbone(), out);
  }
  
  //! Insert all the c_alpha coordinates from the list in to the output iterator.
  /*!  See the file check_protein.cc for an example of usage. Either
    CGAL::PDB::Residue::Atoms_iterator or CGAL::PDB::Protein::Atoms_iterator
    are appropriate input.
  */
  template <class Ait, class Voit>
  inline void ca_coordinates(Ait ab, Ait ae, Voit out){
    CGAL_PDB_INTERNAL_NS::filtered_coordinates(ab, ae, Is_CA(), out);
  }



  /*!  Insert all the atom coordinates from the protein in to the
    output iterator and an edge corresponding to each bond in the
    second output iterator.  The bonds differ from those returned by
    the bonds returned by the protein in that they are renumbered to
    correspond to the sequence positions of the atom coordinates.
  */
  template <class Voit, class Boit>
  inline void coordinates_and_bonds(const Protein &p,
				    Voit out, Boit bout){
    CGAL_PDB_INTERNAL_NS::filtered_coordinates_and_bonds(p, Yes(), out, bout);
  }



  /*! Insert the backbone atom coordinates from the protein in to the
    output iterator and an edge corresponding to each bond in the
    second output iterator. The bonds differ from those returned by
    the bonds returned by the protein in that they are renumbered to
    correspond to the sequence positions of the atom coordinates.
  */
  template <class Voit, class Boit>
  inline void backbone_coordinates_and_bonds(const Protein &p, 
					     Voit out, Boit bout){
    CGAL_PDB_INTERNAL_NS::filtered_coordinates_and_bonds(p, Is_backbone(), out, bout);
  }
 


  //! Inserts the backbone atoms and the sidechain centroids and appropriate bonds into the iterators
  /*!  Note that this now returns all the backbone atoms and puts all
    the sidechains after allthe backbone atoms. The number of backbone
    atoms is returned. 
  */
  template <class Voit, class Boit>
  inline int simplified_coordinates_and_bonds(const Protein &p,
					      Voit out, Boit bout){
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
	
	Point spt= rit->sidechain_point();
	if (spt != rit->atom(Residue::AL_CA).cartesian_coords()) {
	  *out = rit->sidechain_point();
	  ++out;
	  
	  ++index;
	  
	  *bout = std::pair<int,int>(ca_indices[rindex], index);
	  ++bout;
	}
	++rindex;
      }
    }
    return nbackbone;
  }

CGAL_PDB_END_NAMESPACE
#endif
