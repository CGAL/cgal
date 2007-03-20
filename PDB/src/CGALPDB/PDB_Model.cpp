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

#include <CGAL/PDB/Model.h>
#include <cassert>
#include <CGAL/PDB/internal/pdb_utils.h>
CGAL_PDB_BEGIN_NAMESPACE

  Model::Model():index_(-1){}
  Model::Model(unsigned int i): index_(i){}
  size_t Model::number_of_chains() const {
    return chains_.size();
  }
  Protein &Model::chain(unsigned int i) {
    assert(i < chains_.size());
    return chains_[i];
  }
  const Protein &Model::chain(unsigned int i) const {
    assert(i < chains_.size());
    return chains_[i];
  }
  void Model::new_chain(const Protein &p){
    chains_.push_back(p);
  }
  

  void Model::process_line(const char *line) {
    CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
    if (lt== CGAL_PDB_INTERNAL_NS::ATOM) {
      char chain=line[21];
      /*int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::atom_line_iformat_,
			  &snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
			  &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);
			  assert(numscan >5);*/
      if (chains_.empty() || chains_.back().chain() != chain){
	chains_.push_back(Protein());
	//std::cout << "New chain " << chain << std::endl;
      }
      chains_.back().process_line(line);
    } else if (lt == CGAL_PDB_INTERNAL_NS::TER) {
      assert(!chains_.empty());
      chains_.back().process_line(line);
    } else if (lt== CGAL_PDB_INTERNAL_NS::HETATM){
      int snum=-1;
      char name[5]={'\0','\0','\0','\0','\0',};
      char alt='\0';
      char resname[4]={'\0','\0','\0','\0'};
      char chain;
      int resnum=-1;
      char insertion_residue_code;
      float x,y,z;
      float occupancy, tempFactor;
      char segID[5]={'\0','\0','\0','\0','\0'};
      char element[3]={'\0','\0','\0'};
      char charge[3]={'\0','\0','\0'};
      // What is a field is missing?
      int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::hetatom_line_iformat_,
			  //"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
			  &snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
			  &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);

      Atom::Index sindex(snum);
      Atom a;
      //a.set_label(Residue::atom_label(al));
      a.set_cartesian_coords(Point(x,y,z));
      a.set_index(sindex);
      if (numscan >10) {
	a.set_occupancy(occupancy);
      }
      if (numscan >11) {
	a.set_temperature_factor(tempFactor);
      }
      a.set_segment_id(segID);
      a.set_element(element);
      a.set_charge(charge);
      a.set_type(Atom::string_to_type(name));
      
      hetatoms_.push_back(std::pair<Hetatom_data, Atom>(Hetatom_data(resname, name, resnum, chain), a));
	
	//residue(resnum-1)->set_coords (al, Point(x,y,z));
	//residue(resnum-1)->set_index(al, snum);
	//++cur_atom_index;

    } else if (lt== CGAL_PDB_INTERNAL_NS::ENDMDL){
      
    }
  }

  void Model::write(std::ostream &out) const {
    char line[81];
  
    sprintf(line, "MODEL %8d         ", index_);
    out << line << std::endl;
    for (unsigned int i=0; i< chains_.size(); ++i){
      chains_[i].write(out);
    }
    for (unsigned int i=0; i< hetatoms_.size(); ++i){
      //Point pt= res->cartesian_coords(al);
      const Atom &a= hetatoms_[i].second;
      Point pt = a.cartesian_coords();
      char alt=' ';
      char insertion_residue_code=' ';
      sprintf(line, CGAL_PDB_INTERNAL_NS::hetatom_line_oformat_,
	      a.index().to_index(), 
	      hetatoms_[i].first.atom_name(), alt,
	      hetatoms_[i].first.molecule_name(), 
	      hetatoms_[i].first.chain(), 
	      static_cast<unsigned int>(hetatoms_[i].first.molecule_number()), 
	      insertion_residue_code,
	      pt.x(), pt.y(), pt.z(), 
	      a.occupancy(), a.temperature_factor(), a.segment_id(),
	      a.element(), a.charge());
      out << line << std::endl;
    }
    for (unsigned int i=0; i< extra_.size(); ++i){
      out << extra_[i] << std::endl;
    }
    out << "ENDMDL                       " << std::endl;
  }

CGAL_PDB_END_NAMESPACE
