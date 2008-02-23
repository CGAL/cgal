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

Model::Model(){}

void Model::swap_with(Model &o) {
  std::swap(extra_, o.extra_);
  std::swap(chains_, o.chains_);
  std::swap(heterogens_, o.heterogens_);
}
  

void Model::process_line(const char *line) {
  CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
  if (lt== CGAL_PDB_INTERNAL_NS::ATOM) {
    char chain=line[21];
    /*int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::atom_line_iformat_,
      &snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
      &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);
      assert(numscan >5);*/
    if (chains_.empty() || chains_.find(Chain_key(chain)) == chains_.end()){
      chains_[Chain_key(chain)]=Chain();
      //std::cout << "New chain " << chain << std::endl;
    }
    chains_[Chain_key(chain)].process_line(line);
  } else if (lt == CGAL_PDB_INTERNAL_NS::TER) {
    assert(!chains_.empty());
    //chains_.back().process_line(line);
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
    // What field is missing?
    int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::hetatom_line_iformat_,
			//"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
			&snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
			&x,&y,&z, &occupancy, &tempFactor, segID, element, charge);

    Atom::Index sindex(snum);
    Atom a;
    //a.set_label(Residue::atom_label(al));
    a.set_point(Point(x,y,z));
    a.set_index(sindex);
    if (numscan >10) {
      a.set_occupancy(occupancy);
    }
    if (numscan >11) {
      a.set_temperature_factor(tempFactor);
    }
    a.set_segment_id(segID);
    a.set_element(element);
    //a.set_charge(charge);
    a.set_type(Atom::string_to_type(name));
    
    std::string rname(resname);
    Heterogen_key hk(resname, resnum);
    if (heterogens_[hk].find(name) != heterogens_[hk].atoms_end()) {
      CGAL_LOG(Log::LOTS, "Duplicate atom in heterogen " << name << std::endl);
    } else {
      heterogens_[hk].insert(name, a);
    }
    heterogens_[hk].set_chain(chain);
    //residue(resnum-1)->set_coords (al, Point(x,y,z));
    //residue(resnum-1)->set_index(al, snum);
    //++cur_atom_index;

  } else if (lt== CGAL_PDB_INTERNAL_NS::ENDMDL){
      
  }
}

void Model::write(int model_index, std::ostream &out) const {
  char line[81];
  
  sprintf(line, "MODEL %8d         ", model_index);
  out << line << std::endl;
  int index=1;
  for (Chain_const_iterator it= chains_.begin(); it != chains_.end(); ++it){
    index= it->chain().write(it->key().index(), index, out);
  }
  for (unsigned int i=0; i< extra_.size(); ++i){
    out << extra_[i] << std::endl;
  }
  for (Heterogen_const_iterator it= heterogens_begin(); it != heterogens_end();
       ++it){
    index= it->heterogen().write(it->key().name(), it->key().number(),
                                 index, out);
  }
  out << "ENDMDL                       " << std::endl;
}

CGAL_PDB_END_NAMESPACE
