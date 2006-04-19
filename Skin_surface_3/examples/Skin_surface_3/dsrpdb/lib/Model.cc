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

#include "dsrpdb/Model.h"
#include <cassert>
#include "pdb_utils.h"


namespace dsrpdb {
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
    dsrpdb_internal::Line_type lt= dsrpdb_internal::line_type(line);
    if (lt== dsrpdb_internal::ATOM) {

      int snum=-1;
      char name[5]={'\0'};
      char alt='\0';
      char resname[4]={'\0'};
      char chain;
      int resnum=-1;      char insertion_residue_code;
      float x,y,z;
      float occupancy, tempFactor;
      char segID[5]={'\0'}, element[3]={'\0'}, charge[3]={'\0'};
      int numscan= sscanf(line, dsrpdb_internal::atom_line_iformat_,
			  //"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
			  &snum, name, &alt, resname, &chain, &resnum, &insertion_residue_code,
			  &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);
      assert(numscan >5);
      if (chains_.empty() || chains_.back().chain() != chain){
	chains_.push_back(Protein());
      }
      chains_.back().process_line(line);
    } else if (lt == dsrpdb_internal::TER) {
      assert(!chains_.empty());
      chains_.back().process_line(line);
    } else if (lt== dsrpdb_internal::HETATM){
      extra_.push_back(line);
    } else if (lt== dsrpdb_internal::ENDMDL){
      
    }
  }

  void Model::write(std::ostream &out) const {
    char line[81];
  
    sprintf(line, "MODEL %8d         ", index_);
    out << line << std::endl;
    for (unsigned int i=0; i< chains_.size(); ++i){
      chains_[i].write(out);
    }
    for (unsigned int i=0; i< extra_.size(); ++i){
      out << extra_[i] << std::endl;
    }
    out << "ENDMDL                       " << std::endl;
  }
};
