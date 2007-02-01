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

#include <dsrpdb/PDB.h>
#include "pdb_utils.h"
#include <cassert>
#include <dsrpdb_internal/Error_logger.h>

namespace dsrpdb {
 
  PDB::PDB(std::istream &in, bool print_errors) {
    load(in, print_errors);
  }
  PDB::PDB(){}
  PDB::~PDB(){}
  
  void PDB::load(std::istream &in, bool print_errors){
    char line[600];
  
    dsrpdb_internal::error_logger.set_is_output(print_errors);
    
    while (in.getline (line, 600)) {
      
      dsrpdb_internal::Line_type lt= dsrpdb_internal::line_type(line);
      if (lt== dsrpdb_internal::HEADER) {
	header_.push_back(std::string(line));
      } else if (lt== dsrpdb_internal::DBREF){
	header_.push_back(std::string(line));
      } else if (lt== dsrpdb_internal::SEQRES){
	header_.push_back(std::string(line));
      } else if (lt== dsrpdb_internal::MODEL) {
	int mnum=0;
	char buf[81];
	sscanf(line, "%s %d", buf, &mnum);
	new_model(Model(mnum));
      } else if (lt== dsrpdb_internal::HETATM || lt== dsrpdb_internal::ATOM 
		 || lt== dsrpdb_internal::TER || lt== dsrpdb_internal::ENDMDL){
	if (models_.empty()){
	  new_model(Model(0));
	}
	models_.back().process_line(line);
      } else if (lt== dsrpdb_internal::MASTER){
      } else if (lt == dsrpdb_internal::END){
      }
    }
    for (unsigned int i=0; i< models_.size(); ++i){
      for (unsigned int j=0; j< models_[i].number_of_chains(); ++j){
	models_[i].chain(j).set_has_bonds(true);
      }
    }
    
    dsrpdb_internal::error_logger.dump();
  }

  
  void PDB::new_model(const Model &m){
    models_.push_back(m);
  }
  Model &PDB::model(unsigned int i) {
    assert(i < models_.size());
    return models_[i];
  };

  const Model &PDB::model(unsigned int i) const {
    assert(i < models_.size());
    return models_[i];
  };
  size_t PDB::number_of_models() const {
    return models_.size(); 
  }

  void PDB::write(std::ostream &out) const {
    for (unsigned int i=0; i< header_.size(); ++i){
      out << header_[i] << std::endl;
    }
    for (unsigned int i=0; i< models_.size(); ++i){
      models_[i].write(out);
    }
    out << "END   \n";
  }


  PDB::Header_iterator PDB::header_begin() const {
    return header_.begin();
  }
  PDB::Header_iterator PDB::header_end() const {
    return header_.end();
  }
};
