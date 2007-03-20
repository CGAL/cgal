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

#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/internal/pdb_utils.h>
#include <cassert>
#include <CGAL/PDB/internal/Error_logger.h>
CGAL_PDB_BEGIN_NAMESPACE

  PDB::PDB(std::istream &in, bool print_errors) {
    load(in, print_errors);
  }
  PDB::PDB(){}
  PDB::~PDB(){}
  
  void PDB::load(std::istream &in, bool print_errors){
    char line[600];
  
    CGAL_PDB_INTERNAL_NS::error_logger.set_is_output(print_errors);
    
    while (in.getline (line, 600)) {
      
      CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
      if (lt== CGAL_PDB_INTERNAL_NS::HEADER) {
	header_.push_back(std::string(line));
      } else if (lt== CGAL_PDB_INTERNAL_NS::DBREF){
	header_.push_back(std::string(line));
      } else if (lt== CGAL_PDB_INTERNAL_NS::SEQRES){
	header_.push_back(std::string(line));
      } else if (lt== CGAL_PDB_INTERNAL_NS::MODEL) {
	int mnum=0;
	char buf[81];
	sscanf(line, "%s %d", buf, &mnum);
	new_model(Model(mnum));
      } else if (  lt== CGAL_PDB_INTERNAL_NS::HETATM || lt== CGAL_PDB_INTERNAL_NS::ATOM 
		 || lt== CGAL_PDB_INTERNAL_NS::TER || lt== CGAL_PDB_INTERNAL_NS::ENDMDL){
	if (models_.empty()){
	  new_model(Model(0));
	}
	models_.back().process_line(line);
      } else if (lt== CGAL_PDB_INTERNAL_NS::MASTER){
      } else if (lt == CGAL_PDB_INTERNAL_NS::END){
      }
    }
    for (unsigned int i=0; i< models_.size(); ++i){
      for (unsigned int j=0; j< models_[i].number_of_chains(); ++j){
	models_[i].chain(j).set_has_bonds(true);
      }
    }
    
    CGAL_PDB_INTERNAL_NS::error_logger.dump();
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

CGAL_PDB_END_NAMESPACE
