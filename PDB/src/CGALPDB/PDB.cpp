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
  std::map<char, std::string> names;
  std::string last_name;
  Model_key cur_model(0);
  while (in.getline (line, 600)) {
      
    CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
    if (lt== CGAL_PDB_INTERNAL_NS::HEADER) {
      header_.push_back(std::string(line));
    } else if (lt == CGAL_PDB_INTERNAL_NS::COMPND) {
      header_.push_back(std::string(line));
      int ti;
      char chain;
      if (sscanf(line, "COMPND %d CHAIN: %c", &ti, &chain) == 2) {
	names[chain]=last_name;
      } else if (sscanf(line, "COMPND %d MOLECULE:", &ti) ==1) {
	last_name= std::string(line+20);
      }
    } else if (lt== CGAL_PDB_INTERNAL_NS::DBREF){
      header_.push_back(std::string(line));
    } else if (lt== CGAL_PDB_INTERNAL_NS::SEQRES){
      header_.push_back(std::string(line));
    } else if (lt== CGAL_PDB_INTERNAL_NS::MODEL) {
      int mnum=0;
      char buf[81];
      sscanf(line, "%s %d", buf, &mnum);
      //new_model(Model_key(mnum), Model());
      cur_model= Model_key(mnum);
    } else if (  lt== CGAL_PDB_INTERNAL_NS::HETATM 
		 || lt== CGAL_PDB_INTERNAL_NS::ATOM 
		 || lt== CGAL_PDB_INTERNAL_NS::TER 
		 || lt== CGAL_PDB_INTERNAL_NS::ENDMDL){
      models_[cur_model].process_line(line);
    } else if (lt== CGAL_PDB_INTERNAL_NS::MASTER){
    } else if (lt == CGAL_PDB_INTERNAL_NS::END){
    }
  }

  for (Model_iterator it= models_begin(); it != models_end(); ++it){
    for (Model::Chain_iterator cit = it->model().chains_begin(); cit != it->model().chains_end(); ++cit){
      cit->chain().set_has_bonds(true);
      if (names.find(cit->key().to_index()) != names.end()) {
	cit->chain().set_name(names[cit->key().to_index()]);
      }
    }
  }
    
  CGAL_PDB_INTERNAL_NS::error_logger.dump();
}


void PDB::swap_with(PDB &o) {
  swap(header_, o.header_);
  swap(models_, o.models_);
}


std::ostream& PDB::write(std::ostream &out) const {
  for (unsigned int i=0; i< header_.size(); ++i){
    out << header_[i] << std::endl;
  }
  for (Model_const_iterator it = models_begin(); it != models_end(); ++it){
    it->model().write(it->key().to_index(), out);
  }
  out << "END   \n";
  return out;
}


PDB::Model_key PDB::push_back(const Model &m) {
  Model_key k(0);
  if (!empty()) k= Model_key((--models_end())->key().to_index()+1);
  insert(k, m);
  return k;
}

CGAL_PDB_END_NAMESPACE
