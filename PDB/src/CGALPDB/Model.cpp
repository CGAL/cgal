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
#include <CGAL/PDB/internal/pdb_utils.h>
#include <CGAL/PDB/internal/Error_logger.h>

#include <boost/format.hpp>

#include <cassert>
#include <cstdio>

using std::sscanf;

CGAL_PDB_BEGIN_NAMESPACE

/*static unsigned int  getSequenceNumber (const char* line)
{

  char s[5];

  strncpy (s, line + 22, 4); s[4] = '\0';
  return atoi (s);

  }*/



Model::Model(){}

void Model::swap_with(Model &o) {
  std::swap(extra_, o.extra_);
  std::swap(chains_, o.chains_);
  std::swap(heterogens_, o.heterogens_);
}
 

void Model::process_hetatom(const char *line) {
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
  
  add_hetatom(numscan, snum, name, alt, resname, chain, resnum,
                insertion_residue_code,
                x,y,z, occupancy, tempFactor, segID, element, charge);
 
}

void Model::add_hetatom(int numscan, int snum, char name[], char /* alt */, 
                        char resname[], char chain, int resnum,
                        char /* insertion_residue_code */, 
                        float x, float y, float z,
                        float occupancy, float tempFactor,
                        char segID[], char element[], char /* charge */ []) {
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
}



void Model::process_atom(const char *line) {

  // the read values are not zero padded so we must fill the buffers for strings with 0s.
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
  int numscan= sscanf(line, CGAL_PDB_INTERNAL_NS::atom_line_iformat_,
                      //"ATOM  %5d%4s%1c%3s%1c%4d%1c%8f%8f%8f%6f%6f%4s%2s%2s",
                      &snum, name, &alt, resname, &chain, &resnum,
                      &insertion_residue_code,
                      &x,&y,&z, &occupancy, &tempFactor, segID, element, charge);

  // recognize non-standard residues posing as standard ones and process them as heterogens
  if (Monomer::type(resname) == Monomer::INV) {
    CGAL_LOG(Log::LOTS, "Treating unknown ATOM residue type as heterogen: "
             << line << std::endl);
    add_hetatom(numscan,  snum, name, alt, resname, chain, resnum,
                insertion_residue_code,
                x,y,z, occupancy, tempFactor, segID, element, charge);
    return;
  }

  if (chains_.empty() || chains_.find(Chain_key(chain)) == chains_.end()){
    chains_[Chain_key(chain)]=Chain();
    //std::cout << "New chain " << chain << std::endl;
  }
  Chain &curchain= chains_[Chain_key(chain)];
  /*if (chain_==' ') chain_=chain;
    if (!(chain_==' ' || chain== chain_)){
    std::ostringstream oss;
    oss << "Confusion over chain numbers. Expected " << chain_ << " got " << chain 
    << " on line:\n" << line;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return;
    }*/

  if (resnum < 0) {
    std::ostringstream oss;
    oss << "Got negative residue index on line " << line;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return;
  }


  Chain::Monomer_key resindex(resnum);
      
  Monomer::Atom_key al= Monomer::atom_key(name);
      
  if (al != Monomer::AL_OTHER) {
    Monomer *cur_residue=NULL;
  
    if (insertion_residue_code != ' '){
      if (curchain.insert_residues_.find(resindex) 
          == curchain.insert_residues_.end()
          || curchain.insert_residues_[resindex].find(Chain::IR_key(insertion_residue_code)) 
          == curchain.insert_residues_[resindex].end()) {
        Monomer::Type rl= Monomer::type(resname);
        curchain.insert_residues_[resindex][Chain::IR_key(insertion_residue_code)]=Monomer(rl); 
      } 
      cur_residue= &curchain.insert_residues_[resindex][Chain::IR_key(insertion_residue_code)];
    } else {
      if (curchain.residues_.find(resindex) == curchain.residues_.end()) {
        std::string nm(line, 17, 3);
        Monomer::Type rl= Monomer::type(resname);
        Monomer m(rl);
        Chain::Monomer_iterator mit= curchain.residues_.insert(Chain::Container::value_type(resindex, m));
        cur_residue= &mit->monomer();
        CGAL_assertion(mit->key()== resindex);
        //residues_[resindex]=m;
        CGAL_postcondition(m.type()== cur_residue->type());
      } else {
        CGAL_assertion(curchain.residues_.find(resindex) < curchain.residues_.end());
        CGAL_assertion(curchain.residues_.find(resindex) >= curchain.residues_.begin());
        cur_residue =&curchain.residues_.find(resindex)->monomer();
        CGAL_assertion(curchain.residues_.find(resindex)->key()== resindex);
      }
    }    
      
    Atom a;
    //a.set_label(Residue::atom_label(al));
    a.set_point(Point(x,y,z));

    if (numscan >10) {
      a.set_occupancy(occupancy);
    }
    if (numscan >11) {
      a.set_temperature_factor(tempFactor);
    }
    a.set_segment_id(segID);
    a.set_element(element);
    a.set_charge(charge);
      
    /*if (cur_residue->index().index() 
      != static_cast<unsigned int>(resnum)){
      std::ostringstream oss;
      oss << "Confusion over residue numbers. Expected" << cur_residue->index()
      << " got " << resnum 
      << " on line:\n" << line << std::endl;
      CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error(oss.str().c_str());
      return;
      }*/
    if (cur_residue->type() != Monomer::type(resname)){
      std::ostringstream oss;
      oss << "Confusion over residue types. Expected " 
          << Monomer::type_string(cur_residue->type())
          << " got " << Monomer::type_string(Monomer::type(resname))
          << " on line:\n" << line << std::endl;
      CGAL_PDB_INTERNAL_NS::error_logger.new_fatal_error(oss.str().c_str());
    }
    //assert(cur_residue->type() == Residue::type(resname));
    Monomer::Atom_iterator rit= cur_residue->insert(al, a);
    if (rit == cur_residue->atoms_end()) {
      std::ostringstream oss;
      oss << "Error adding atom to residue " <<  resindex << std::endl;
      CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    }
    
    //residue(resnum-1)->set_coords (al, Point(x,y,z));
    //residue(resnum-1)->set_index(al, snum);
    //++cur_atom_index;
  } else {
    /*std::ostringstream out;
      out << "Unhandled atom type: " << name;
      CGAL_PDB_INTERNAL_NS::error_logger.new_warning(out.str().c_str());*/
  }
 
}

void Model::process_line(const char *line) {
  CGAL_PDB_INTERNAL_NS::Line_type lt= CGAL_PDB_INTERNAL_NS::line_type(line);
  if (lt== CGAL_PDB_INTERNAL_NS::ATOM) {
    process_atom(line);
  } else if (lt == CGAL_PDB_INTERNAL_NS::TER) {
    //CGAL_assertion(!chains_.empty());
    //chains_.back().process_line(line);
  } else if (lt== CGAL_PDB_INTERNAL_NS::HETATM){
    process_hetatom(line);

  } else if (lt== CGAL_PDB_INTERNAL_NS::ENDMDL){
      
  }
}

void Model::write(int model_index, std::ostream &out) const {
  
  out << boost::format("MODEL %8d         ") % model_index << std::endl;
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
