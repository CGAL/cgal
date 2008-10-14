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

#include <CGAL/PDB/Monomer.h>
#include <CGAL/PDB/internal/Monomer_data.h>
#include <CGAL/PDB/internal/pdb_utils.h>
#include <CGAL/PDB/internal/Error_logger.h>

#include <boost/format.hpp>

#include <iostream>
#include <limits>
#include <sstream>
#include <cstdio>

using std::sscanf;

CGAL_PDB_BEGIN_NAMESPACE

void Monomer::set_has_bonds(bool tf) {
  typedef Monomer_data::Possible_bond Possible_bond;

  if (!tf) {
    bonds_.clear();
  } else {
    const std::vector<Possible_bond >& bonds= Monomer_data::amino_acid_data_[type()].bonds;
    for (unsigned int i=0; i< bonds.size(); ++i){
      Atom_const_iterator itf= find(bonds[i].first);
      Atom_const_iterator its= find(bonds[i].second);
      if ( itf != atoms_.end() && its != atoms_.end()){
	bonds_.push_back(Bond(Bond_endpoint(itf),
			      Bond_endpoint(its)));
      } 
    }
  }
}










bool Monomer::can_have_atom(Monomer::Atom_key ial) const {
  if (ial== AL_INVALID) return false; // hack to avoid some circularity
  Monomer::Atom_key al= Monomer_data::fix_atom_label(label_, ial);
  //assert(Monomer_data::amino_acid_data_.find(label_)!= Monomer_data::amino_acid_data_.end());
  for (unsigned int i=0; i< Monomer_data::amino_acid_data_[label_].atoms.size(); ++i){
    if (Monomer_data::amino_acid_data_[label_].atoms[i] == al) return true;
  }
  return false;
}



void Monomer::copy_from(const Monomer &o) {
  atoms_= o.atoms_;
  label_= o.label_;
  set_has_bonds(false);
  set_has_bonds(o.has_bonds());
}
    


void Monomer::swap_with(Monomer &o) {
  std::swap(atoms_, o.atoms_);
  std::swap(bonds_, o.bonds_);
  std::swap(label_, o.label_);
}



void Monomer::erase_atom(Monomer::Atom_key al) {
  Monomer::Atom_key fal= Monomer_data::fix_atom_label(label_, al);
  //int fa= find_atom(al);
  atoms_.erase(fal);
}




Monomer::Atom_iterator Monomer::insert_internal(Atom_key ial, const Atom &a) {
  Monomer::Atom_key al= Monomer_data::fix_atom_label(label_, ial);
  if (!can_have_atom(al)){
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning((std::string("Trying to set invalid atom ")
						    + atom_key_string(ial) 
						    + " on a residue of type "
						    + type_string(label_)).c_str());
    return atoms_end();
  }
  if (al == AL_INVALID) {
    return atoms_end();
  }
  if (atoms_.find(al) != atoms_.end()) {
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning((std::string("Duplicate atoms ")
						    + atom_key_string(ial) 
						    + " on a residue of type "
						    + type_string(label_)).c_str());
    return atoms_end();
  }
  //assert(atoms_.find(al)== atoms_.end());
  //bonds_.clear();

  Atom_iterator ret= atoms_.insert(Atoms::value_type(al,a));
  // This is a bit evil, I haven't figured out a better way though.
  ret->atom().set_type(element(al));
  if (has_bonds()){
    set_has_bonds(false);
    set_has_bonds(true);
  }
  //atoms_.push_back(a);
  return ret;
}





Monomer::Atom_key Monomer::fix_atom_key(Monomer::Atom_key ial) const {
  return Monomer_data::fix_atom_label(label_, ial);
}





void Monomer::dump(std::ostream &out) const {
  out << "Type: " << type_string(type()) << std::endl;
  const std::vector<Atom_key> &valid_atoms= Monomer_data::amino_acid_data_[type()].atoms;
  for (unsigned int i=0; i< valid_atoms.size(); ++i){
    Monomer::Atom_key al= valid_atoms[i];
    out << Monomer::atom_key_string(al); //Monomer::write_atom_label(al, out);
    if (find(al) != atoms_end()){
      out << " (" << find(al)->atom().point().x() << ", " 
	  << find(al)->atom().point().y() << ", "
	  << find(al)->atom().point().z() << ") " << std::endl;
    } else {
      out << " X" << std::endl;
    }
  }
}

int Monomer::write(char chain, int monomer_index,
		   char insertion_residue_code, int start_index, std::ostream &out) const {
  
 
  for (Atoms::const_iterator it= atoms_.begin(); it != atoms_.end(); ++it) {
    Atom_key al= it->key();
    //Point pt= res->cartesian_coords(al);
    const Atom &a= it->atom();
    Point pt = a.point();
    char alt=' ';
    //char chain=' ';
    out << boost::format(CGAL_PDB_INTERNAL_NS::atom_line_oformat_)
      % (start_index++) % Monomer::atom_key_string(al).c_str() % alt 
      % Monomer::type_string(type()).c_str() % chain % monomer_index % insertion_residue_code
      % pt.x() % pt.y() % pt.z() % a.occupancy() % a.temperature_factor() % a.segment_id().c_str()
      % a.element().c_str()% a.charge().c_str() << std::endl;
    //++anum;
  }
  return start_index;
}


 


Monomer::Monomer(Type al): atoms_(20), label_(al){
  Monomer_data::initialize();
  CGAL_assertion(al!= INV);
  //std::cout << "Inside size is " << sizeof(CGAL::PDB::Monomer) << std::endl;
}


Point Monomer::sidechain_point() const {
  double x=0; 
  double y=0; 
  double z=0;
  int count=std::distance(Monomer_data::amino_acid_data_[type()].extremes.begin(), 
			  Monomer_data::amino_acid_data_[type()].extremes.end());
  for (std::vector<Atom_key>::const_iterator it = Monomer_data::amino_acid_data_[type()].extremes.begin();
       it != Monomer_data::amino_acid_data_[type()].extremes.end(); ++it){
    x+= find(*it)->atom().point().x();
    y+= find(*it)->atom().point().y();
    z+= find(*it)->atom().point().z();
  }
  if (count == 0) {
    return find(AL_CA)->atom().point();
  } else {
    return Point(x/count, y/count, z/count);
  }
}


/*Monomer::Simplified_atoms_iterator Monomer::simplified_atoms_begin() const {
  return  Monomer_data::amino_acid_data_[type()].extremes.begin();
  }
  
  Monomer::Simplified_atoms_iterator Monomer::simplified_atoms_end() const{
  return  Monomer_data::amino_acid_data_[type()].extremes.end();
  }*/


Monomer::Atom_key Monomer::atom_key(const char *nm) {
  Monomer_data::initialize();
  char nn[5];
  sscanf(nm, "%4s", nn);
  std::string s(nn);
  for (unsigned int i=0; Monomer_data::clean_atom_name_data_[i].l != AL_INVALID; ++i){
    if (s==Monomer_data::clean_atom_name_data_[i].s){
      return Monomer_data::clean_atom_name_data_[i].l;
    }
  }
  CGAL_PDB_INTERNAL_NS::error_logger.new_warning(std::string(s + " is not a known atom type.").c_str());;
  return Monomer::AL_OTHER;
}





Atom::Type Monomer::element(Atom_key al){
  Monomer_data::initialize();
  for (unsigned int i=0; Monomer_data::atom_name_data_[i].l != Monomer::AL_INVALID; ++i){
    if (al==Monomer_data::atom_name_data_[i].l){
      return Monomer_data::atom_name_data_[i].t;
    }
  }
  CGAL_PDB_INTERNAL_NS::error_logger.new_internal_error("Unknown element label ");
  return Atom::INVALID;
}





std::string Monomer::atom_key_string(Atom_key al) {
  Monomer_data::initialize();
  for (unsigned int i=0; Monomer_data::atom_name_data_[i].l != Monomer::AL_INVALID; ++i){
    if (al==Monomer_data::atom_name_data_[i].l){
      return Monomer_data::atom_name_data_[i].s;
    }
  }

  std::ostringstream oss;
  oss << "Unknown atom label: " << al << " returning UNKN";
  CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());

  return "UNKN";
}






Monomer::Type Monomer::type(const std::string &s){
  /*switch(s[0]) {
    case 'A':
    switch(s[1]) {
    case 'C':
    return ACE;
    case 'L':
    return 
    default:
    return INV;
    };
    default:
    return INV;
    }*/
  if (s =="ACE") return ACE;
  if (s =="ALA") return ALA;
  if (s =="ARG") return ARG;
  if (s =="ASN") return ASN;
  if (s =="ASP") return ASP;
  if (s =="CYS") return CYS;
  if (s =="GLN") return GLN;
  if (s =="GLU") return GLU;
  if (s =="HIS") return HIS;
  if (s =="ILE") return ILE;
  if (s =="LEU") return LEU;
  if (s =="LYS") return LYS;
  if (s =="MET") return MET;
  if (s =="NH2") return NH2;
  if (s =="PHE") return PHE;
  if (s =="PRO") return PRO;
  if (s =="SER") return SER;
  if (s =="THR") return THR;
  if (s =="TRP") return TRP;
  if (s =="TYR") return TYR;
  if (s =="VAL") return VAL;
  if (s =="GLY") return GLY;
  
  if (s =="  A" || s =="ADE") return ADE;
  if (s =="  C" || s =="CYT") return CYT;
  if (s =="  G" || s =="GUA") return GUA;
  if (s =="  U" || s =="URA") return URA;
  if (s =="  T" || s =="THY") return THY;
  
  
  return INV;
}

std::string Monomer::type_string(Monomer::Type rl){
  switch (rl) {
  case GLY: return std::string("GLY");
  case ALA: return std::string("ALA");
  case VAL: return std::string("VAL");
  case LEU: return std::string("LEU");
  case ILE: return std::string("ILE");
  case SER: return std::string("SER");
  case THR: return std::string("THR");
  case CYS: return std::string("CYS");
  case MET: return std::string("MET");
  case PRO: return std::string("PRO");
  case ASP: return std::string("ASP");
  case ASN: return std::string("ASN");
  case GLU: return std::string("GLU");
  case GLN: return std::string("GLN");
  case LYS: return std::string("LYS");
  case ARG: return std::string("ARG");
  case HIS: return std::string("HIS");
  case PHE: return std::string("PHE");
  case TYR: return std::string("TYR");
  case TRP: return std::string("TRP");
  case ACE: return std::string("ACE");
  case NH2: return std::string("NH2");
  case ADE: return std::string("ADE");
  case CYT: return std::string("CYT");
  case GUA: return std::string("GUA");
  case URA: return std::string("URA");
  case THY: return std::string("THY");
  default: return std::string("UNK");
  }
    
}
CGAL_PDB_END_NAMESPACE
