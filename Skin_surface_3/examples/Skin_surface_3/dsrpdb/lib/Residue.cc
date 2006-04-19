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

You should have received a copy of the GNU Lesser General Public
License along with the DSR PDB Library; see the file COPYING.LIB.  If
not, write to the Free Software Foundation, Inc., 59 Temple Place -
Suite 330, Boston, MA 02111-1307, USA. */

#include "dsrpdb/Residue.h"
#include "Residue_data.h"
#include <algorithm>
#include <iostream>
#include <cassert>
#include <dsrpdb_internal/Error_logger.h>
#include <limits>
#include <sstream>

namespace dsrpdb {
  

  static Atom dummy_atom_;

  void Residue::set_has_bonds(bool tf) {
    typedef Residue_data::Possible_bond Possible_bond;

    if (!tf) {
      bonds_.clear();
    } else {
      const std::vector<Possible_bond >& bonds= Residue_data::amino_acid_data_[type()].bonds;
      for (unsigned int i=0; i< bonds.size(); ++i){
	Const_atoms_iterator itf= atoms_.find(bonds[i].first);
	Const_atoms_iterator its= atoms_.find(bonds[i].second);
	if ( itf != atoms_.end() && its != atoms_.end()){
	  bonds_.push_back(Bond(itf->second.index(), its->second.index()));
	} 
      }
    }
  }






  bool Residue::has_atom(Residue::Atom_label ial) const {
    Residue::Atom_label al= Residue_data::fix_atom_label(label_, ial);
    assert(can_have_atom(al));
    return atoms_.find(al) != atoms_.end();
  }









  bool Residue::can_have_atom(Residue::Atom_label ial) const {
    if (ial== AL_INVALID) return false; // hack to avoid some circularity
    Residue::Atom_label al= Residue_data::fix_atom_label(label_, ial);
    //assert(Residue_data::amino_acid_data_.find(label_)!= Residue_data::amino_acid_data_.end());
    for (unsigned int i=0; i< Residue_data::amino_acid_data_[label_].atoms.size(); ++i){
      if (Residue_data::amino_acid_data_[label_].atoms[i] == al) return true;
    }
    return false;
  }




 



  const Atom &Residue::atom(Residue::Atom_label al) const {
    Residue::Atom_label fal= Residue_data::fix_atom_label(label_, al);
    //int fa= find_atom(al);
    if (atoms_.find(fal) != atoms_.end()) return atoms_.find(fal)->second;
    else {
      return dummy_atom_;
    }
  }
    





  void Residue::set_atom(Residue::Atom_label ial, const Atom &a) {
    Residue::Atom_label al= Residue_data::fix_atom_label(label_, ial);
    if (!can_have_atom(al)){
      dsrpdb_internal::error_logger.new_warning((std::string("Trying to set invalid atom ") + atom_label_string(ial) 
						 + " on a residue of type " + type_string(label_)).c_str());
    }
    if (al == AL_INVALID) {
      return;
    }
    //assert(atoms_.find(al)== atoms_.end());
    //bonds_.clear();
    atoms_[al]=a;
    // This is a bit evil, I haven't figured out a better way though.
    atoms_[al].set_type(element(al));
    if (min_atom_index_) {
      min_atom_index_= std::min(min_atom_index_, a.index());
    } else {
      min_atom_index_= a.index();
    }
    if (has_bonds()){
      set_has_bonds(false);
      set_has_bonds(true);
    }
    //atoms_.push_back(a);
  }








  Residue::Type Residue::type() const {
    return label_;
  }

  void Residue::dump(std::ostream &out) const {
    out << "Type: " << type_string(type()) << std::endl;
    const std::vector<Residue::Atom_label> &valid_atoms= Residue_data::amino_acid_data_[type()].atoms;
    for (unsigned int i=0; i< valid_atoms.size(); ++i){
      Residue::Atom_label al= valid_atoms[i];
      out << Residue::atom_label_string(al); //Residue::write_atom_label(al, out);
      if (has_atom(al)){
	out << " (" << atom(al).cartesian_coords() << ") " << atom(al).index() << std::endl;
      } else {
	out << " X" << std::endl;
      }
    }
  }

  struct Compare_index{
    bool operator()(const std::pair<Residue::Atom_label, Atom> &a,
		    const std::pair<Residue::Atom_label, Atom>&b) const {
      return a.second.index() < b.second.index();
    }
  };

  void Residue::write(char chain, std::ostream &out) const {
    static const char atom_line_oformat_[]=
      "ATOM  %5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s";
    char line[81];
    std::vector<std::pair<Atom_label, Atom> > atoms(atoms_.begin(), atoms_.end());
    std::sort(atoms.begin(), atoms.end(), Compare_index());

    for (unsigned int i=0; i< atoms.size(); ++i) {
	Atom_label al= atoms[i].first;
	//Point pt= res->cartesian_coords(al);
	const Atom &a= atoms[i].second;
	Point pt = a.cartesian_coords();
	char alt=' ';
	//char chain=' ';
	char insertion_residue_code=' ';
	sprintf(line, atom_line_oformat_,
		static_cast<unsigned int>(a.index()), Residue::atom_label_string(al).c_str(), alt,
		Residue::type_string(type()).c_str(), chain, static_cast<unsigned int>(index()), insertion_residue_code,
		pt.x(), pt.y(), pt.z(), a.occupancy(), a.temperature_factor(), a.segment_id(),
		a.element(), a.charge());
	out << line << std::endl;
	//++anum;
      }
  }
  
  Atom::Index Residue::last_atom_index() const {
    Atom::Index max= atoms_.begin()->second.index();
    for ( Const_atoms_iterator it= atoms_.begin(); it != atoms_.end();
	  ++it) {
      max= std::max(max, it->second.index());
    }
    
    return max;
  }

  unsigned int Residue::number_of_atoms() const {
    return atoms_.size();
  }

  unsigned int Residue::number_of_bonds() const {
    bonds_begin();
    return bonds_.size();
  }

  void Residue::set_index(Residue::Index i) {
    index_=i;
  }

 
  Residue::Atom_label Residue::atom_label(Atom::Index model_index) const{
    Const_atoms_iterator it= atoms_iterator_from_index(model_index);
    if (it != atoms_end()) {
      return it->first;
    } else {
      return AL_INVALID;
    }
  }


 
  // protected
  Atom::Index Residue::index(Residue::Atom_label al) const {
    return atom(al).index();
  }
 
  Residue::Atoms_iterator Residue::atoms_iterator_from_index(Atom::Index ind) {
    for (Atoms_iterator it= atoms_begin(); it != atoms_end(); ++it){
      if (it->second.index()==ind) return it;
    }
    dsrpdb_internal::error_logger.new_warning("Invalid atom index used to request atom from residue.");

    return atoms_end();
  }
  Residue::Const_atoms_iterator Residue::atoms_iterator_from_index(Atom::Index ind) const {
    for (Const_atoms_iterator it= atoms_begin(); it != atoms_end(); ++it){
      if (it->second.index()==ind) return it;
    }
    dsrpdb_internal::error_logger.new_warning("Invalid atom index used to request atom from residue.");

    return atoms_end();
  }

 


  Residue::Residue(Type al): atoms_(20), label_(al){
    Residue_data::initialize();
    assert(al!= INV);
  }


  Point Residue::sidechain_point() const {
    double x=0; 
    double y=0; 
    double z=0;
    int count=std::distance(Residue_data::amino_acid_data_[type()].extremes.begin(), 
			    Residue_data::amino_acid_data_[type()].extremes.end());
    for (std::vector<Atom_label>::const_iterator it = Residue_data::amino_acid_data_[type()].extremes.begin();
	 it != Residue_data::amino_acid_data_[type()].extremes.end(); ++it){
      x+= atom(*it).cartesian_coords().x();
      y+= atom(*it).cartesian_coords().y();
      z+= atom(*it).cartesian_coords().z();
    }
    if (count == 0) {
      return atom(AL_CA).cartesian_coords();
    } else {
      return Point(x/count, y/count, z/count);
    }
  }


  /*Residue::Simplified_atoms_iterator Residue::simplified_atoms_begin() const {
    return  Residue_data::amino_acid_data_[type()].extremes.begin();
  }
  
  Residue::Simplified_atoms_iterator Residue::simplified_atoms_end() const{
    return  Residue_data::amino_acid_data_[type()].extremes.end();
    }*/


  Residue::Atom_label Residue::atom_label(const char *nm) {
    Residue_data::initialize();
    char nn[5];
    sscanf(nm, "%4s", nn);
    std::string s(nn);
    for (unsigned int i=0; Residue_data::clean_atom_name_data_[i].l != AL_INVALID; ++i){
      if (s==Residue_data::clean_atom_name_data_[i].s){
	return Residue_data::clean_atom_name_data_[i].l;
      }
    }
    dsrpdb_internal::error_logger.new_warning(std::string(s + " is not a known atom type.").c_str());;
    return Residue::AL_OTHER;
  }





  Atom::Type Residue::element(Atom_label al){
    Residue_data::initialize();
    for (unsigned int i=0; Residue_data::atom_name_data_[i].l != Residue::AL_INVALID; ++i){
      if (al==Residue_data::atom_name_data_[i].l){
	return Residue_data::atom_name_data_[i].t;
      }
    }
    dsrpdb_internal::error_logger.new_internal_error("Unknown element label ");
    return Atom::INVALID;
  }





  std::string Residue::atom_label_string(Atom_label al) {
    Residue_data::initialize();
    for (unsigned int i=0; Residue_data::atom_name_data_[i].l != Residue::AL_INVALID; ++i){
      if (al==Residue_data::atom_name_data_[i].l){
	return Residue_data::atom_name_data_[i].s;
      }
    }

    std::ostringstream oss;
    oss << "Unknown atom label: " << al << " returning UNKN";
    dsrpdb_internal::error_logger.new_warning(oss.str().c_str());

    return "UNKN";
  }






  Residue::Type Residue::type(const std::string &s){
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
    else return INV;
  }

  std::string Residue::type_string(Residue::Type rl){
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
    default: return std::string("UNK");
    }
    
  }


  
};
