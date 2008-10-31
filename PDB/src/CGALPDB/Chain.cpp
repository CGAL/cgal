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

#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/internal/Error_logger.h>
#include <CGAL/PDB/internal/pdb_utils.h>

#include <boost/format.hpp>

#include <sstream>
#include <cstdio>

CGAL_PDB_BEGIN_NAMESPACE
//Residue dummy_residue;
//Atom dummy_atom;

Chain::Chain(){}
 
void Chain::write_pdb(std::ostream &out) const {
  assert(!residues_.empty());
  for (unsigned int i=0; i< header_.size(); ++i){
    out << header_[i] << std::endl;
  }

  out << boost::format("MODEL %8d         ") % 1 << std::endl;
  
  //    int anum=1;
  write(' ' , 1, out);
  
  out << "ENDMDL                       " << std::endl;
}

int Chain::write(char chain, int start_index, std::ostream &out) const {
  // char line[81];
  //    int anum=1;
  Monomer_key last_resindex;
  Monomer::Type last_type= Monomer::INV;
  for (Monomer_const_iterator it = monomers_begin(); it != monomers_end(); ++it) {
    const Monomer &res= it->monomer();
    //Residue::Label rl =  res.label();
    //residues_[i]->atoms();
    start_index= res.write(chain, it->key().index(), ' ', start_index, out);
    
    /*IR_Map::const_iterator irit= insert_residues_.find(it->key());
    if (irit!= insert_residues_.end()) {
      for (unsigned int i=0; i< irit->data().size(); ++i){
        ir= irit->data().find(IR_key(i));
        CGAL_assertion(ir != irit->data().
	start_index= 
	  ir->data().write(chain, it->key().index(),
                           ir->key().index(), start_index, out);
      }
      }*/
    last_resindex= it->key();
    last_type= it->data().type();
  }
  const char *terformat="TER   %5d      %3s %c %3d%c";
  if (!residues_.empty()) {
    out << boost::format(terformat) % start_index
      % Monomer::type_string(last_type).c_str() %  chain 
      % last_resindex.index() % ' '<< std::endl;
  }
  return start_index+1;
}


std::vector<Monomer::Type> Chain::sequence() const{
  std::vector<Monomer::Type> ret; ret.reserve(residues_.size());
  for (Container::const_iterator it= residues_.begin(); it != residues_.end(); ++it){
    ret.push_back(it->monomer().type());
  }
  return ret;
}

unsigned int Chain::number_of_atoms() const {
  unsigned int ret=0;
  for (Container::const_iterator it= residues_.begin(); it != residues_.end(); ++it){
    ret += it->monomer().number_of_atoms();
  }
  return ret;
}
unsigned int Chain::number_of_bonds() const {
  unsigned int ret=0;
  for (Container::const_iterator it= residues_.begin(); it != residues_.end(); ++it){
    ret += it->monomer().number_of_bonds();
  }
  return ret;
}


void Chain::swap_with(Chain &o) {
  std::swap(residues_, o.residues_);
  std::swap(header_, o.header_);
  std::swap(insert_residues_, o.insert_residues_);
}


void Chain::dump(std::ostream &out) const {
 for (Container::const_iterator it= residues_.begin(); it != residues_.end(); ++it){
   out << "Residue " << it->key() << std::endl;
    it->monomer().dump(out);
  }
}

void Chain::set_has_bonds(bool tf) {
  if (tf) {
    for (Container::iterator it= residues_.begin(); it != residues_.end(); ++it){
      it->monomer().set_has_bonds(tf);
    }
  }
}

Chain::Monomer_iterator Chain::insert_internal(Monomer_key k, const Monomer &m) {
  if (residues_.find(k) != residues_.end()){
    std::ostringstream eout;
    eout << "Warning, newly added residue has index "<< k
	 << " which already exists.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(eout.str().c_str());
  }
  return residues_.insert(Container::value_type(k,m));
}

bool Chain::has_bonds() const {
  for (Monomer_const_iterator it= monomers_begin();
       it != monomers_end(); ++it){
    if (!it->monomer().has_bonds()) return false;
  }
  return true;
}

#if 0
Chain::Graph Chain::graph() const {
  std::vector<Point> points;
  std::vector<std::pair<int,int> > edges;
  //std::map<int,int> index_index_map;
  int max_index=-1;
  for (unsigned int i=0; i< residues_.size(); ++i){
    const Residue &r= residues_[i];
    //dsr::vector<Residue::Atom_label> als= r->atoms();
    //dsr::vector<Residue::Bond> bls= r->bonds();

    /*for (Residue::Atoms_iterator it= r->atoms_begin(); it != r->atoms_end(); ++it){
      index_index_map[it->index()]= index_index_map.size();
      }*/
    for (Residue::Bonds_iterator bit= r.bonds_begin(); bit != r.bonds_end(); ++bit){
      edges.push_back(std::pair<int,int>(bit->key(),
					 bit->second));
      max_index= (std::max)(max_index, bit->key());
      max_index= (std::max)(max_index, bit->second);
    }
    if (i!= 0){
      edges.push_back(std::pair<int,int>(r.atom(Residue::AL_N).index(),
					 residues_[i-1].atom(Residue::AL_C).index()));
    }
  }
      
  points.resize(max_index);
  for (unsigned int i=0; i< residues_.size(); ++i){
    const Residue &r= residues_[i];
    for (Residue::Atoms_iterator it= r.atoms_begin(); it != r.atoms_end(); ++it){
      points[it->monomer().index()]= it->atom().cartesian_coords();
    }
  }
  return Graph(points,edges);
}
  
std::vector<Point> Chain::backbone() const {
  // skip ACE, probably should do something more clever
  std::vector<Point> pts;
  for (Const_residues_iterator it= residues_begin(); it != residues_.end(); ++it) {
    const Residue &aa= *it;
		
    for (Residue::Atoms_iterator it= aa.atoms_begin(); it != aa.atoms_end(); ++it){
      if (it->first == Residue::AL_C || it->first== Residue::AL_CA || it->first== Residue::AL_N) {
	Point pt=it->atom().cartesian_coords();
	pts.push_back(pt);
      }
    }
  }
  return pts;
}
#endif


/*const Residue& Chain::residue(Residue::Index i) const{
  unsigned int cur= residue_offset(i);
  if (cur == residues_.size()){
    std::ostringstream oss;
    oss << "residue(int) called with index that does not correspond to a valid residue: " << i;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return dummy_residue_;
  } else {
    return residues_[cur];
  }
  }*/


/*unsigned int Chain::residue_offset(Residue::Index i) const {
  unsigned int cur= residues_.size();
  if (!residues_.empty()){
    cur = std::min BOOST_PREVENT_MACRO_SUBSTITUTION(i.index(), cur-1);
    if (residues_[cur].index() > i) {
      do {
	--cur;
      } while (cur > 0 && residues_[cur].index() > i);
    } else if (residues_[cur].index() < i) {
      do {
	--cur;
      } while (cur < residues_.size() && residues_[cur].index() < i);
    }
    if (residues_[cur].index() != i) cur= residues_.size();
  }
  return cur;
  }*/

/*unsigned int Chain::residue_offset_of_atom_index(Atom::Index index) const {
  for (int i= residues_.size()-1; i >=0; --i) {
    if (residues_[i].min_atom_index() <= index) return i;
  }
  return residues_.size();
  }*/

/*void Chain::set_atom(Atom::Index index, const Atom &a) {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
    std::ostringstream oss;
    oss << "set_atom called with index " << index << " which does not corresponding to an existing atom.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
  } else {
    residues_[ind].atoms_iterator_from_index(index)->second=a;
  }
  //if (!min_atom_index_ || a.index() < min_atom_index_) min_atom_index_=a.index();
  }*/

/*const Atom& Chain::atom(Atom::Index index) const {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
    std::ostringstream oss;
    oss << "set_atom called with index " << index << " which does not corresponding to an existing atom.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return internal::dummy_atom;
  } else {
    return residues_[ind].atoms_iterator_from_index(index)->second;
  }
}

void Chain::erase_atom(Atom::Index index) {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
  } else {
    residues_[ind].erase_atom(residues_[ind].atom_label(index));
  }
  }*/

/*Atom::Index Chain::parent_atom(Atom::Index index) const {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
    std::ostringstream oss;
    oss << "parent_atom called with index " << index << " which does not corresponding to an existing atom.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return Atom::Index();
  } else {
    Residue::Atom_label al= residues_[ind].atoms_iterator_from_index(index)->first;
    if (al != Residue::AL_CA && al != Residue::AL_C && al != Residue::AL_N) {
      std::ostringstream oss;
      oss << "parent_atom called with an atom which was not an N, CA or C.";
      CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
      return Atom::Index();
    }
    if (al == Residue::AL_C) {
      return residues_[ind].atom(Residue::AL_CA).index();
    } else if (al == Residue::AL_CA) {
      return residues_[ind].atom(Residue::AL_N).index();
    } else {
      if (ind==0) return Atom::Index();
      else return residues_[ind-1].atom(Residue::AL_C).index();
    }
  }
}

Spherical_point Chain::spherical_coordinates(Atom::Index index) const {
  Atom::Index p= parent_atom(index);
  Point pc(100000, 100000, 100000);
  Point popc(100000, -100000, 100000);
  Point popopc(100000, -100000, -100000);
  if (p != Atom::Index()) {
    pc= atom(p).cartesian_coords();
    Atom::Index pop= parent_atom(p);
    if (pop != Atom::Index()) {
      popc= atom(pop).cartesian_coords();
      Atom::Index popop= parent_atom(pop);
      if (popop != Atom::Index()){
	popopc= atom(popop).cartesian_coords();
      }
    }
  }
  Construct_spherical_point csp(pc, popc, popopc);
  return csp(atom(index).cartesian_coords());
  }*/
  

/*const Residue& Chain::residue_containing_atom(Atom::Index atom_index) const{
  unsigned int rindex= residue_offset_of_atom_index(atom_index);
  if (rindex == residues_.size() || residues_[rindex].atom_label(atom_index) == Residue::AL_INVALID) {
    std::ostringstream oss;
    oss << "Chain::atom_label_of_atom() called with uninitialized atom " << atom_index;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return internal::dummy_residue;
  }
  return residues_[rindex];
}
  
Residue& Chain::residue_containing_atom(Atom::Index atom_index) {
  unsigned int rindex= residue_offset_of_atom_index(atom_index);
  if (rindex == residues_.size()  || residues_[rindex].atom_label(atom_index) == Residue::AL_INVALID) {
    std::ostringstream oss;
    oss << "Chain::atom_label_of_atom() called with uninitialized atom " << atom_index;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return internal::dummy_residue;
  }
  return residues_[rindex];
  }*/


CGAL_PDB_END_NAMESPACE
