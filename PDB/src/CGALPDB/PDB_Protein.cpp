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

#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/Residue.h>
#include <CGAL/PDB/internal/Error_logger.h>
#include <sstream>
CGAL_PDB_BEGIN_NAMESPACE

static Residue dummy_residue_;
static Atom dummy_atom_;

Protein::Protein(): chain_(' '){}

char Protein::chain() const {return chain_;}
  
void Protein::set_chain(char c) {chain_=c;}

/*bProtein::Protein(const dsr::vector<Residue_label> &seq) {
  for (unsigned int i=0; i< seq.size(); ++i){
  residues_.push_back(Residue::new_residue(seq[i]));
  }
  model_=1;
  };*/


std::vector<Residue::Type> Protein::sequence() const{
  std::vector<Residue::Type> ret(residues_.size());
  for (unsigned int i=0; i< residues_.size(); ++i){
    ret[i]= residues_[i].type();
  }
  return ret;
}

unsigned int Protein::number_of_atoms() const {
  unsigned int ret=0;
  for (unsigned int i=0; i< residues_.size(); ++i){
    ret += residues_[i].number_of_atoms();
  }
  return ret;
}
unsigned int Protein::number_of_bonds() const {
  unsigned int ret=0;
  for (unsigned int i=0; i< residues_.size(); ++i){
    ret += residues_[i].number_of_bonds();
  }
  return ret;
}



void Protein::dump(std::ostream &out) const {
  for (unsigned int i=0; i< residues_.size(); ++i){
    out << "Residue " << residues_[i].index() << std::endl;
    residues_[i].dump(out);
  }
}

void Protein::new_residue(const Residue &res){
  if (!residues_.empty() && res.index() <= residues_.back().index()){
    std::ostringstream eout;
    eout << "Warning, newly added residue has index "<< res.index() 
	 << " while previous residue has index " << residues_.back().index();
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(eout.str().c_str());
  }
  if (!residues_.empty() && !residues_.back().has_atom(Residue::AL_C)) {
    std::ostringstream eout;
    eout << "Warning, newly added residue " << residues_.back().index() 
	 << " either missing atom C or atoms out of order in pdb.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(eout.str().c_str());
  }
  if (!residues_.empty() && !residues_.back().has_atom(Residue::AL_N)) {
    std::ostringstream eout;
    eout << "Warning, newly added residue " << residues_.back().index() 
	 << " either missing atom N or atoms out of order in pdb.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(eout.str().c_str());
  }
  if (!residues_.empty() && !residues_.back().has_atom(Residue::AL_CA)) {
    std::ostringstream eout;
    eout << "Warning, newly added residue " << residues_.back().index() 
	 << " either missing atom CA or atoms out of order in pdb. ";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(eout.str().c_str());
  }
  residues_.push_back(res);
  //residues_.back().write('t', std::cout);
}

#if 0
Protein::Graph Protein::graph() const {
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
      edges.push_back(std::pair<int,int>(bit->first,
					 bit->second));
      max_index= (std::max)(max_index, bit->first);
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
      points[it->second.index()]= it->second.cartesian_coords();
    }
  }
  return Graph(points,edges);
}
  
std::vector<Point> Protein::backbone() const {
  // skip ACE, probably should do something more clever
  std::vector<Point> pts;
  for (Const_residues_iterator it= residues_begin(); it != residues_.end(); ++it) {
    const Residue &aa= *it;
		
    for (Residue::Atoms_iterator it= aa.atoms_begin(); it != aa.atoms_end(); ++it){
      if (it->first == Residue::AL_C || it->first== Residue::AL_CA || it->first== Residue::AL_N) {
	Point pt=it->second.cartesian_coords();
	pts.push_back(pt);
      }
    }
  }
  return pts;
}
#endif

Protein::Atoms_iterator Protein::atoms_begin() {
  return Atoms_iterator(residues_.begin(), residues_.end());
}
Protein::Atoms_iterator Protein::atoms_end() {
  return Atoms_iterator(residues_.end(), residues_.end());
}

Protein::Const_atoms_iterator Protein::atoms_begin() const{
  return Const_atoms_iterator(residues_.begin(), residues_.end());
}
Protein::Const_atoms_iterator Protein::atoms_end() const{
  return Const_atoms_iterator(residues_.end(), residues_.end());
}
Protein::Bonds_iterator Protein::bonds_begin() const{
  CGAL_PDB_INTERNAL_NS::error_logger.new_warning("bonds_begin() called without has_bonds() being true.\n");
  return Bonds_iterator(residues_.begin(), residues_.end());
}
Protein::Bonds_iterator Protein::bonds_end() const{
  return Bonds_iterator(residues_.end(), residues_.end());
}

const Residue& Protein::residue(Residue::Index i) const{
  unsigned int cur= residue_offset(i);
  if (cur == residues_.size()){
    std::ostringstream oss;
    oss << "residue(int) called with index that does not correspond to a valid residue: " << i;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return dummy_residue_;
  } else {
    return residues_[cur];
  }
}

bool Protein::has_residue(Residue::Index i) const {
  return residue_offset(i) != residues_.size();
}

unsigned int Protein::residue_offset(Residue::Index i) const {
  unsigned int cur= residues_.size();
  if (!residues_.empty()){
    cur = std::min BOOST_PREVENT_MACRO_SUBSTITUTION(i.to_index(), cur-1);
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
}

unsigned int Protein::residue_offset_of_atom_index(Atom::Index index) const {
  for (int i= residues_.size()-1; i >=0; --i) {
    if (residues_[i].min_atom_index() <= index) return i;
  }
  return residues_.size();
}

void Protein::set_atom(Atom::Index index, const Atom &a) {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
    std::ostringstream oss;
    oss << "set_atom called with index " << index << " which does not corresponding to an existing atom.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
  } else {
    residues_[ind].atoms_iterator_from_index(index)->second=a;
  }
  //if (!min_atom_index_ || a.index() < min_atom_index_) min_atom_index_=a.index();
}

const Atom& Protein::atom(Atom::Index index) const {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
    std::ostringstream oss;
    oss << "set_atom called with index " << index << " which does not corresponding to an existing atom.";
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return dummy_atom_;
  } else {
    return residues_[ind].atoms_iterator_from_index(index)->second;
  }
}

void Protein::erase_atom(Atom::Index index) {
  unsigned int ind= residue_offset_of_atom_index(index);
  if (ind == residues_.size()) {
  } else {
    residues_[ind].erase_atom(residues_[ind].atom_label(index));
  }
}

Atom::Index Protein::parent_atom(Atom::Index index) const {
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

Spherical_point Protein::spherical_coordinates(Atom::Index index) const {
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
}
  

const Residue& Protein::residue_containing_atom(Atom::Index atom_index) const{
  unsigned int rindex= residue_offset_of_atom_index(atom_index);
  if (rindex == residues_.size() || residues_[rindex].atom_label(atom_index) == Residue::AL_INVALID) {
    std::ostringstream oss;
    oss << "Protein::atom_label_of_atom() called with uninitialized atom " << atom_index;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return dummy_residue_;
  }
  return residues_[rindex];
}
  
Residue& Protein::residue_containing_atom(Atom::Index atom_index) {
  unsigned int rindex= residue_offset_of_atom_index(atom_index);
  if (rindex == residues_.size()  || residues_[rindex].atom_label(atom_index) == Residue::AL_INVALID) {
    std::ostringstream oss;
    oss << "Protein::atom_label_of_atom() called with uninitialized atom " << atom_index;
    CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
    return dummy_residue_;
  }
  return residues_[rindex];
}

/*Residue::Atom_label Protein::atom_label_of_atom(int atom_index) const {
  int rindex= residue_index_of_atom_index(atom_index);
  Residue::Const_atoms_iterator it= residues_[rindex].atoms_iterator_from_index(atom_index);
  if (it == residues_[rindex].atoms_end()) {
  std::ostringstream oss;
  oss << "Protein::atom_label_of_atom() called with uninitialized atom " << atom_index;
  CGAL_PDB_INTERNAL_NS::error_logger.new_warning(oss.str().c_str());
  return Residue::AL_INVALID;
  } else {
  return it->first;
  }
  }*/

#if 0
Protein::Backbone_coordinates_iterator Protein::backbone_coordinates_begin() const{
  if (residues_.empty()){
    Backbone_coordinates_iterator b(residues_.begin(), residues_.end());
    Backbone_coordinates_iterator e(residues_.begin(), residues_.end());
    assert(b==e);
  }
  return Backbone_coordinates_iterator(residues_.begin(), residues_.end());
}
Protein::Backbone_coordinates_iterator Protein::backbone_coordinates_end() const{
  return Backbone_coordinates_iterator(residues_.end(), residues_.end());
}
#endif
CGAL_PDB_END_NAMESPACE
