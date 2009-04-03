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

#include <CGAL/PDB/Heterogen.h>
#include <CGAL/PDB/internal/Error_logger.h>
#include <CGAL/PDB/internal/pdb_utils.h>
#include <CGAL/basic.h>

#include <boost/format.hpp>

#include <iostream>
#include <limits>

#include <sstream>
#include <cstdio>
namespace CGAL { namespace PDB {


void Heterogen::copy_from(const Heterogen &o) {
  atoms_= o.atoms_;
  type_= o.type_;
  bonds_.clear();
  for (unsigned int i=0; i< o.bonds_.size(); ++i) {
    std::string nma=o.bonds_[i].first.key();
    std::string nmb=o.bonds_[i].second.key();
    Atom_consts::iterator ita=find(nma);
    Atom_consts::iterator itb=find(nmb);
    CGAL_assertion(ita != atoms().end());
    CGAL_assertion(itb != atoms().end());
    // needed due to idiotic C++ syntax
    Bond_endpoint bea(ita);
    Bond_endpoint beb(itb);
    Bond nb(bea, beb);
    bonds_.push_back(nb);
  }
  chain_=o.chain_;
}


bool Heterogen::connect(Atom::Index a, Atom::Index b) {
  Atoms::iterator ita=atoms().end(), itb= atoms().end();
  if (b < a) std::swap(a,b);
  for (Atoms::iterator it= atoms().begin(); it != atoms().end(); ++it) {
    if (it->atom().index() == a) {
      ita= it;
    }
    if (it->atom().index() == b) {
      itb= it;
    }
  }
  if (ita == atoms().end()) return false;
  if (itb == atoms().end()) return false;
  for (unsigned int i=0; i< bonds_.size(); ++i) {
    if (bonds_[i].first.key() == ita->key() 
        && bonds_[i].second.key() == itb->key()){
      return true;
    }
  }
  bonds_.push_back(Bond(Bond_endpoint(ita),
                        Bond_endpoint(itb)));
  return true;
}

void Heterogen::swap_with(Heterogen &o) {
  std::swap(atoms_, o.atoms_);
  std::swap(bonds_, o.bonds_);
  std::swap(type_, o.type_);
  std::swap(chain_, o.chain_);
}




void Heterogen::dump(std::ostream & /* out */ ) const {
  
}

int Heterogen::write(std::string name, int num, 
                     int start_index, std::ostream &out) const {
  CGAL_PDB_FOREACH(const Atom_pair &aa, atoms()) {
    Atom_key al= aa.key();
    //Point pt= res->cartesian_coords(al);
    const Atom &a= aa.atom();
    a.set_index(Atom::Index(start_index));
    Point pt = a.point();
    //char chain=' ';

    //"HETATM%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s";
    out << boost::format(internal::hetatom_line_oformat_)
      % (start_index++) % al.c_str()
      % name.c_str() % num
      % pt.x() % pt.y() % pt.z() % a.occupancy() % a.temperature_factor()
      % a.element().c_str() % "  " << std::endl;
    //++anum;
  }
  std::map<int, std::vector<int> > connects;
  for (unsigned int i=0; i< bonds_.size(); ++i) {
    int ia= bonds_[i].first.atom().index().index();
    int ib= bonds_[i].second.atom().index().index();
    connects[ia].push_back(ib);
    connects[ib].push_back(ia);
  }
  for (std::map<int, std::vector<int> >::const_iterator it= connects.begin();
       it != connects.end(); ++it) {
    out << "CONECT " << it->first;
    for (unsigned int i=0; i< it->second.size(); ++i) {
      out << " " << it->second[i];
    }
    out << std::endl;
  }
  return start_index;
}


 


Heterogen::Heterogen(std::string name): atoms_(20), type_(name){
}




}}
