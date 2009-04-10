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
#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/range.h>
#include <fstream>
#include <cassert>

template <class Vector, class Range>
void fill(Vector &v, const Range &r) {
v.insert(v.end(), r.begin(), r.end());
}

int main() {

 //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  using namespace CGAL::PDB;
  typedef Atom::Index Index;
  typedef std::pair<Index,Index> Bond;
  {
    std::ifstream in("data/check_pdb.pdb");
    PDB pdb(in);
    assert(pdb.models().size() != 0);
    Model &model= pdb.models().begin()->model();
    assert(model.chains().size() != 0);
    Chain &p= model.chains().begin()->chain();
    unsigned int na= p.atoms().size();
    //p.write(std::cout);
    std::cout << "There are " << p.monomers().size() << " residues." 
	      << std::endl;
  
    std::vector<Atom> atoms;
	fill(atoms, make_atom_range(p.atoms()));
    assert(na==atoms.size());
    std::vector<Point> points;
	fill(points,make_point_range(make_atom_range(p.atoms())));
 
    index_atoms(make_backbone_range(p.atoms()));
    std::vector<Index> bbi;
	fill(bbi, make_index_range(make_atom_range(make_backbone_range(p.atoms()))));
    assert(bbi.size() < atoms.size());
    std::cout << bbi.size() << std::endl;
    std::vector<Bond> bbs;
	fill(bbs, make_bond_indices_range(make_filtered_bond_range(Is_backbone(), p.bonds())));
    std::cout << bbs.size() << std::endl;
    for (unsigned int i=0; i< bbs.size(); ++i) {
      assert(bbs[i].first < Index(bbi.size()));
      assert(bbs[i].second < Index(bbi.size()));
    }
  }
  
  return 0;
}
