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
#include <CGAL/PDB/iterator.h>
#include <fstream>

int main() {

 //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  using namespace CGAL::PDB;
  typedef Atom::Index Index;
  typedef std::pair<Index,Index> Bond;
  {
    std::ifstream in("data/check_protein.pdb");
    CGAL_PDB_NS::Chain p(in);
    unsigned int na= p.number_of_atoms();
    //p.write(std::cout);
    std::cout << "There are " << p.number_of_monomers() << " residues." 
	      << std::endl;
  
    std::vector<Atom> atoms (make_atom_iterator(p.atoms_begin()),
			    make_atom_iterator(p.atoms_end()));
    CGAL_assertion(na==atoms.size());
    std::vector<Point> points (make_point_iterator(make_atom_iterator(p.atoms_begin())),
			      make_point_iterator(make_atom_iterator(p.atoms_end())));
 
    index_atoms(make_backbone_iterator(p.atoms_begin(), p.atoms_end()),
		make_backbone_iterator(p.atoms_end(), p.atoms_end()));
    std::vector<Index> bbi (make_index_iterator(make_atom_iterator(make_backbone_iterator(p.atoms_begin(), p.atoms_end()))),
			   make_index_iterator(make_atom_iterator(make_backbone_iterator(p.atoms_end(), p.atoms_end()))));
    CGAL_assertion(bbi.size() < atoms.size());
    std::cout << bbi.size() << std::endl;
    std::vector<Bond> bbs (make_bond_indices_iterator(make_ok_bond_iterator(Is_backbone(), p.bonds_begin(), p.bonds_end())),
			  make_bond_indices_iterator(make_ok_bond_iterator(Is_backbone(), p.bonds_end(), p.bonds_end())));
    std::cout << bbs.size() << std::endl;
    for (unsigned int i=0; i< bbs.size(); ++i) {
      CGAL_assertion(bbs[i].first < Index(bbi.size()));
      CGAL_assertion(bbs[i].second < Index(bbi.size()));
    }
  }
  
  return 0;
}
