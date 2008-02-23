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

#ifndef CGAL_DSR_PDB_HETEROGEN_H
#define CGAL_DSR_PDB_HETEROGEN_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/PDB/Atom.h>
#include <CGAL/PDB/small_map.h>
#include <CGAL/Tools/Label.h>
#include <vector>
#include <map>

CGAL_PDB_BEGIN_NAMESPACE


//! This class defines the value_type of the Atom_iterator
class Heterogen_Atom_iterator_value_type:
  public small_map_value_type<std::string,
                              Atom> {
  typedef small_map_value_type<std::string, Atom> P;
public:
  Heterogen_Atom_iterator_value_type(std::string k,
                                     const Atom &a=Atom()): P(k,a){}
  Heterogen_Atom_iterator_value_type(){}
  const Atom &atom() const {return P::data();}
  Atom &atom() {return P::data();}
};

/*bool operator<(const std::string &s, 
               const Heterogen_Atom_iterator_value_type &vt) {
  return s < vt.key();
  }*/

//! The class representing a residue or nucleotide.
/*!  All the information concerning atoms and bonds for each monomer
  is stored here. To add atoms to monomers, new monomers, or bonds
  to monomers, look in Monomer_data.cpp. There is documentation there
  of what you need to do.
*/
class Heterogen {
  struct Heterogen_type_tag{};
  struct Atom_type_tag{};
public:
  typedef std::string Atom_key;

 
  typedef small_map<Heterogen_Atom_iterator_value_type> Atoms;


  //! Default constructor. Makes and invalid monomer.
  Heterogen(){}

  //! Make a monomer of a given type
  Heterogen(std::string name);

  CGAL_COPY_CONSTRUCTOR(Heterogen);

  void copy_from(const Heterogen &o);

  void swap_with(Heterogen &o);

  //! The label for the monomer
  std::string type() const{
    return type_;
  }

  CGAL_GETNR(char, chain, return chain_;)
  CGAL_SET(char, chain, chain_=k;)
  

  CGAL_CONST_ITERATOR(Atom, atom, Atoms::const_iterator, 
			  return atoms_.begin(),
			  return atoms_.end());
  CGAL_ITERATOR(Atom, atom, Atoms::iterator, 
			  return atoms_.begin(),
			  return atoms_.end());
  
  CGAL_SIZE(atoms, return atoms_.size());

  CGAL_FIND(Atom, {
      return atoms_.find(k);
    });

  CGAL_INSERT(Atom, {bonds_.clear(); 
      return atoms_.insert(Atoms::value_type(k,m));});

   //! Write it for debugging
  void dump(std::ostream &out) const;

  //! Write it for debugging
  std::ostream& write(std::ostream &out) const {dump(out);return out;}


  //! Write the lines for a pdb file
  /*!
    Indices start at the start_index and the new start_index is returned.
  */
  int write(std::string name, int num, int start_index, std::ostream &out) const;    


  //! The Bond_iterator value_type is a pair of these
  class Bond_endpoint {
    Atom_const_iterator aci_;
    Bond_endpoint(Atom_const_iterator aci): aci_(aci){}
    friend class Heterogen;
  public:
    Bond_endpoint(){}
    Atom_key key() const {return aci_->key();}
    const Atom &atom() const {return aci_->atom();}
  };



  //! A bond between two atoms in a monomer.
  /*!
    The ints refer the the atom index.
  */
  typedef std::pair<Bond_endpoint, Bond_endpoint> Bond;

  //! Return a list of all the bonds in the monomer
  CGAL_CONST_ITERATOR(Bond, bond, 
                      std::vector<Bond>::const_iterator,
                      return bonds_.begin(),
                      return bonds_.end());
  

  //! The number of atoms present in the monomer
  CGAL_SIZE(bonds, return bonds_.size());

  bool connect(Atom::Index a, Atom::Index b);

private:

  /*--------------- Data members-----------------------------------*/
  Atoms atoms_;
  std::vector<Bond> bonds_;
  char chain_;
  std::string type_;
};


//! Assign unique indices to all atoms in the monomer, starting at optional start value
/*!
  This returns the next unused index. 
*/
inline int index_atoms(const Heterogen &m, int start=0) {
  for (Heterogen::Atom_const_iterator it= m.atoms_begin();
       it != m.atoms_end(); ++it) {
    it->atom().set_index(Atom::Index(start++));
  }
  return start;
}


CGAL_SWAP(Heterogen);
CGAL_OUTPUT(Heterogen);

CGAL_PDB_END_NAMESPACE
#endif
