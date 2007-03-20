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

#ifndef DSR_MODEL_H_
#define DSR_MODEL_H_
#include <CGAL/PDB/Protein.h>
#include <vector>
#include <string>

CGAL_PDB_BEGIN_NAMESPACE

  class PDB;

  //! A class representing a single model from a PDB file.
  /*!
    You can iterator through the chains and soon the heterogens. 
  */
  class Model {
    friend class PDB;
  public:
    //! Construct an empty model
    Model();
    //! Construct an empty model with a number
    Model(unsigned int i);
    
    //! The number of strands
    size_t number_of_chains() const;
    //! get the ith strand
    Protein &chain(unsigned int i);
    //! get the ith strand
    const Protein &chain(unsigned int i) const;    

    //! add a strand
    void new_chain(const Protein &p);

    //! write to a pdb
    void write(std::ostream &out) const;

    //! return the index or -1 if it is not valid
    int index() const{return index_;}

    //! set the index 
    void set_index(int ind) {index_=ind;}

    typedef std::vector<Protein>::const_iterator Const_chains_iterator;
    //! Begin iterating through the chains.
    Const_chains_iterator chains_begin() const {
      return chains_.begin();
    }
    //! End iterating through the chains. 
    Const_chains_iterator chains_end() const {
      return chains_.end();
    }

    class Hetatom_data {
    public:
      Hetatom_data(const char *rnm, 
		   const char *anm, int rn, char ch): resname_(rnm),
						      atomname_(anm),
						      rnum_(rn), chain_(ch){
      }
      const char *molecule_name() const {
	return resname_.c_str();
      }
      const char *atom_name() const {
	return atomname_.c_str();
      }
      int molecule_number() const {
	return rnum_;
      }
      char chain() const {
	return chain_;
      }
    protected:
      std::string resname_;
      std::string atomname_;
      int rnum_;
      char chain_;
      
    };

    //! An iterator through CGAL::PDB::Atom values for the HETATM records.
    typedef std::vector<std::pair<Hetatom_data, Atom> >::const_iterator Const_hetatoms_iterator;
    //! Begin iterating through CGAL:PDB::Atom values for the HETATM records.
    Const_hetatoms_iterator hetatoms_begin() const {
      return hetatoms_.begin();
    }
    Const_hetatoms_iterator hetatoms_end() const {
      return hetatoms_.end();
    }

    //! The number of hetatoms
    unsigned int hetatoms_size() const {
      return hetatoms_.size();
    }
  protected:
    void process_line(const char *c);

    std::vector<std::string> extra_;
    std::vector<Protein> chains_;
    std::vector<std::pair<Hetatom_data, Atom> > hetatoms_;
    int index_;
  };

CGAL_PDB_END_NAMESPACE
#endif
