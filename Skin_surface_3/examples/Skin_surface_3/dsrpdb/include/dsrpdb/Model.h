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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_MODEL_H_
#define DSR_MODEL_H_
#include <dsrpdb/Protein.h>
#include <vector>
#include <string>

namespace dsrpdb {

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

    typedef std::vector<Protein>::const_iterator Const_chains_iterator;
    Const_chains_iterator chains_begin() const {
      return chains_.begin();
    }
    Const_chains_iterator chains_end() const {
      return chains_.end();
    } 

    //! add a strand
    void new_chain(const Protein &p);

    //! write to a pdb
    void write(std::ostream &out) const;

    //! return the index or -1 if it is not valid
    int index() const{return index_;}

    //! set the index 
    void set_index(int ind) {index_=ind;}
  protected:
    void process_line(const char *c);

    std::vector<std::string> extra_;
    std::vector<Protein> chains_;
    int index_;
  };
}
#endif
