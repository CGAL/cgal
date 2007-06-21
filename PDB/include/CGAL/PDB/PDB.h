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

#ifndef CGAL_DSR_PDB_H
#define CGAL_DSR_PDB_H

#include <CGAL/auto_link/PDB.h>

#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Model.h>
#include <iostream>
#include <vector>

CGAL_PDB_BEGIN_NAMESPACE

  //! A class for representing a whole PDB file with possibly several models.
  /*!  See pdb_split.cc for an example manipulating a PDB by splitting
    it into parts.
  */
  class PDB {
  public:
    //! Read a pdb file from the stream
    /*!  The optional bool controls whether errors (such as unparsable
      PDB lines). Set it to false to disable printing errors.
    */
    PDB(std::istream &in, bool print_errors=false);
    //! Construct a empty PDB
    PDB();
    ~PDB();

    //! Write a pdb file to the stream
    void write(std::ostream &out) const;

    //! get the ith model
    const Model &model(unsigned int i) const;
    //! get the ith model
    Model &model(unsigned int i);
    //! set a model
    void new_model(const Model &m);
    //! how many models are there?
    size_t number_of_models() const;

    typedef std::vector<std::string>::const_iterator Header_iterator;
    //! iterator through the lines of the header
    Header_iterator header_begin() const;
    //! end iterator
    Header_iterator header_end() const;

    //! Set the header
    template <class It>
    void set_header(It b, It e){
      header_.clear();
      header_.insert(header_.end(), b,e);
    }

    //! iterator for models
    typedef std::vector<Model>::const_iterator Const_models_iterator;
    //! non-const iterator for models
    typedef std::vector<Model>::iterator Models_iterator;
    
    Const_models_iterator models_begin() const {
      return models_.begin();
    }
    Const_models_iterator models_end() const {
      return models_.end();
    }
    Models_iterator models_begin() {
      return models_.begin();
    }
    Models_iterator models_end() {
      return models_.end();
    }
  protected:
    void load(std::istream &in, bool print_errors);

    std::vector<std::string> header_;
    std::vector<Model> models_;
  };

CGAL_PDB_END_NAMESPACE
#endif
