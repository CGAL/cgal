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
#include <CGAL/PDB/small_map.h>
#include <iostream>
#include <vector>
#include <CGAL/PDB/internal/Nested_iterator.h>


CGAL_PDB_BEGIN_NAMESPACE

//! A class for representing a whole PDB file with possibly several models.
/*!  See pdb_split.cpp for an example manipulating a PDB by splitting
  it into parts.
*/
class PDB {
  CGAL_SMALL_MAP_VALUE_TYPE(Model_vt, Label<PDB>, Model, model);
  typedef small_map<Model_vt> Models;
public:
  //! Read a pdb file from the stream
  /*!  The optional bool controls whether errors (such as unparsable
    PDB lines). Set it to false to disable printing errors.
  */
  PDB(std::istream &in, bool print_errors=false);
  //! Construct a empty PDB
  PDB();
  ~PDB();

  typedef Label<PDB> Model_key;

  //! Write a pdb file to the stream
  std::ostream& write(std::ostream &out) const;

 
  
  //! Return the number of models. 
  CGAL_SIZE(models, return models_.size());

  //! add a model with an automatically chosen number
  Model_key push_back(const Model &m);

  //! check if there are no models
  bool empty() const {return models_.empty();}

  void swap_with(PDB &o);
 
  //! Set the header
  template <class It>
  void set_header(It b, It e){
    header_.clear();
    header_.insert(header_.end(), b,e);
  }

  //! A class which identifies a Chain in the PDB
  class Chain_key {
    Model_key model_;
    Model::Chain_key chain_;
  public:
    Chain_key(){}
    Chain_key(Model_key m, Model::Chain_key c): model_(m), chain_(c){}
    Model_key model() const {return model_;}
    Model::Chain_key chain() const {return chain_;}
  };


  //! An iterator through the unparsed std::string lines of the header of the PDB.
  CGAL_CONST_ITERATOR(Header, header, 
		      std::vector<std::string>::const_iterator,
		      return header_.begin(),
		      return header_.end());

  //! An iterator through the Model objects in the PDB
  CGAL_ITERATOR(Model, model, Models::iterator, 
		return models_.begin(),
		return models_.end());

  //! An iterator through the Model objects in the PDB
  CGAL_CONST_ITERATOR(Model, model, Models::const_iterator, 
		return models_.begin(),
		return models_.end());

  //! Find a Model with the given key, return models_end() if none is found
  CGAL_FIND(Model, return models_.find(k));

  //! Add a model (or change an existing one).
  CGAL_INSERT(Model, return models_.insert(Models::value_type(k,m)));


  class Chain_iterator_value_type {
      Chain_key index_;
      Chain *chain_;
    public:
      Chain_iterator_value_type(Chain_key f, Chain* s): index_(f), chain_(s){}
      Chain_key key() const {return index_;}
      Chain &chain() const {return *chain_;}
      Chain_iterator_value_type():chain_(NULL){}
    };
  class Chain_const_iterator_value_type {
      Chain_key index_;
      const Chain *chain_;
    public:
      Chain_const_iterator_value_type(Chain_key f, const Chain* s): index_(f), chain_(s){}
      Chain_key key() const {return index_;}
      const Chain &chain() const {return *chain_;}
      Chain_const_iterator_value_type():chain_(NULL){}
    };
protected:
  //! \cond
  struct Iterator_traits {
    typedef Model_iterator Outer_it;
    typedef Model::Chain_iterator Inner_it;
    typedef Chain_iterator_value_type value_type;
    struct Inner_range{
      std::pair<Inner_it, Inner_it> operator()(Outer_it it) const {
	return std::make_pair(it->model().chains_begin(), it->model().chains_end());
      }
    };
    struct Make_value{
      value_type operator()(Outer_it oit, Inner_it iit) const {
	return value_type(Chain_key(oit->key(), iit->key()), &iit->chain());
      }
    };
  };


  struct Iterator_const_traits {
    typedef Model_const_iterator Outer_it;
    typedef Model::Chain_const_iterator Inner_it;
    typedef Chain_const_iterator_value_type value_type;
    struct Inner_range{
      std::pair<Inner_it, Inner_it> operator()(Outer_it it) const {
	return std::make_pair(it->model().chains_begin(), it->model().chains_end());
      }
    };
    struct Make_value{
      value_type operator()(Outer_it oit, Inner_it iit) const {
	return value_type(Chain_key(oit->key(), iit->key()), &iit->chain());
      }
    };
  };
  //! \endcond
public:
 
  //! An iterator through the CGAL::PDB::Chain objects contained in the PDB.
  CGAL_ITERATOR(Chain, chain,
		internal::Nested_iterator<Iterator_traits>,  
		return Chain_iterator(models_begin(), models_end()),
		return Chain_iterator(models_begin(), models_end()));
  
  //! An iterator through the CGAL::PDB::Chain objects contained in the PDB.
  CGAL_CONST_ITERATOR(Chain, chain,
		      internal::Nested_iterator<Iterator_const_traits>,  
		      return Chain_const_iterator(models_begin(), models_end()),
		      return Chain_const_iterator(models_begin(), models_end()));


 

private:
  void load(std::istream &in, bool print_errors);
  void build_heterogens();

  std::vector<std::string> header_;
  Models models_;
  std::vector<std::pair<int, int> > connections_;
};


//! Assign unique indices to all atoms in the PDB, starting at optional start value
/*!
  This returns the next unused index. 
*/
inline int index_atoms(const PDB &c, int start=0) {
  for (PDB::Model_const_iterator it= c.models_begin(); it != c.models_end(); ++it) {
    start= index_atoms(it->model(), start);
  }
  return start;
}

CGAL_OUTPUT(PDB);

CGAL_SWAP(PDB);

CGAL_PDB_END_NAMESPACE
#endif
