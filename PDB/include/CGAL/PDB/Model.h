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
#include <CGAL/PDB/Chain.h>
#include <CGAL/PDB/Heterogen.h>
#include <vector>
#include <string>
#include <boost/tuple/tuple.hpp>

namespace CGAL { namespace PDB {

class PDB;

class Heterogen_key {
  typedef Heterogen_key This;
  std::string name_;
  int num_;
public:
  Heterogen_key(std::string name,
                int num): name_(name), num_(num) {
  }
  Heterogen_key(): num_(-1){}
  CGAL_COMPARISONS2(name_, num_);
  const std::string &name() const {
    return name_;
  }
  int number() const {return num_;}
};

//! A class representing a single model from a PDB file.
/*!
  You can iterator through the chains and soon the heterogens. 
*/
class Model {
  typedef Model This;
  friend class PDB;
public:
  CGAL_SMALL_MAP_VALUE_TYPE(Chain_pair, CGAL::Label<Model>, Chain, chain); 
  CGAL_SMALL_MAP_VALUE_TYPE(Heterogen_pair, CGAL::PDB::Heterogen_key,
                            Heterogen, heterogen); 
private:
  typedef small_map<Chain_pair> ChainsMap; 
  typedef small_map<Heterogen_pair> HeterogensMap; 
public:
  //! Construct an empty model
  Model();

  typedef CGAL::PDB::Heterogen_key Heterogen_key;
  typedef CGAL::Label<Model> Chain_key;

  void swap_with(Model &o);

  //! write to a pdb
  void write(int model_index, std::ostream &out) const;
  
  //! write screen
  std::ostream& write(std::ostream &out) const{
    write(0,out);
    return out;
  }


  //! Iterator through the chains
  CGAL_ITERATOR(Chain, chain, ChainsMap::const_iterator, ChainsMap::iterator, 
                chains_.begin(), chains_.end());


  //! get the chain identified by k
  CGAL_FIND(Chain, chains_.find(k), chains_.end());

  CGAL_INSERT(Chain, chains_.insert(ChainsMap::value_type(k, m)););

  //! An iterator through the heterogen objects.
  CGAL_ITERATOR(Heterogen, heterogen, 
                HeterogensMap::const_iterator, HeterogensMap::iterator,
                heterogens_.begin(),
                heterogens_.end());

  //! get the chain identified by k
  CGAL_FIND(Heterogen, heterogens_.find(k), heterogens_.end());

  CGAL_INSERT(Heterogen, heterogens_.insert(HeterogensMap::value_type(k, m)););


  //! A unique identified of an atom in the Model
  struct Atom_key: public boost::tuple<Chain_key, Chain::Monomer_key, 
				       Monomer::Atom_key>{
    typedef boost::tuple<Chain_key, Chain::Monomer_key, Monomer::Atom_key> P;
    Atom_key(Chain_key ck, Chain::Monomer_key mk, Monomer::Atom_key at): P(ck, mk, at){}
    Atom_key(){}
    operator Chain_key() const {
      return chain_key();
    }
    operator Monomer::Atom_key() const {
      return atom_key();
    }
    operator Chain::Monomer_key() const {
      return monomer_key();
    }
    Monomer::Atom_key atom_key() const {
      return P::get<2>();
    }
    Chain::Monomer_key monomer_key() const {
      return P::get<1>();
    }
    Chain_key chain_key() const {
      return P::get<0>();
    }
  };



  //! The endpoint of a bond
  class Bond_endpoint {
    const Atom *atom_;
    Atom_key ak_;
    friend class Bond_it;
 public:
    Bond_endpoint(Atom_key ak, const Atom *atom): atom_(atom), ak_(ak){}    
    Bond_endpoint():atom_(NULL){}
    const Atom &atom() const {return *atom_;}
    Atom_key key() const {return ak_;}
  };


  //! A chemical bond within the protein
  typedef std::pair<Bond_endpoint, Bond_endpoint> Bond; 
  class Atom_iterator_value_type {
      Atom_key index_;
      Atom *atom_;
  public:
      Atom_iterator_value_type(Atom_key f, Atom* s): index_(f), atom_(s){}
      Atom_key key() const {return index_;}
      Atom &atom() const {return *atom_;}
      Atom_iterator_value_type():atom_(NULL){}
    };
  class Atom_const_iterator_value_type {
      Atom_key index_;
      const Atom *atom_;
  public:
      Atom_const_iterator_value_type(Atom_key f, const Atom* s): index_(f), atom_(s){}
      Atom_key key() const {return index_;}
      const Atom &atom() const {return *atom_;}
      Atom_const_iterator_value_type():atom_(NULL){}
    };
protected:
  //! \cond
 struct Iterator_traits {
   typedef Chains Outer;
   typedef Chain::Atoms Inner;
   typedef Atom_iterator_value_type value_type;
    struct Get_inner {
      Inner operator()(Outer::iterator it) const {
	return it->chain().atoms();
      }
    };
    struct Make_value{
      value_type operator()(Outer::iterator oit, Inner::iterator iit) const {
	return value_type(Atom_key(oit->key(), iit->key().monomer_key(),
				   iit->key().atom_key()), &iit->atom());
      }
    };
  };


  struct Iterator_const_traits {
    typedef Chain_consts Outer;
    typedef Chain::Atom_consts Inner;
    typedef Atom_const_iterator_value_type value_type;
   struct Get_inner{
     Inner operator()(Outer::iterator it) const {
       return it->chain().atoms();
      }
    };
    struct Make_value{
      value_type operator()(Outer::iterator oit, Inner::iterator iit) const {
	return value_type(Atom_key(oit->key(), iit->key().monomer_key(),
				   iit->key().atom_key()), &iit->atom());
      }
    };
  };
  //! \endcond
public:
 //! An iterator to iterate through all the atoms of the protein  
  CGAL_ITERATOR(Atom, atom, 
                internal::Nested_iterator<Iterator_const_traits >,
                internal::Nested_iterator<Iterator_traits >,
                boost::make_iterator_range(chains_.begin(), chains_.end()),
                boost::make_iterator_range(chains_.end(), chains_.end()));

  //! \cond
  class Bond_it {
    friend class Model;
  public:
    //! The value_type is a CGAL::PDB::Chain::Bond
    typedef Bond value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef const Bond& reference;
    typedef const Bond* pointer;
    
    reference operator*() const {
      return ret_;
    }
    pointer operator->() const {
      return &ret_;
    }
    Bond_it operator++() {
      ++ait_;
      while (ait_== rit_->chain().bonds().end()) {
	++rit_;
	if (rit_!= rend_) {
	  ait_= rit_->chain().bonds().begin();
	} else {
	  return *this;
	}
      }
      make_bond();
      return *this;
    }
    bool operator==(const Bond_it& o) const {
      if (rit_ == rend_) return rit_==o.rit_;
      else return rit_== o.rit_ && ait_ == o.ait_;
    }
  bool operator!=(const Bond_it& o) const {
      return !operator==(o);
    }
    
    Bond_it(){}

    CGAL_COPY_CONSTRUCTOR(Bond_it);

    Bond_it(Chain_consts r): rit_(r.begin()), rend_(r.end()){
      if (!r.empty()) {
	ait_= rit_->chain().bonds().begin();
	make_bond();
      }
    }

  
  protected:
    void copy_from(const Bond_it &o) {
      rit_= o.rit_;
      rend_= o.rend_;
      if (rit_!= rend_) {
	ait_= o.ait_;
	ret_= o.ret_;
      }
    }

   
    void make_bond() {
      ret_= Bond(Bond_endpoint(Atom_key(rit_->key(), 
					ait_->first.key().monomer_key(),
					ait_->first.key().atom_key()),
			       &ait_->first.atom()),
		 Bond_endpoint(Atom_key(rit_->key(), 
					ait_->second.key().monomer_key(),
					ait_->second.key().atom_key()),
			       &ait_->second.atom())); 
    }

    Chain_consts::iterator rit_, rend_;
    Chain::Bonds::iterator ait_;
    Bond ret_;
  };
  //! \endcond

  CGAL_CONST_ITERATOR(Bond, bond, Bond_it,
                      boost::make_iterator_range(chains_.begin(),
                                                 chains_.end()),
                      boost::make_iterator_range(chains_.end(), 
                                                 chains_.end()));
 

private:
  void process_line(const char *c);
  void process_atom(const char *c);
  void process_hetatom(const char *c);
  void add_hetatom(int numscan, int snum, char name[], char alt, 
                   char resname[], char chain, int resnum,
                   char insertion_residue_code, 
                   float x, float y, float z,
                   float occupancy, float tempFactor,
                   char segID[], char element[], char charge[]);

  std::vector<std::string> extra_;
  ChainsMap chains_;
  HeterogensMap heterogens_;
};

CGAL_SWAP(Model)

CGAL_OUTPUT(Model)

//! Assign unique indices to all atoms in the Model, starting at optional start value
/*!
  This returns the next unused index. 
*/
inline int index_atoms(const Model &c, int start=0) {
  CGAL_PDB_FOREACH(std::iterator_traits<Model::Chain_consts::iterator>::reference cc, c.chains()) {
    start= index_atoms(cc.chain(), start);
  }
  return start;
}

}}
#endif
