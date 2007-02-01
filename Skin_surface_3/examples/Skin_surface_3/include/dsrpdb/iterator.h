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

#ifndef DSR_PDB_ITERATOR_H
#define DSR_PDB_ITERATOR_H

#include <dsrpdb/Point.h>
#include <dsrpdb/Protein.h>

namespace dsrpdb {
  //! Take an atom label, atom pair and return the coordinates
  struct Atom_coordinates {
    typedef Point result_type;
    typedef const std::pair<Residue::Atom_label, Atom>& argument_type;
    const result_type& operator()(argument_type a) const {
      return a.second.cartesian_coords();
    }
    
  };

  struct Atom_index {
    typedef Atom::Index result_type;
    typedef const std::pair<Residue::Atom_label, Atom>& argument_type;
    const result_type& operator()(argument_type a) const {
      static result_type rt;
      rt= a.second.index();
      return rt;
    }
    
  };


  //! Take an iterator which returns atoms and return their coordinates instead
  template <class It, class Projector>
  class Projection_iterator {
    typedef Projection_iterator<It, Projector> This;
  public:
    typedef typename Projector::result_type value_type;
    typedef typename It::iterator_category iterator_category;
    typedef typename It::difference_type difference_type;
    typedef const typename Projector::result_type& reference;
    typedef const typename Projector::result_type* pointer;

    reference operator*() {
      //std::cout << ait_->second.index() << std::endl;
      return p_(*cur_);
    }
    pointer operator->() {
      //std::cout << ait_->second.index() << std::endl;
      return &p_(*cur_);
    }
    This operator++(int) {
      This t=*this;
      ++cur_;
      return t;
    }
    This operator++() {
      ++cur_;
      return *this;
    }
    
    bool operator==(const This& o) const {
      return cur_ == o.cur_;
    }
    bool operator!=(const This& o) const {
      return cur_!= o.cur_;
    }

    Projection_iterator(){}
    Projection_iterator(It c, Projector p): cur_(c), p_(p){}
    
  protected:
    It cur_;
    Projector p_;
  };





  //! Create an atom coordinates iterator
  template <class It, class Proj>
  Projection_iterator<It, Proj> projection_iterator(It a, Proj proj) {
    return Projection_iterator<It, Proj>(a, proj);
  }

  

  //! Filters the iterator output
  template <class It, class Filter>
  class Filter_iterator {
    typedef Filter_iterator<It, Filter> This;
  public:
    typedef typename It::value_type value_type;
    typedef typename It::iterator_category iterator_category;
    typedef typename It::difference_type difference_type;
    typedef typename It::reference reference;
    typedef typename It::pointer pointer;

    reference operator*() {
      return cur_.operator*();
    }

    pointer operator->() {
      return cur_.operator->();
    }

    This operator++(int) {
      This t=*this;
      ++cur_;
      fix();
      return t;
    }
    This operator++() {
      ++cur_;
      fix();
      return *this;
    }
    
    bool operator==(const This& o) const {
      return cur_ == o.cur_;
    }
    bool operator!=(const This& o) const {
      return cur_!= o.cur_;
    }

    Filter_iterator(){}
    Filter_iterator(It c, It e, Filter f): cur_(c), end_(e), f_(f){fix();}
    
  protected:
    void fix(){
      while (cur_ != end_ && !f_(*cur_)) ++cur_;
    }
    It cur_, end_;
    Filter f_;
  };

  //! Create a filter iterator by passing the functor and the iterator.
  /*!
    See Yes, Is_backbone and Is_CA for example filters.
  */
  template <class It, class F>
  Filter_iterator<It, F> filter_iterator(It b, It e, F f){
    return Filter_iterator<It, F>(b,e,f);
  }

  //! Functor that returns true.
  struct Yes {
    template <class A>
    bool operator()(A) {return true;}
  };
  
  //! Return true if an atom is a backbone atoms.
  struct Is_backbone {
    template <class A>
    bool operator()(const A &a) {
      return a.first == Residue::AL_C || a.first == Residue::AL_CA 
	|| a.first == Residue::AL_N;
    }
  };

  //! Returns true if an atom is a CA atom.
  struct Is_CA {
    template <class A>
    bool operator()(const A &a) {
      return a.first == Residue::AL_CA;
    }
  };

  //! An iterator for iterating through the coordinates of backbone atoms of a protein
  typedef Projection_iterator<Filter_iterator<Protein::Const_atoms_iterator, Is_backbone>, Atom_coordinates>  Protein_backbone_coordinates_iterator;
  
  //! a begin iterator for iterating through the coordinates of backbone atoms of a protein
  /*!
    This is a good example of how to use the various iterator classes.
  */
  inline Protein_backbone_coordinates_iterator backbone_coordinates_begin(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_begin(), p.atoms_end(), Is_backbone()), Atom_coordinates());
  }

  //! an end iterator for iterating through the coordinates of backbone atoms of a protein
  inline Protein_backbone_coordinates_iterator  backbone_coordinates_end(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_end(), p.atoms_end(), Is_backbone()), Atom_coordinates());
  }


  //! An iterator for iterating through the coordinates of C_alpha atoms of a protein
  typedef Projection_iterator<Filter_iterator<Protein::Const_atoms_iterator, Is_CA>, Atom_coordinates >  Ca_backbone_coordinates_iterator;
  
  //! a begin iterator for iterating through the coordinates of C_alpha atoms of a protein
  inline Ca_backbone_coordinates_iterator ca_coordinates_begin(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_begin(), p.atoms_end(), Is_CA()), Atom_coordinates());
  }

  //! an end iterator for iterating through the coordinates of C_alpha atoms of a protein
  inline Ca_backbone_coordinates_iterator  ca_coordinates_end(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_end(), p.atoms_end(), Is_CA()), Atom_coordinates());
  }


  //! An iterator for iterating through the coordinates of all atoms of a protein
  typedef Projection_iterator<Filter_iterator<Protein::Const_atoms_iterator, Yes>, Atom_coordinates>  All_coordinates_iterator;
  
  //! a begin iterator for iterating through the coordinates of all the atoms of a protein
  inline All_coordinates_iterator all_coordinates_begin(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_begin(), p.atoms_end(), Yes()), Atom_coordinates());
  }

  //! an end iterator for iterating through the coordinates of all the atoms of a protein
  inline All_coordinates_iterator  all_coordinates_end(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_end(), p.atoms_end(), Yes()), Atom_coordinates());
  }


  //! An iterator through the indices of the backbone atoms
  typedef Projection_iterator<Filter_iterator<Protein::Const_atoms_iterator, Is_backbone>, Atom_index> Protein_backbone_indices_iterator;

  //! Begin iterating through the indices of backbone atoms
  inline Protein_backbone_indices_iterator protein_backbone_indices_begin(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_begin(), p.atoms_end(), Is_backbone()), Atom_index());
  }

  //! an end iterator for iterating through the indices of all the backbone atoms of a protein
  inline Protein_backbone_indices_iterator  protein_backbone_indices_end(const Protein &p) {
    return projection_iterator(filter_iterator(p.atoms_end(), p.atoms_end(), Is_backbone()), Atom_index());
  }

}

#endif
