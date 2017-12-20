// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_CELL_ITERATORS_H
#define CGAL_CELL_ITERATORS_H 1

#include <CGAL/Dart_iterators.h>
#include <CGAL/Combinatorial_map_basic_operations.h>

#include <boost/type_traits/is_same.hpp>

// TODO do all the orbit iterator of any orbit ?

namespace CGAL {

  /** @file Cell_iterators.h
   * All the cell iterators.
   * There are 3 classes:
   *  - CMap_cell_iterator<Map,Ite,i,dim>: "tools" class used for the
   *    two following iterators.
   * - CMap_one_dart_per_incident_cell_iterator<Map,Ite,i,dim>
   * - CMap_one_dart_per_cell_iterator<Map,Ite,i,dim>
   */

  //****************************************************************************
  template <typename Map_,typename Ite, unsigned int i,
            unsigned int dim=Map_::dimension,bool Const=false,
            typename Use_mark=typename Ite::Use_mark>
  class CMap_cell_iterator;
  //****************************************************************************
  /* Class CMap_cell_iterator<Map,Ite,i,dim,Tag_true>: to iterate onto the
   * cells incident to the given iterator which uses mark.
   */
  template <typename Map_,typename Ite,unsigned int i,
            unsigned int dim,bool Const>
  class CMap_cell_iterator<Map_,Ite,i,dim,Const,Tag_true>: public Ite
  {
  public:
    typedef CMap_cell_iterator<Map_,Ite,i,dim,Const,Tag_true> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

  protected:
    /// Unmark all the marked darts during the iterator.
    void unmark_treated_darts()
    {
      if (this->mmap->is_whole_map_unmarked(mcell_mark_number)) return;

      this->mmap->negate_mark(this->mmark_number);
      this->mmap->negate_mark(mcell_mark_number);

      Ite::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this),
                           mcell_mark_number);

      while (this->mmap->number_of_unmarked_darts(mcell_mark_number) > 0 ||
             this->mmap->number_of_unmarked_darts(this->mmark_number) > 0)
      {
        this->operator++();
      }

      this->mmap->negate_mark(mcell_mark_number);
      this->mmap->negate_mark(this->mmark_number);

      CGAL_assertion(this->mmap->is_whole_map_unmarked(this->mmark_number));
      CGAL_assertion(this->mmap->is_whole_map_unmarked(mcell_mark_number));
    }

  public:
    /// Main constructor.
    CMap_cell_iterator(Map& amap, Dart_handle adart):
      Ite(amap, adart, amap.get_new_mark()),
      mcell_mark_number(amap.get_new_mark())
    {
      CGAL_static_assertion( (boost::is_same<typename Ite::Basic_iterator,
                              Tag_true>::value) );
      CGAL_assertion(amap.is_whole_map_unmarked(mcell_mark_number));

      mark_cell<Map,i,dim>(amap, adart, mcell_mark_number);
    }

    /// Destructor.
    ~CMap_cell_iterator()
    {
      if (this->mmap->get_number_of_times_mark_reserved(mcell_mark_number)==1)
        unmark_treated_darts();
      this->mmap->free_mark(mcell_mark_number);
      this->mmap->free_mark(this->mmark_number);
      this->mcell_mark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
      this->mmark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
    }

    /// Copy constructor.
    CMap_cell_iterator(const Self& aiterator):
      Ite(aiterator),
      mcell_mark_number(aiterator.mcell_mark_number)
    {
      this->mmap->share_a_mark(this->mmark_number);
      this->mmap->share_a_mark(this->mcell_mark_number);
    }

    /// Assignment operator.
    Self& operator=(const Self& aiterator)
    {
      if (this != &aiterator)
      {
        Ite::operator=(aiterator);
        this->mmap->share_a_mark(this->mmark_number);
        this->mmap->share_a_mark(mcell_mark_number);
      }
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      unmark_treated_darts();
      Ite::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this), mcell_mark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      do
      {
        Ite::operator++();
      }
      while (this->cont() &&
             this->mmap->is_marked((*this), mcell_mark_number));

      if (this->cont())
      {
        mark_cell<Map,i,dim>(*this->mmap, (*this),
                             mcell_mark_number);
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  private:
    /// A mark used to mark treated cells.
    typename Map::size_type mcell_mark_number;
  };
  //****************************************************************************
  /* Class CMap_cell_iterator<Map,Ite,i,dim,Tag_false>: to iterate onto the
   * cells incident to the given iterator which does not use mark.
   */
  template <typename Map_,typename Ite, unsigned int i,
            unsigned int dim,bool Const>
  class CMap_cell_iterator<Map_,Ite,i,dim,Const,Tag_false>: public Ite
  {
  public:
    typedef CMap_cell_iterator<Map_,Ite,i,dim,Const,Tag_false> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

  protected:
    /// Unmark all the marked darts during the iterator.
    void unmark_treated_darts()
    {
      if (this->mmap->is_whole_map_unmarked(mmark_number)) return;

      this->mmap->negate_mark(mmark_number);

      if (this->mmap->is_whole_map_unmarked(mmark_number)) return;

      Ite::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this),
                           mmark_number);
      while (this->mmap->number_of_unmarked_darts(mmark_number) > 0)
        this->operator++();
      this->mmap->negate_mark(mmark_number);
      CGAL_assertion(this->mmap->is_whole_map_unmarked(mmark_number));
    }

  public:
    /// Main constructor.
    CMap_cell_iterator(Map& amap, Dart_handle adart):
      Ite(amap, adart),
      mmark_number(amap.get_new_mark())
    {
      CGAL_static_assertion( (boost::is_same<typename Ite::Basic_iterator,
                              Tag_true>::value) );
      CGAL_assertion(amap.is_whole_map_unmarked(mmark_number));
      mark_cell<Map,i,dim>(amap, adart, mmark_number);
    }

    /// Destructor.
    ~CMap_cell_iterator()
    {
      if (this->mmap->get_number_of_times_mark_reserved(mmark_number)==1)
        unmark_treated_darts();
      this->mmap->free_mark(mmark_number);
      this->mmark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
    }

    /// Copy constructor.
    CMap_cell_iterator(const Self& aiterator):
      Ite(aiterator),
      mmark_number(aiterator.mmark_number)
    { this->mmap->share_a_mark(mmark_number); }

    /// Assignment operator.
    Self& operator=(const Self & aiterator)
    {
      if (this != &aiterator)
      {
        Ite::operator=(aiterator);
        mmark_number = aiterator.mmark_number;
        this->mmap->share_a_mark(mmark_number);
      }
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      unmark_treated_darts();
      Ite::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this), mmark_number);
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      do
      {
        Ite::operator++();
      }
      while (this->cont() &&
             this->mmap->is_marked((*this), mmark_number));

      if (this->cont())
        mark_cell<Map,i,dim>(*this->mmap, (*this), mmark_number);
      return *this;
    }

  private:
    /// A mark used to mark treated cells.
    typename Map::size_type mmark_number;
  };
  //****************************************************************************
  /* Class CMap_cell_iterator<Map,CMap_dart_iterator_basic_of_all<Map>,
     i,dim,Tag_false>: specialization to iterate onto
     * all the cells of the map.
     */
  template <typename Map_,unsigned int i,unsigned int dim,bool Const>
  class CMap_cell_iterator<Map_,CMap_dart_iterator_basic_of_all<Map_,Const>,
                           i,dim,Const,Tag_false>:
    public CMap_dart_iterator_basic_of_all<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_all<Map_,Const> Base;
    typedef CMap_cell_iterator<Map_,Base,i,dim,Const,Tag_false> Self;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

  protected:
    /// Unmark all the marked darts during the iterator.
    void unmark_treated_darts()
    {
      if (this->mmap->is_whole_map_unmarked(mmark_number)) return;

      this->mmap->negate_mark(mmark_number);

      if (this->mmap->is_whole_map_unmarked(mmark_number)) return;

      Base::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this),
                           mmark_number);
      while (this->mmap->number_of_unmarked_darts(mmark_number) > 0)
        this->operator++();
      this->mmap->negate_mark(mmark_number);
      CGAL_assertion(this->mmap->is_whole_map_unmarked(mmark_number));
    }

  public:
    /// Main constructor.
    CMap_cell_iterator(Map& amap):
      Base(amap),
      mmark_number(amap.get_new_mark())
    {
      CGAL_static_assertion( (boost::is_same<typename Base::Basic_iterator,
                              Tag_true>::value) );
      CGAL_assertion(amap.is_whole_map_unmarked(mmark_number));
      mark_cell<Map,i,dim>(amap, (*this), mmark_number);
    }

   /// Constructor with a dart in parameter (for end iterator).
    CMap_cell_iterator(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mmark_number(amap.get_new_mark())
    {
      if (adart!=this->mmap->null_handle)
        mark_cell<Map,i,dim>(amap, (*this), mmark_number);
    }

    /// Destructor.
    ~CMap_cell_iterator()
    {
      if (this->mmap->get_number_of_times_mark_reserved(mmark_number)==1)
        unmark_treated_darts();
      this->mmap->free_mark(mmark_number);
      this->mmark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
    }

    /// Copy constructor.
    CMap_cell_iterator(const Self& aiterator):
      Base(aiterator),
      mmark_number(aiterator.mmark_number)
    { this->mmap->share_a_mark(mmark_number); }

    /// Assignment operator.
    Self& operator=(const Self& aiterator)
    {
      if (this != &aiterator)
      {
        Base::operator=(aiterator);
        mmark_number = aiterator.mmark_number;
        this->mmap->share_a_mark(mmark_number);
      }
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      unmark_treated_darts();
      Base::rewind();
      mark_cell<Map,i,dim>(*this->mmap, (*this), mmark_number);
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      do
      {
        Base::operator++();
      }
      while (this->cont() &&
             this->mmap->is_marked((*this), mmark_number));

      if (this->cont())
        mark_cell<Map,i,dim>(*this->mmap, (*this), mmark_number);
      return *this;
    }

  private:
    /// A mark used to mark treated cells.
    typename Map::size_type mmark_number;
  };
  //****************************************************************************
  /* Class CMap_one_dart_per_incident_cell_iterator<Map,i,j,dim>: to iterate
   * onto one dart per i-cell incident to the given j-cell.
   */
  template <typename Map_,unsigned int i,unsigned int j,
            unsigned int dim=Map_::dimension,bool Const=false>
  class CMap_one_dart_per_incident_cell_iterator:
    public CMap_cell_iterator<Map_,
                              CMap_dart_iterator_basic_of_cell
                              <Map_,j,dim,Const>, i,dim,Const>
  {
  public:
    typedef CMap_one_dart_per_incident_cell_iterator<Map_,i,j,dim,Const> Self;
    typedef CMap_cell_iterator<Map_,
                               CMap_dart_iterator_basic_of_cell<Map_,j,
                                                                dim,Const>,
                               i,dim,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;
    typedef Tag_false Basic_iterator;

    /// Main constructor.
    CMap_one_dart_per_incident_cell_iterator(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  /* Class CMap_one_dart_per_cell_iterator<Map,i,dim>: to iterate onto the
   * i-cells incident of the map (one dart by each i-cell).
   */
  template <typename Map_,unsigned int i,unsigned int dim=Map_::dimension,
            bool Const=false>
  class CMap_one_dart_per_cell_iterator:
    public CMap_cell_iterator<Map_,CMap_dart_iterator_basic_of_all<Map_,Const>,
                              i,dim,Const>
  {
  public:
    typedef CMap_one_dart_per_cell_iterator<Map_,i,dim,Const> Self;
    typedef CMap_cell_iterator<Map_,
                               CMap_dart_iterator_basic_of_all<Map_,Const>,
                               i,dim,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;
    typedef Tag_false Basic_iterator;

    /// Main constructor.
    CMap_one_dart_per_cell_iterator(Map& amap): Base(amap)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    CMap_one_dart_per_cell_iterator(Map& amap, Dart_handle adart):
	     Base(amap, adart)
    {}
  };
//****************************************************************************
//****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_CELL_ITERATORS_H
//******************************************************************************
