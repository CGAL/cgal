// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_GMAP_CELL_ITERATORS_H
#define CGAL_GMAP_CELL_ITERATORS_H 1

#include <CGAL/GMap_dart_iterators.h>
#include <CGAL/Cell_iterators.h>

// TODO do all the orbit iterator of any orbit ?

namespace CGAL {

  /** @file Cell_iterators.h
   * Cell iterators. There are 3 classes:
   * - GMap_cell_iterator<Map,i,dim>: one dart per each i-cell
   * - GMap_one_dart_per_incident_cell_iterator<Map,Ite,i,dim>
   * - GMap_one_dart_per_cell_iterator<Map,Ite,i,dim>
   * - one specialisation of the CMap_cell_iterator for the
   *    GMap_dart_iterator_basic_of_all iterator
   */

  //****************************************************************************
  /* Class CMap_cell_iterator<Map,GMap_dart_iterator_basic_of_all<Map>,
     i,dim,Tag_false>: specialization to iterate onto
     * all the cells of the gmap.
     */
  template <typename Map_,unsigned int i,unsigned int dim,bool Const>
  class CMap_cell_iterator<Map_,GMap_dart_iterator_basic_of_all<Map_,Const>,
                           i,dim,Const,Tag_false>:
    public GMap_dart_iterator_basic_of_all<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_all<Map_,Const> Base;
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
  /* Class GMap_cell_iterator<Map,i,dim,Tag_false>: to iterate onto
     * all the cells of the gmap.
     */
  template <typename Map_,unsigned int i,unsigned int dim,bool Const>
  class GMap_cell_iterator: public GMap_dart_iterator_basic_of_all<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_all<Map_,Const> Base;
    typedef GMap_cell_iterator<Map_,i,dim,Const> Self;

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
    GMap_cell_iterator(Map& amap):
      Base(amap),
      mmark_number(amap.get_new_mark())
    {
      CGAL_static_assertion( (boost::is_same<typename Base::Basic_iterator,
                              Tag_true>::value) );
      CGAL_assertion(amap.is_whole_map_unmarked(mmark_number));
      mark_cell<Map,i,dim>(amap, (*this), mmark_number);
    }

   /// Constructor with a dart in parameter (for end iterator).
    GMap_cell_iterator(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mmark_number(amap.get_new_mark())
    {
      if (adart!=this->mmap->null_handle)
        mark_cell<Map,i,dim>(amap, (*this), mmark_number);
    }

    /// Destructor.
    ~GMap_cell_iterator()
    {
      if (this->mmap->get_number_of_times_mark_reserved(mmark_number)==1)
        unmark_treated_darts();
      this->mmap->free_mark(mmark_number);
      this->mmark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
    }

    /// Copy constructor.
    GMap_cell_iterator(const Self& aiterator):
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
  /* Class GMap_one_dart_per_incident_cell_iterator<Map,i,j,dim>: to iterate
   * onto one dart per i-cell incident to the given j-cell.
   */
  template <typename Map_,unsigned int i,unsigned int j,
            unsigned int dim=Map_::dimension,bool Const=false>
  class GMap_one_dart_per_incident_cell_iterator:
    public CMap_cell_iterator<Map_,
                              GMap_dart_iterator_basic_of_cell
                              <Map_,j,dim,Const>, i,dim,Const>
  {
  public:
    typedef GMap_one_dart_per_incident_cell_iterator<Map_,i,j,dim,Const> Self;
    typedef CMap_cell_iterator<Map_,
                               GMap_dart_iterator_basic_of_cell<Map_,j,
                                                                dim,Const>,
                               i,dim,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;
    typedef Tag_false Basic_iterator;

    /// Main constructor.
    GMap_one_dart_per_incident_cell_iterator(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  /* Class GMap_one_dart_per_cell_iterator<Map,i,dim>: to iterate onto the
   * i-cells of the map (one dart by each i-cell).
   */
  template <typename Map_,unsigned int i,unsigned int dim=Map_::dimension,
            bool Const=false>
  class GMap_one_dart_per_cell_iterator:
    public CMap_cell_iterator<Map_,GMap_dart_iterator_basic_of_all<Map_,Const>,
                              i,dim,Const>
  {
  public:
    typedef GMap_one_dart_per_cell_iterator<Map_,i,dim,Const> Self;
    typedef CMap_cell_iterator<Map_,
                               GMap_dart_iterator_basic_of_all<Map_,Const>,
                               i,dim,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;
    typedef Tag_false Basic_iterator;

    /// Main constructor.
    GMap_one_dart_per_cell_iterator(Map& amap): Base(amap)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    GMap_one_dart_per_cell_iterator(Map& amap, Dart_handle adart):
        Base(amap, adart)
    {}
  };
//****************************************************************************
//****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_GMAP_CELL_ITERATORS_H
//******************************************************************************
