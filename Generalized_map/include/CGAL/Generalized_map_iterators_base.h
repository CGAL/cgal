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
#ifndef CGAL_GENERALIZED_MAP_ITERATORS_BASE_HH
#define CGAL_GENERALIZED_MAP_ITERATORS_BASE_HH 1

 // to get OperationState type, some OP, and CMap_dart_iterator
#include <CGAL/Combinatorial_map_iterators_base.h>

// Other includes
#include <CGAL/Compact_container.h>
#include <queue>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {

  /** @file generalized_map_iterators_base.h
   * Basic classes that serve as tools for definition of iterators.
   * There is 1 class:
   *  - GMap_extend_iterator<GMap,Ite,Ai> to extend the given iterator by
   *    adding the involution Ai.
   */
    //****************************************************************************
  /// Enum of all the possible operations used by the ++ operator.
  /// Extension of the enum for combinatorial maps.
  enum
    {
    OP_ALPHAI=OP_END+1, ///< Previous op was the first alpha.
    OP_ALPHAJ,     ///< Previous op was the second alpha.
    OP_ALPHAK,     ///< Previous op was the third alpha.
    OP_ALPHAIJ,    ///< Previous op was the composition of two alpha.
    OP_ALPHAJI,    ///< Previous op was the composition of two alpha.
    };
  //****************************************************************************
  /* Class GMap_extend_iterator<Map,Ite,Ai> which extend a given iterator by
   * adding Ai and by using a stack and a mark.
   * General case when Ite does not have already a stack.
   */
  template <typename Map_,typename Ite,int Ai,
            typename Ite_has_stack=typename Ite::Use_mark>
  class GMap_extend_iterator: public Ite
  {
  public:
    typedef GMap_extend_iterator<Map_,Ite,Ai, Ite_has_stack> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Base::size_type size_type;

    typedef Tag_true Use_mark;

    CGAL_static_assertion( (Ai<=Map::dimension &&
                            boost::is_same<Ite_has_stack,Tag_false>::value) );

  public:
    /// Main constructor.
    GMap_extend_iterator(Map& amap, Dart_handle adart, size_type amark):
      Base(amap, adart),
      mmark_number(amark),
      minitial_dart(adart)
    {
      if ( adart!=amap.null_handle )
      {
        this->mmap->mark(adart, mmark_number);
        if (!this->mmap->template is_free<Ai>(adart))
        {
          mto_treat.push(this->mmap->template alpha<Ai>(adart));
        }
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      Base::operator= ( Base(*this->mmap,minitial_dart) );
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark(minitial_dart, mmark_number);
      if (!this->mmap->template is_free<Ai>(minitial_dart))
      {
        mto_treat.push(this->mmap->template alpha<Ai>(minitial_dart));
      }
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      CGAL_assertion(this->cont());

      do
      {
        Base::operator++();
      }
      while ( this->cont() &&
              this->mmap->is_marked(*this, mmark_number) );

      if ( !this->cont() )
      {
        while ( !mto_treat.empty() &&
                this->mmap->is_marked(mto_treat.front(), mmark_number))
        {
          mto_treat.pop();
        }

        if ( !mto_treat.empty() )
        {
          Base::operator= ( Base(*this->mmap,mto_treat.front()) );
          mto_treat.pop();
          this->mprev_op = OP_POP;
          CGAL_assertion( !this->mmap->is_marked((*this), mmark_number) );
          this->mmap->mark((*this), mmark_number);

          if (!this->mmap->template is_free<Ai>(*this) &&
              !this->mmap->is_marked(this->mmap->template alpha<Ai>(*this), mmark_number) )
          {
            mto_treat.push(this->mmap->template alpha<Ai>(*this));
          }
        }
      }
      else
      {
        this->mmap->mark((*this), mmark_number);
        if (!this->mmap->template is_free<Ai>(*this) &&
            !this->mmap->is_marked(this->mmap->alpha(*this, Ai), mmark_number) )
        {
          mto_treat.push(this->mmap->template alpha<Ai>(*this));
        }
      }

      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    size_type mmark_number;

    /// Initial dart
    Dart_handle minitial_dart;
  };
  //****************************************************************************
  /* Class GMap_extend_iterator<Map,Ite,Ai> which extend a given iterator by
   * adding Ai and by using a stack and a mark.
   * Specialization when Ite has already a stack.
   */
  template <typename Map_,typename Ite,int Ai>
  class GMap_extend_iterator<Map_,Ite,Ai,Tag_true>: public Ite
  {
  public:
    typedef GMap_extend_iterator<Map_,Ite,Ai,Tag_true> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

    typedef typename Map::size_type size_type;

    /// Main constructor.
    GMap_extend_iterator(Map& amap, Dart_handle adart, size_type amark):
      Base(amap, adart, amark)
    {
      if (adart!=amap.null_handle)
      {
        if (!this->mmap->is_free(adart, Ai) &&
            !this->mmap->is_marked(this->mmap->alpha(amap, Ai), this->mmark_number))
        {
          this->mto_treat.push(adart->alpha(Ai));
        }
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != Map::INVALID_MARK);
      Base::rewind();
      if ( !this->mmap->is_free(this->minitial_dart, Ai) )
      {
        this->mto_treat.push(this->mmap->alpha(this->minitial_dart, Ai));
      }
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      Base::operator++();

      if ( this->cont() )
      {
        CGAL_assertion( this->mmap->is_marked(*this, this->mmark_number) );

        if (!this->mmap->is_free(*this, Ai) &&
            !this->mmap->is_marked(this->mmap->alpha(*this, Ai),
                                   this->mmark_number))
        {
          this->mto_treat.push(this->mmap->alpha(*this, Ai));
        }
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_GENERALIZED_MAP_ITERATORS_BASE_HH
//******************************************************************************
