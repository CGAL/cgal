// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_ITERATORS_BASE_HH
#define CGAL_COMBINATORIAL_MAP_ITERATORS_BASE_HH 1

#include <CGAL/disable_warnings.h>

#include <CGAL/Compact_container.h>
#include <queue>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {

  /** @file Combinatorial_map_iterators_base.h
   * Basic classes that serve as tools for definition of iterators.
   There are 3 classes:
   *  - CMap_dart_iterator<Map,Const> is the basic generic class defining
   *    what is an interator on darts.
   *  - CMap_extend_iterator<Map,Ite,Bi> to extend the given iterator by adding
   *    the involution Bi.
   *  - CMap_non_basic_iterator<Map_,Ite> to transform the basic iterator Ite
   *    into the corresponding non basic iterator.
   *  - CMap_range<Map,It,ConstIt,BasicIt> generic definition of a range
   *    given an iterator and its const version
   *  - CMap_const_range<Map,It,ConstIt,BasicIt> generic definition of a const
   *    range given an iterator and its const version
   */
  //****************************************************************************
  /// OperationState: type to keep the last operation used by the previous ++.
  typedef char OperationState;

  /// Enum of all the possible operations used by the ++ operator.
  enum
  {
    OP_NONE = -1, ///< Beginning of the iterator (there is not yet operator++).
    OP_BETAI,     ///< Previous op was the first beta.
    OP_BETAI_INV, ///< Previous op was the inverse of the first beta.
    OP_BETAJ,     ///< Previous op was the second beta.
    OP_BETAK,     ///< Previous op was the third beta.
    OP_BETA0I,    ///< Previous op was beta0 o the first beta.
    OP_BETAI1,    ///< Previous op was the first beta o beta1.
    OP_BETAIJ,    ///< Previous op was the composition of two beta.
    OP_BETAJI,    ///< Previous op was the composition of two beta.
    OP_BETA21,     ///< Previous op was beta21.
    OP_JUMP,      ///< Previous op was a jump .
    OP_POP,       ///< Previous op pop a dart from a stack or a queue.
    OP_END        ///< Previous op go out of the iterator.
  };
  //****************************************************************************
  /** Generic class of iterator onto darts.
   * Class CMap_dart_iterator is a generic iterator. All the combinatorial
   * map iterator classes inherit from this class (or one of its subclass).
   */
  template < typename Map_,bool Const=false >
  class CMap_dart_iterator: public boost::mpl::if_c< Const,
      typename Map_::Dart_container::const_iterator,
      typename Map_::Dart_container::iterator>::type
    //public internal::CC_iterator<typename Map_::Dart_container,Const>
  {
  public:
    typedef CMap_dart_iterator<Map_,Const> Self;

    typedef typename boost::mpl::if_c< Const,
          typename Map_::Dart_container::const_iterator,
          typename Map_::Dart_container::iterator>::type Base;

    typedef typename boost::mpl::if_c< Const,
                                       typename Map_::Dart_const_handle,
                                       typename Map_::Dart_handle>::type
                                     Dart_handle;
    typedef typename boost::mpl::if_c< Const, const Map_,
                                       Map_>::type Map;

    typedef std::input_iterator_tag iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::difference_type difference_type;
    typedef typename Base::pointer pointer;
    typedef typename Base::reference reference;

    /// true iff this iterator is basic
    typedef Tag_true Basic_iterator;

    typedef typename Map::size_type size_type;

  public:
    /// Main constructor.
    CMap_dart_iterator(Map& amap, Dart_handle adart):
      Base(adart),
      mmap(&amap),
      mfirst_dart(adart),
      mprev_op(OP_NONE)
    {}

    /// == operator.
    bool operator==(const Self& aiterator) const
    {
      return ( ((*this==mmap->null_handle) && (aiterator==mmap->null_handle)) ||
               (mfirst_dart == aiterator.mfirst_dart &&
                static_cast<const Base&>(*this)==
                static_cast<const Base&>(aiterator)) );
    }

    /// != operator.
    bool operator!=(const Self& aiterator) const
    { return !operator==(aiterator); }

    /// Accessor to the initial dart of the iterator.
    Dart_handle get_first_dart() const { return mfirst_dart; }

    /// Accessor to the combinatorial map.
    Map* get_combinatorial_map() const { return mmap; }

    /// Rewind of the iterator to its beginning.
    void rewind()
    { set_current_dart(mfirst_dart); mprev_op = OP_NONE; }

    /// Test if the iterator is at its end.
    bool cont() const
    { return static_cast<const Base&>(*this)!=mmap->null_handle; }

    /// Get the previous operation used for the last ++.
    OperationState prev_operation()  const { return mprev_op; }

  protected:
    /// Set the current dart to a given dart
    void set_current_dart(Dart_handle adart)
    { Base::operator=(adart); }

  private:
    /// operator -- in private to invalidate the base operator.
    Self& operator--()
    { return *this; }
    /// operator -- in private to invalidate the base operator.
    Self operator--(int)
    { return *this; }

  protected:
    /// test if adart->beta(ai) exists and is not marked for amark
    bool is_unmarked(Dart_handle adart, unsigned int ai, size_type amark) const
    { return
        !mmap->is_marked(mmap->beta(adart,ai), amark);
    }

    /// test if adart->beta(ai)->beta(aj) exists
    bool exist_betaij(Dart_handle adart, unsigned int ai, unsigned int aj) const
    { return !mmap->is_free(adart, ai) && !mmap->is_free(mmap->beta(adart, ai), aj); }

    /// test if adart->beta(ai)->beta(aj) exists and is not marked for amark
    bool is_unmarked2(Dart_handle adart, unsigned int ai, unsigned int aj,
                      typename Map::size_type amark) const
    { return
        !mmap->is_marked(mmap->beta(adart, ai, aj), amark);
    }

  protected:
    /// The map containing the darts to iterate on.
    Map* mmap;

    /// The initial dart of the iterator.
    Dart_handle mfirst_dart;

    /// The last operation used for the ++ operator.
    OperationState mprev_op;
  };
  //****************************************************************************
  /* Class CMap_extend_iterator<Map,Ite,Bi> which extend a given iterator by
   * adding Bi and by using a stack and a mark.
   * General case when Ite does not have already a stack.
   */
  template <typename Map_,typename Ite,int Bi,
            typename Ite_has_stack=typename Ite::Use_mark>
  class CMap_extend_iterator: public Ite
  {
  public:
    typedef CMap_extend_iterator<Map_,Ite,Bi, Ite_has_stack> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

    typedef typename Map::size_type size_type;

    CGAL_static_assertion( (Bi<=Map::dimension &&
                            boost::is_same<Ite_has_stack,Tag_false>::value) );

  public:
    /// Main constructor.
    CMap_extend_iterator(Map& amap, Dart_handle adart, size_type amark):
      Base(amap, adart),
      mmark_number(amark),
      minitial_dart(adart)
    {
      if ( minitial_dart!=amap.null_handle )
      {
        this->mmap->mark_null_dart(mmark_number);
        this->mmap->mark(minitial_dart, mmark_number);
        if (!this->mmap->is_free(minitial_dart, Bi) &&
            this->mmap->beta(minitial_dart, Bi)!=minitial_dart )
        {
          mto_treat.push(this->mmap->beta(minitial_dart, Bi));
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
      this->mmap->mark_null_dart(mmark_number);
      if (!this->mmap->is_free(minitial_dart, Bi) &&
          this->mmap->beta(minitial_dart, Bi)!=minitial_dart)
      {
        mto_treat.push(this->mmap->beta(minitial_dart, Bi));
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

          if (
               !this->mmap->is_marked(this->mmap->beta(*this, Bi), mmark_number) )
          {
            mto_treat.push(this->mmap->beta(*this, Bi));
          }
        }
      }
      else
      {
        this->mmap->mark((*this), mmark_number);
        if (
             !this->mmap->is_marked(this->mmap->beta(*this, Bi), mmark_number) )
        {
          mto_treat.push(this->mmap->beta(*this, Bi));
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
  /* Class CMap_extend_iterator<Map,Ite,Bi> which extend a given iterator by
   * adding Bi and by using a stack and a mark.
   * Specialization when Ite has already a stack.
   */
  template <typename Map_,typename Ite,int Bi>
  class CMap_extend_iterator<Map_,Ite,Bi,Tag_true>: public Ite
  {
  public:
    typedef CMap_extend_iterator<Map_,Ite,Bi,Tag_true> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

    typedef typename Map::size_type size_type;

    /// Main constructor.
    CMap_extend_iterator(Map& amap, Dart_handle adart, size_type amark):
      Base(amap, adart, amark)
    {
      if ( this->minitial_dart!=amap.null_handle &&
           !this->mmap->is_free(this->minitial_dart, Bi) &&
           this->mmap->beta(this->minitial_dart, Bi)!=this->minitial_dart )
      {
        this->mto_treat.push(this->mmap->beta(this->minitial_dart, Bi));
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != Map::INVALID_MARK);
      Base::rewind();
      if ( !this->mmap->is_free(this->minitial_dart, Bi) &&
           this->mmap->beta(this->minitial_dart, Bi)!=this->minitial_dart )
      {
        this->mto_treat.push(this->mmap->beta(this->minitial_dart, Bi));
      }
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      Base::operator++();

      if ( this->cont() )
      {
        CGAL_assertion( this->mmap->is_marked(*this, this->mmark_number) );

        if (
            !this->mmap->is_marked(this->mmap->beta(*this, Bi), this->mmark_number) )
        {
          this->mto_treat.push(this->mmap->beta(*this, Bi));
        }
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  //* Class CMap_non_basic_iterator allows to transform a basic_iterator onto
  //* a non basic one, depending if the basic iterator uses mark or not.
  template <typename Map_,typename Basic_iterator,
            typename Use_mark=typename Basic_iterator::Use_mark>
  class CMap_non_basic_iterator;
  //****************************************************************************
  template <typename Map_,typename Base_>
  class CMap_non_basic_iterator<Map_,Base_,Tag_true>:
    public Base_
  {
  public:
    typedef CMap_non_basic_iterator<Map_,Base_,Tag_true> Self;
    typedef Base_ Base;

    typedef typename Base::Map Map;
    typedef typename Base::Dart_handle Dart_handle;

    /// True iff this iterator is basic
    typedef Tag_false Basic_iterator;

    typedef typename Map::size_type size_type;

    CGAL_static_assertion( (boost::is_same<typename Base::Basic_iterator,
                            Tag_true>::value) );

    /// Main constructor.
    CMap_non_basic_iterator(Map& amap, Dart_handle adart1):
      Base(amap, adart1, amap.get_new_mark())
    {}

    /// Destructor.
    ~CMap_non_basic_iterator() CGAL_NOEXCEPT(CGAL_NO_ASSERTIONS_BOOL)
    {
      CGAL_destructor_assertion( this->mmark_number!=Map::INVALID_MARK );
      if (this->mmap->get_number_of_times_mark_reserved
          (this->mmark_number)==1)
        unmark_treated_darts();
      this->mmap->free_mark(this->mmark_number);
      this->mmark_number = Map::INVALID_MARK; // To avoid basic class to try to unmark darts.
    }

    /// Copy constructor.
    CMap_non_basic_iterator(const Self& aiterator):
      Base(aiterator)
    { this->mmap->share_a_mark(this->mmark_number); }

    /// Assignment operator.
    Self& operator=(const Self& aiterator)
    {
      if (this != &aiterator)
      {
        if (this->mmap->get_number_of_times_mark_reserved
            (this->mmark_number)==1)
          unmark_treated_darts();
        this->mmap->free_mark(this->mmark_number);
        this->mmark_number = Map::INVALID_MARK;

        Base::operator=(aiterator);
        this->mmap->share_a_mark(this->mmark_number);
      }
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      unmark_treated_darts();
      Base::rewind();
    }

    using Base::operator++;

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Unmark all the marked darts during the iterator.
    void unmark_treated_darts()
    {
      if (this->mmap->is_whole_map_unmarked(this->mmark_number)) return;

      this->mmap->negate_mark(this->mmark_number);

      if (this->mmap->is_whole_map_unmarked(this->mmark_number)) return;

      Base::rewind();
      while (this->mmap->number_of_unmarked_darts(this->mmark_number) > 0)
        this->operator++();
      this->mmap->negate_mark(this->mmark_number);
      CGAL_assertion(this->mmap->is_whole_map_unmarked(this->mmark_number));
    }
  };
  //****************************************************************************
  template <typename Map_,typename Base_>
  class CMap_non_basic_iterator<Map_,Base_,Tag_false>:
    public Base_
  {
  public:
    typedef CMap_non_basic_iterator<Map_,Base_,Tag_false> Self;
    typedef Base_ Base;

    typedef typename Base::Map Map;
    typedef typename Base::Dart_handle Dart_handle;

    /// True iff this iterator is basic
    typedef Tag_false Basic_iterator;

    CGAL_static_assertion( (boost::is_same<typename Base::Basic_iterator,
                            Tag_true>::value) );

    /// Main constructor.
    CMap_non_basic_iterator(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  template <typename Map_, typename It, typename Const_it,
            typename Basic_iterator=typename It::Basic_iterator>
  struct CMap_range
  {
    typedef It iterator;
    typedef Const_it const_iterator;
    CMap_range(Map_ &amap, typename Map_::Dart_handle adart) :
      mmap(amap), mdart(adart), msize(0)
    {}
    iterator begin()             { return iterator(mmap,mdart); }
    iterator end()               { return iterator(mmap,mmap.null_handle); }
    const_iterator begin() const { return const_iterator(mmap,mdart); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
    typename Map_::size_type size() const
    {
      if (msize==0)
        for ( const_iterator it=begin(), itend=end(); it!=itend; ++it)
          ++msize;
      return msize;
    }
    bool empty() const
    { return mmap.is_empty(); }
  private:
    Map_ & mmap;
    typename Map_::Dart_handle mdart;
    mutable typename Map_::size_type msize;
  };
  //****************************************************************************
  template <typename Map_, typename It, typename Const_it>
  struct CMap_range<Map_,It,Const_it,Tag_true>
  {
    typedef CMap_range<Map_,It,Const_it,Tag_true> Base_cmap_range;
    typedef It iterator;
    typedef Const_it const_iterator;

    typedef typename Map_::size_type size_type;

    CMap_range(Map_ &amap, typename Map_::Dart_handle adart, size_type amark=Map_::INVALID_MARK):
      mmap(amap), mdart(adart), msize(0), mmark(amark)
    {}
    iterator begin()             { return iterator(mmap,mdart,mmark); }
    iterator end()               { return iterator(mmap,mmap.null_handle,mmark); }
    const_iterator begin() const { return const_iterator(mmap,mdart,mmark); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle,mmark); }
    typename Map_::size_type size() const
    {
      if (msize==0)
        for ( CMap_non_basic_iterator<Map_,const_iterator> it(mmap,mdart);
              it.cont(); ++it )
          ++msize;
      return msize;
    }
    bool empty() const
    { return mmap.is_empty(); }
  private:
    Map_ & mmap;
    typename Map_::Dart_handle mdart;
    mutable typename Map_::size_type msize;
    size_type mmark;
  };
  //****************************************************************************
  template <typename Map_, typename Const_it,
            typename Basic_iterator=typename Const_it::Basic_iterator>
  struct CMap_const_range
  {
    typedef Const_it const_iterator;
    CMap_const_range(const Map_ &amap, typename Map_::Dart_const_handle adart):
      mmap(amap), mdart(adart), msize(0)
    {}
    const_iterator begin() const { return const_iterator(mmap,mdart); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
    typename Map_::size_type size() const
    {
      if (msize==0)
        for ( const_iterator it=begin(), itend=end(); it!=itend; ++it)
          ++msize;
      return msize;
    }
    bool empty() const
    { return mmap.is_empty(); }
  private:
    const Map_ & mmap;
    typename Map_::Dart_const_handle mdart;
    mutable typename Map_::size_type msize;
  };
  //****************************************************************************
  template <typename Map_, typename Const_it>
  struct CMap_const_range<Map_,Const_it,Tag_true>
  {
    typedef Const_it const_iterator;
    typedef typename Map_::size_type size_type;
    CMap_const_range(const Map_ &amap, typename Map_::Dart_const_handle adart,
                     size_type amark=Map_::INVALID_MARK):
      mmap(amap), mdart(adart), msize(0), mmark(amark)
    {}
    const_iterator begin() const { return const_iterator(mmap,mdart,mmark); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle,mmark); }
    typename Map_::size_type size() const
    {
      if (msize==0)
        for ( CMap_non_basic_iterator<Map_,const_iterator> it(mmap,mdart);
              it.cont(); ++it )
          ++msize;
      return msize;
    }
    bool empty() const
    { return mmap.is_empty(); }
  private:
    const Map_ & mmap;
    typename Map_::Dart_const_handle mdart;
    mutable typename Map_::size_type msize;
    size_type mmark;
  };
  //****************************************************************************
} // namespace CGAL

#include <CGAL/enable_warnings.h>

//******************************************************************************
#endif // CGAL_COMBINATORIAL_MAP_ITERATORS_BASE_HH
//******************************************************************************
