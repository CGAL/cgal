// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#ifndef CGAL_DART_ITERATORS_HH
#define CGAL_DART_ITERATORS_HH 1

#include <queue>

namespace CGAL {

  /** @file Dart_iterators.h
   * Definition of dart iterators. There are 9 iterators:
   *  - CMap_dart_iterator_basic_of_orbit<Map,Beta...>
   *  - CMap_dart_iterator_basic_of_cell<Map,i,d>
   *  - CMap_dart_iterator_basic_of_all
   *  - CMap_dart_iterator_of_orbit<Map,Beta...>
   *  - CMap_dart_iterator_of_cell<Map,i,d>
   *  - CMap_dart_iterator_basic_of_involution<Map,i,d>
   *  - CMap_dart_iterator_of_involution<Map,i,d>
   *  - CMap_dart_iterator_basic_of_involution_inv<Map,i,d>
   *  - CMap_dart_iterator_of_involution_inv<Map,i,d>
   * but many specializations to optimize specific cases.
   *  
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
  //****************************************************************************
  /** Generic class of iterator onto darts.
   * Class CMap_dart_iterator is a pure virtual generic iterator. This
   * class defines what is an iterator. All the iterator classes inherit
   * from this class.
   */
  template < typename Map_,bool Const=false >
  class CMap_dart_iterator: 
    public internal::CC_iterator<typename Map_::Dart_container,Const>
  {
  public:
    typedef CMap_dart_iterator<Map_,Const> Self;
    typedef internal::CC_iterator<typename Map_::Dart_container,Const> Base;

    typedef Base Dart_handle;
    typedef typename boost::mpl::if_c< Const, const Map_,
                                       Map_>::type Map;

    typedef std::input_iterator_tag iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::difference_type difference_type;
    typedef typename Base::pointer pointer;
    typedef typename Base::reference reference;    

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
      return ( ((*this==NULL) && (aiterator==NULL)) ||
	       (mfirst_dart == aiterator.mfirst_dart && 
		((const Base&)*this==(const Base&)aiterator)) );
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
    bool cont() const { return *this != NULL; }

    /// Get the previous operation used for the last ++.
    OperationState prev_operation()  const { return mprev_op; }

    /// Return true iff this iterator is basic
    static bool is_basic_iterator()
    { return true; }
    
  protected:
    /// Set the current dart to a given dart
    void set_current_dart(Dart_handle adart)
    { Base::operator=(adart); }

  private:
    /// operator -- in private to invalidate the base operator.
    Self& operator--()
    { return *this; }
    /// operator -- in private to invalidate the base operator.
    void operator--(int)
    {}

  protected:
    /// test if adart->beta(ai) exists and is not marked for amark
    bool is_unmarked(Dart_handle adart, unsigned int ai, unsigned amark) const
    { return !adart->is_free(ai) &&					
	!mmap->is_marked(adart->beta(ai), amark); }

    /// test if adart->beta(ai)->beta(aj) exists
    bool exist_betaij(Dart_handle adart, unsigned int ai, unsigned int aj) const
    { return !adart->is_free(ai) && !adart->beta(ai)->is_free(aj); }

    /// test if adart->beta(ai)->beta(aj) exists and is not marked for amark
    bool is_unmarked2(Dart_handle adart, unsigned int ai, unsigned int aj,
		      unsigned amark) const
    { return exist_betaij(adart,ai,aj) &&
	!mmap->is_marked(adart->beta(ai)->beta(aj), amark); }

  protected:
    /// The map containing the darts to iterate on.
    Map* mmap;

    /// The initial dart of the iterator.
    Dart_handle mfirst_dart;

    /// The last operation used for the ++ operator.
    OperationState mprev_op;
  };
  //****************************************************************************
  //**********************BASIC ITERATORS***************************************
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map, Beta...>: to iterate
   * on the darts of the orbit <Beta...>
   */
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template<typename Map,bool Const,int... Beta>
  class CMap_dart_iterator_basic_of_orbit_generic;
  #else
  template <typename Map,bool Const,int B1=-1,int B2=-1,int B3=-1,int B4=-1,
	    int B5=-1,int B6=-1,int B7=-1,int B8=-1,int B9=-1>
  class CMap_dart_iterator_basic_of_orbit_generic;

  template <typename Map,bool Const,int B1=-1,int B2=-1,int B3=-1,int B4=-1,
	    int B5=-1,int B6=-1,int B7=-1,int B8=-1,int B9=-1>
  struct Get_CMap_dart_iterator_basic_of_orbit;

  template<typename Map,bool Const,int B1,int B2,int B3,int B4,int B5,int B6,
	   int B7,int B8,int B9>
  struct Get_CMap_dart_iterator_basic_of_orbit
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2,B3,B4,
						     B5,B6,B7,B8,B9> type;
  };

  template<typename Map,bool Const,int B1,int B2,int B3,int B4,int B5,int B6,
	   int B7,int B8>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,B4,B5,B6,B7,B8,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2,B3,B4,B5,
						     B6,B7,B8> type;
  };

  template<typename Map,bool Const,int B1,int B2,int B3,int B4,int B5,int B6,
	   int B7>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,B4,B5,B6,B7,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,
					     B1,B2,B3,B4,B5,B6,B7> type;
  };

  template<typename Map,bool Const,int B1,int B2,int B3,int B4,int B5,int B6>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,B4,B5,B6,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2,B3,B4,
						     B5,B6> type;
  };

  template<typename Map,bool Const,int B1,int B2,int B3,int B4,int B5>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,B4,B5,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,B1,B2,B3,B4,
						     B5,Const> type;
  };

  template<typename Map,bool Const,int B1,int B2,int B3,int B4>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,B4,-1,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2,B3,
						     B4> type;
  };

  template<typename Map, int B1,int B2,int B3,bool Const>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,B3,-1,-1,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2,B3> type;
  };

  template<typename Map, int B1,int B2,bool Const>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,B2,-1,-1,-1,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1,B2> type;
  };

  template<typename Map, int B1,bool Const>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      B1,-1,-1,-1,-1,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const,B1> type;
  };

  template<typename Map,bool Const>
  struct Get_CMap_dart_iterator_basic_of_orbit<Map,Const,
					      -1,-1,-1,-1,-1,-1,-1,-1,-1>
  {
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,Const> type;
  };
  
  template <typename Map_,bool Const,int B1,int B2,int B3,int B4,int B5,int B6,
	    int B7,int B8,int B9>
  class CMap_dart_iterator_basic_of_orbit_generic: 
    public Get_CMap_dart_iterator_basic_of_orbit<Map_,Const,B2,B3,B4,B5,
						B6,B7,B8,B9>::type
  {
  public:
    typedef typename Get_CMap_dart_iterator_basic_of_orbit<Map_,Const,B1,B2,B3,
							  B4,B5,B6,B7,B8,
							  B9>::type Self;
    typedef typename Get_CMap_dart_iterator_basic_of_orbit<Map_,Const,B2,B3,B4,
							  B5,B6,B7,B8,
							  B9>::type Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
					     int amark):
      Base(amap, adart, amark)
    { 
      CGAL_assertion( B1>=0 && B1<=Map::dimension );

      if (adart!=NULL)
	{
	  if (!adart->is_free(B1) &&
	      !this->mmap->is_marked(adart->beta(B1), this->mmark_number))
	    this->mto_treat.push(adart->beta(B1));
	}
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != -1);
      Base::rewind();
      if (!(*this)->is_free(B1) &&
	  !this->mmap->is_marked((*this)->beta(B1), this->mmark_number))
	this->mto_treat.push((*this)->beta(B1));
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      Base::operator++();

      if (this->cont())
	{
	  if (!(*this)->is_free(B1) &&
	      !this->mmap->is_marked((*this)->beta(B1), 
				     this->mmark_number))
	    this->mto_treat.push((*this)->beta(B1));
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };  
  #endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  //****************************************************************************
  // Case when Beta... is empty: iterator of self
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int /*amark*/):
      Base(amap, adart)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      this->set_current_dart(NULL);
      this->mprev_op = OP_END;
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,0>: iterate onto orbit <beta0>.
   * Begin by turning around the facet with beta0, then turn if
   * necessary in the second direction by using beta1.
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int /*amark*/):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Assignment operator.
    Self& operator= (const Self & aiterator)
    {
      if (this != &aiterator)
	{
	  Base::operator=(aiterator);
	  mfirst_dir = aiterator.mfirst_dir;
	}
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mfirst_dir && (*this)->is_free(0))
	{
	  this->set_current_dart(this->mfirst_dart);
	  mfirst_dir = false;
	  this->mprev_op = OP_JUMP;
	}
      else
	{
	  this->mprev_op = OP_BETAI;
	}

      if (mfirst_dir)
	{
	  CGAL_assertion(!(*this)->is_free(0));
	  this->set_current_dart((*this)->beta(0));

	  if ((*this)==this->mfirst_dart)
	    {
	      this->set_current_dart(NULL);
	      this->mprev_op = OP_END;
	    }
	}
      else
	{
	  if ((*this)->is_free(1))
	    {
	      this->set_current_dart(NULL);
	      this->mprev_op = OP_END;
	    }
	  else
	    {
	      this->set_current_dart((*this)->beta(1));
	      this->mprev_op = OP_BETAI_INV;
	    }
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Boolean: true iff we turn in the first direction (i.e. using beta0).
    bool mfirst_dir;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,1>: iterate onto orbit <beta1>.
   * Begin by turning around the facet with beta1, then turn if
   * necessary in the second direction by using beta0.
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int /*amark*/):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mfirst_dir && (*this)->is_free(1))
	{
	  this->set_current_dart(this->mfirst_dart);
	  mfirst_dir = false;
	  this->mprev_op = OP_JUMP;
	}
      else
	{
	  this->mprev_op = OP_BETAI;
	}

      if (mfirst_dir)
	{
	  CGAL_assertion(!(*this)->is_free(1));
	  this->set_current_dart((*this)->beta(1));

	  if ((*this)==this->mfirst_dart)
	    {
	      this->set_current_dart(NULL);
	      this->mprev_op = OP_END;
	    }
	}
      else
	{
	  if ((*this)->is_free(0))
	    {
	      this->set_current_dart(NULL);
	      this->mprev_op = OP_END;
	    }
	  else
	    {
	      this->set_current_dart((*this)->beta(0));
	      this->mprev_op = OP_BETAI_INV;
	    }
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Boolean: true iff we turn in the first direction (i.e. using beta0).
    bool mfirst_dir;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Bi>: to iterate
   * on the darts of the orbit <Bi> (2<=Bi<=dimension)
   * (not for beta0 and beta1 which are special cases).
   */
  template <typename Map_,bool Const,int Bi>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart)
    { CGAL_assertion( Bi>=2 && Bi<=Map::dimension ); }

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int /*amark*/):
      Base(amap, adart)
    { CGAL_assertion( Bi>=2 && Bi<=Map::dimension ); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if ((*this)!=this->mfirst_dart || (*this)->is_free(Bi))
	{
	  this->set_current_dart(NULL);
	  this->mprev_op = OP_END;
	}
      else
	{
	  this->set_current_dart((*this)->beta(Bi));
	  this->mprev_op = OP_BETAI;
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  /* Class CMap_extend_iterator<Map,Ite,Bi> which extend a given iterator by 
   * adding Bi and by using a stack and a mark.
   */
  template <typename Map_,typename Ite,int Bi>
  class CMap_extend_iterator: public Ite
  {
  public:
    typedef CMap_extend_iterator<Map_,Ite,Bi> Self;
    typedef Ite Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_extend_iterator(Map& amap, Dart_handle adart, int amark):
      Base(amap, adart),
      mmark_number(amark),
      minitial_dart(adart)
    {
      if ( adart!=NULL )
	{
	  this->mmap->mark(adart, mmark_number);
	  if (!(*this)->is_free(Bi))
	    mto_treat.push((*this)->beta(Bi));
	}
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::operator= ( Base(*this->mmap,minitial_dart) );
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      if (!(*this)->is_free(Bi))
	mto_treat.push((*this)->beta(Bi));
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Base::operator++();

      if ( !this->cont() )
	{
	  if ( !mto_treat.empty() )
	    {
	      Dart_handle res=NULL;
	      do
		{
		  res = mto_treat.front();
		  mto_treat.pop();
		}
	      while (!mto_treat.empty() &&
		     this->mmap->is_marked(res, mmark_number));

	      if (!this->mmap->is_marked(res, mmark_number))
		{
		  Base::operator= ( Base(*this->mmap,res) );
		  this->mprev_op = OP_POP;
		}
	    }
	}

      if ( this->cont() )
	{
	  CGAL_assertion( !this->mmap->is_marked((*this), 
						 mmark_number) );
	  this->mmap->mark((*this), mmark_number);
	  
	  if (!(*this)->is_free(Bi) &&
	      !this->mmap->is_marked((*this)->beta(Bi),
				     mmark_number))
            mto_treat.push((*this)->beta(Bi));
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
    int mmark_number;

    /// Initial dart
    Dart_handle minitial_dart;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Bi,Bi>: to iterate
   * on the darts of the orbit <Bi,Bj>: Bi<Bj<=dimension and (Bi,Bj)!=(0,1)
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_orbit.
   */
  template <typename Map_,bool Const,int Bi,int Bj>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj>: 
    public CMap_extend_iterator<Map_,CMap_dart_iterator_basic_of_orbit_generic
				<Map_,Const,Bi>, Bj>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj> Self;
    typedef CMap_extend_iterator<Map_,CMap_dart_iterator_basic_of_orbit_generic
				 <Map_,Const,Bi>, Bj> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
					     int amark):
      Base(amap, adart, amark)
    { CGAL_assertion( Bi<Bj && Bj!=1 && Bj<=Map::dimension ); }
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,0,3>: to iterate onto the
   * darts of the orbit <beta0, beta3> (i.e. orbit facet in 3D).
   * Specialized here since we do not need queue nor mark.
   */
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,3>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,3> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mit(amap, adart),
      mexist_beta3(false),
      mprev_beta3(false),
      mfirst_border(true)
    { if (adart!=NULL) mexist_beta3=!adart->is_free(3); }

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
					     int /*amark*/):
      Base(amap, adart),
      mit(amap, adart),
      mexist_beta3(false),
      mprev_beta3(false),
      mfirst_border(true)
    { if (adart!=NULL) mexist_beta3=!adart->is_free(3); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if (mexist_beta3 && !mprev_beta3)
	{
	  mprev_beta3 = true;
	  mfirst_border = ! mfirst_border;
	  this->set_current_dart((*this)->beta(3));
	  this->mprev_op = OP_BETAJ;
	}
      else
	{
	  mprev_beta3 = false;
	  ++mit;
	  this->mprev_op = mit.prev_operation();
	  if ( !mit.cont() ) 
	    this->set_current_dart(NULL);
	  else
	    {	  
	      if ( !mfirst_border ) 
		this->set_current_dart(mit->beta(3));
	      else
		this->set_current_dart(*mit);
	    }
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mit.rewind();
      mprev_beta3   = false;
      mfirst_border = true;
    }

  private:
    /// Iterator on beta0
    CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> mit;

    /// Boolean: true iff there are two half facets.
    bool mexist_beta3;

    /// Boolean: true iff the last ++ used beta3.
    bool mprev_beta3;

    /// Boolean: true iff the current dart is on the first border.
    bool mfirst_border;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,1,3>: to iterate onto the
   * darts of the orbit <beta1, beta3> (i.e. orbit facet in 3D).
   * Specialized here since we do not need queue nor mark.
   */
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,3>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,3> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mit(amap, adart),
      mexist_beta3(false),
      mprev_beta3(false),
      mfirst_border(true)
    { if (adart!=NULL) mexist_beta3=!adart->is_free(3); }

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
					     int /*amark*/):
      Base(amap, adart),
      mit(amap, adart),
      mexist_beta3(false),
      mprev_beta3(false),
      mfirst_border(true)
    { if (adart!=NULL) mexist_beta3=!adart->is_free(3); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if (mexist_beta3 && !mprev_beta3)
	{
	  mprev_beta3 = true;
	  mfirst_border = ! mfirst_border;
	  this->set_current_dart((*this)->beta(3));
	  this->mprev_op = OP_BETAJ;
	}
      else
	{
	  mprev_beta3 = false;
	  ++mit;
	  this->mprev_op = mit.prev_operation();
	  if ( !mit.cont() ) 
	    this->set_current_dart(NULL);
	  else
	    {	  
	      if ( !mfirst_border ) 
		this->set_current_dart(mit->beta(3));
	      else
		this->set_current_dart(mit);
	    }
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mit.rewind();
      mprev_beta3   = false;
      mfirst_border = true;
    }

  private:
    /// Iterator on beta1
    CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1> mit;

    /// Boolean: true iff there are two half facets.
    bool mexist_beta3;

    /// Boolean: true iff the last ++ used beta3.
    bool mprev_beta3;

    /// Boolean: true iff the current dart is on the first border.
    bool mfirst_border;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,2,3>: to iterate onto the
   * darts of the orbit <beta2, beta3> (i.e. orbit edge in 3D).
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_beta2(true)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int /*amark*/):
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_beta2(true)
    {}

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
      mnext_try_beta2   = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mfirst_dir)
	{
	  if (mnext_try_beta2)
	    {
	      if ((*this)->is_free(2))
		{
		  mfirst_dir = false;
		  if (this->mfirst_dart->is_free(3))
		    {
		      this->mprev_op = OP_END;
		      this->set_current_dart(NULL);
		    }
		  else
		    {
		      this->set_current_dart(this->mfirst_dart->beta(3));
		      this->mprev_op = OP_JUMP;
		    }
		}
	      else
		{
		  this->set_current_dart((*this)->beta(2));
		  mnext_try_beta2 = false;
		  this->mprev_op = OP_BETAI;
		}
	    }
	  else
	    {
	      if ((*this)->is_free(3))
		{
		  mfirst_dir = false;
		  if (this->mfirst_dart->is_free(3))
		    {
		      this->mprev_op = OP_END;
		      this->set_current_dart(NULL);
		    }
		  else
		    {
		      this->set_current_dart(this->mfirst_dart->beta(3));
		      mnext_try_beta2 = true;
		      this->mprev_op = OP_JUMP;
		    }
		}
	      else
		{
		  this->set_current_dart((*this)->beta(3));
		  if ((*this)==this->mfirst_dart)
		    {
		      this->mprev_op = OP_END;
		      this->set_current_dart(NULL);
		    }
		  else
		    {
		      mnext_try_beta2 = true;
		      this->mprev_op = OP_BETAJ;
		    }
		}
	    }
	}
      else
	{
	  if (mnext_try_beta2)
	    {
	      if ((*this)->is_free(2))
		{
		  this->mprev_op = OP_END;
		  this->set_current_dart(NULL);
		}
	      else
		{
		  this->set_current_dart((*this)->beta(2));
		  mnext_try_beta2 = false;
		  this->mprev_op = OP_BETAI;
		}
	    }
	  else
	    {
	      if ((*this)->is_free(3))
		{
		  this->mprev_op = OP_END;
		  this->set_current_dart(NULL);
		}
	      else
		{
		  this->set_current_dart((*this)->beta(3));
		  mnext_try_beta2 = true;
		  this->mprev_op = OP_BETAJ;
		}
	    }
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  private:
    /// Boolean: true iff we turn in the first direction (i.e. using beta2).
    bool mfirst_dir;

    /// Boolean: true iff the next ++ must use beta2.
    bool mnext_try_beta2;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Map,Bi,Bj,Bk>: to iterate onto 
   * the darts of the orbit <Bi,Bj,Bk>, Bi<Bj<Bk<=dimension 
   * Basic classes do not guaranty correct marks (i.e. do not unmark
   * darts in the destructor, possible problem with the rewind).
   * not for <B0,B3,Bk>, <B1,B3,Bk>, <B2,B3,Bk> which are specific cases.
   */
  template <typename Map_,int Bi,int Bj,int Bk,bool Const>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj,Bk>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj,Bk> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
					     int amark):
      Base(amap, adart, amark)
    {
      CGAL_assertion( Bi<Bj && Bj<Bk && Bj!=1 && Bk<=Map::dimension );

      if (adart!=NULL)
	{
	  if (!adart->is_free(Bk) &&
	      !this->mmap->is_marked(adart->beta(Bk), this->mmark_number))
	    this->mto_treat.push(adart->beta(Bk));
	}
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != -1);
      Base::rewind();
      if (!(*this)->is_free(Bk) &&
	  !this->mmap->is_marked((*this)->beta(Bk),
				 this->mmark_number))
	this->mto_treat.push((*this)->beta(Bk));
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      Base::operator++();

      if ( this->cont() )
	{
	  if (!(*this)->is_free(Bk) &&
	      !this->mmap->is_marked((*this)->beta(Bk),
				     this->mmark_number))
            this->mto_treat.push((*this)->beta(Bk));
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  template <typename Map_,bool Const,int Bi,int Bk>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,3,Bk>:
    public CMap_extend_iterator<Map_,CMap_dart_iterator_basic_of_orbit_generic
				<Map_,Const,Bi,3>, Bk>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,3,Bk> Self;
    typedef CMap_extend_iterator<Map_,CMap_dart_iterator_basic_of_orbit_generic
				 <Map_,Const,Bi,3>, Bk> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
					     int amark):
      Base(amap, adart, amark)
    { CGAL_assertion( Bi<3 && 3<Bk && Bk<=Map::dimension ); }
  };
  //****************************************************************************
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template<typename Map,int...Beta>
  class CMap_dart_iterator_basic_of_orbit: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map,false,Beta...>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit<Map,Beta...> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map,false,Beta...> Base;

    typedef typename Map::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart):
      Base(amap,adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart,int amark):
      Base(amap,adart,amark)
    {}
  };
  #else
  //****************************************************************************
  template<typename Map,int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1, 
	   int B6=-1,int B7=-1,int B8=-1,int B9=-1>
  class CMap_dart_iterator_basic_of_orbit: 
    public Get_CMap_dart_iterator_basic_of_orbit<Map,false,B1,B2,B3,B4,
						B5,B6,B7,B8,B9>::type
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit<Map,B1,B2,B3,B4,B5,B6,B7,B8,B9> 
    Self;
    typedef typename Get_CMap_dart_iterator_basic_of_orbit<Map,false,B1,B2,B3,B4,
							  B5,B6,B7,B8,B9>::type
    Base;

    typedef typename Map::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart):
      Base(amap,adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart,int amark):
      Base(amap,adart,amark)
    {}
  };
  #endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  //****************************************************************************
  // i-Cell iterator in combinatorial map of dimension d, i>1
  // i<=Map::dimension+1 (for i==Map::dimension+1, iterate on the connected
  // component)
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_basic_of_cell: public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,i,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      CGAL_assertion( i>=2 && i<=Map::dimension+1 );
      if (adart!=NULL)
	this->mmap->mark(adart, mmark_number); 
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark(*this, mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());
      Dart_handle nd = NULL;

      for ( unsigned int k=0; k<=d; ++k )
	{
	  if ( k!=i && this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // 0-Cell iterator in combinatorial map of dimension d
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,0,d,Const>:
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,0,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart),
      mmark_number(amark)
    { 
      if (adart!=NULL) this->mmap->mark(adart, mmark_number); 
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( unsigned int k=2; k<=d; ++k )
	{
	  if ( this->is_unmarked2((*this), 0, k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(0)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETA0I;
		}
	      else
		{
		  mto_treat.push((*this)->beta(0)->beta(k));
		  this->mmap->mark((*this)->beta(0)->beta(k), mmark_number);
		}
	      
	    }
	  if ( this->is_unmarked2((*this), k, 1, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k)->beta(1); 
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI1;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k)->beta(1));
		  this->mmap->mark((*this)->beta(k)->beta(1), mmark_number);
		}
	    }
	  for ( unsigned int l=k+1; l<=d; ++l )
	    {
	      if ( this->is_unmarked2((*this), k, l, mmark_number) )
		{
		  if (nd == NULL)
		    {
		      nd = (*this)->beta(k)->beta(l); 
		      CGAL_assertion(nd!=Map::null_dart_handle);
		      this->mmap->mark(nd, mmark_number);
		      this->mprev_op = OP_BETAIJ;
		    }
		  else
		    {
		      mto_treat.push((*this)->beta(k)->beta(l));
		      this->mmap->mark((*this)->beta(k)->beta(l), mmark_number);
		    }
		}
	      if ( this->is_unmarked2((*this), l, k, mmark_number) )
		{
		  if (nd == NULL)
		    {
		      nd = (*this)->beta(l)->beta(k); 
		      CGAL_assertion(nd!=Map::null_dart_handle);
		      this->mmap->mark(nd, mmark_number);
		      this->mprev_op = OP_BETAJI;
		    }
		  else
		    {
		      mto_treat.push((*this)->beta(l)->beta(k));
		      this->mmap->mark((*this)->beta(l)->beta(k), mmark_number);
		    }
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      CGAL_assertion(nd!=Map::null_dart_handle);
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // 1-Cell iterator in combinatorial map of dimension d
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,1,d,Const>:
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,1,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart),
      mmark_number(amark)
    { if (adart!=NULL) this->mmap->mark(adart, mmark_number); }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( unsigned int k=2; k<=d; ++k )
	{
	  if ( this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
   
  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // Specialization for edge in 2D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,1,2,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,1,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int /*amark*/):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for facet in 2D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,2,2,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,2,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int /*amark*/):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for cc in 2D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,3,2,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,3,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  // Specialization for edge in 3D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,1,3,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,1,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int /*amark*/): Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for facet in 3D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,2,3,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,3>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,2,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,3> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int /*amark*/): Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for volume in 3D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,3,3,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,3,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  // Specialization for cc in 3D 
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,4,3,Const>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2,3>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,4,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
				     Dart_handle adart,
				     int amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_cell<Map,0,2>: to iterate onto the
   * darts of the orbit vertex in 2D.
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_cell<Map_,0,2,Const>: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,0,2,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap, 
				     Dart_handle adart):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap, 
				     Dart_handle adart,
				     int /*amark*/):
      Base(amap, adart),
      mfirst_dir(true)
    {}

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mfirst_dir)
	{
	  this->set_current_dart((*this)->beta(0)->beta(2));
	  if ((*this)==Map::null_dart_handle)
	    {
	      mfirst_dir = false;
	      this->set_current_dart(this->mfirst_dart->beta(2)->beta(1));
	      if ((*this)==Map::null_dart_handle)
		{
		  this->mprev_op = OP_END;
		  this->set_current_dart(NULL);
		}
	      else
		{
		  this->mprev_op = OP_BETAI1;
		}
	    }
	  else
	    {
	      if ((*this)==this->mfirst_dart)
		{
		  this->mprev_op = OP_END;
		  this->set_current_dart(NULL);
		}
	      else
		this->mprev_op = OP_BETA0I;
	    }
	}
      else
	{
	  this->set_current_dart((*this)->beta(2)->beta(1));
	  if ((*this) == Map::null_dart_handle)
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	  else
            this->mprev_op = OP_BETA21;
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Boolean: true iff we turn in the first direction (i.e. using beta02).
    bool mfirst_dir;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_all: to iterate onto all the
   * darts of the map.
   */
  template <typename Map_,bool Const=false>
  class CMap_dart_iterator_basic_of_all: public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_all Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_all(Map& amap):
      Base(amap, amap.darts().begin())
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_all(Map& amap, int /*amark*/):
      Base(amap, amap.darts().begin())
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      Base::operator++();
      if ( (*this) != this->mmap->darts().end())
	{ this->mprev_op = OP_POP; }
      else
	{
	  this->set_current_dart(NULL);
	  this->mprev_op = OP_END;
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //**************************************************************************
  /* Generic nD version. Here we are sure that all the bases classes use mark
   * and queue.
   */
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename Map_,bool Const,int Bi,int... Beta>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Beta...>: 
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Beta...>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Beta...> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Beta...> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, 
					     Dart_handle adart, 
					     int amark):
      Base(amap, adart, amark)
    { 
      CGAL_assertion( Bi>=0 && Bi<=Map::dimension );

      if (adart!=NULL)
	{
	  if (!adart->is_free(Bi) &&
	      !this->mmap->is_marked(adart->beta(Bi), this->mmark_number))
	    this->mto_treat.push(adart->beta(Bi));
	}
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != -1);
      Base::rewind();
      if (!(*this)->is_free(Bi) &&
	  !this->mmap->is_marked((*this)->beta(Bi), 
				 this->mmark_number))
	this->mto_treat.push((*this)->beta(Bi));
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      Base::operator++();

      if (this->cont())
	{
	  if (!(*this)->is_free(Bi) &&
	      !this->mmap->is_marked((*this)->beta(Bi), 
				     this->mmark_number))
	    this->mto_treat.push((*this)->beta(Bi));
	}
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  //****************************************************************************
  //*************************ITERATORS*NON*BASIC*********************************
  //****************************************************************************
  //* Class CMap_non_basic_iterator allows to transform a basic_iterator onto
  //* a non basic one, depending if the basic iterator uses mark or not.
  template <typename Map_,typename Basic_iterator, 
	    typename Use_mark=typename Basic_iterator::Use_mark>
  class CMap_non_basic_iterator;
  //****************************************************************************
  template <typename Map_,typename Basic_iterator>
  class CMap_non_basic_iterator<Map_,Basic_iterator,Tag_true>: 
    public Basic_iterator
  {
  public:
    typedef CMap_non_basic_iterator<Map_,Basic_iterator,Tag_true> Self;

    typedef typename Basic_iterator::Map Map;
    typedef typename Basic_iterator::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_non_basic_iterator(Map& amap, Dart_handle adart1):
      Basic_iterator(amap, adart1, amap.get_new_mark())
    { CGAL_assertion( Basic_iterator::is_basic_iterator() ); }

    /// Destructor.
    ~CMap_non_basic_iterator()
    {
      if (this->mmark_number != -1)
	{
	  unmark_treated_darts();
	  CGAL_assertion( this->mmap->is_whole_map_unmarked
			  (this->mmark_number) );
	  this->mmap->free_mark(this->mmark_number);
	}
    }

    /// Copy constructor.
    CMap_non_basic_iterator(const Self& aiterator):
      Basic_iterator(aiterator)
    { this->mmark_number = -1; }

    /// Assignment operator.
    Self& operator=(const Self& aiterator)
    {
      if (this != &aiterator)
	{
	  Basic_iterator::operator=(aiterator);
	  this->mmark_number = -1;
	}
      return *this;
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(this->mmark_number != -1);
      unmark_treated_darts();
      Basic_iterator::rewind();
    }

    using Basic_iterator::operator++;

    /// Postfix ++ operator.
    void operator++(int)
    { operator ++(); }

    /// Return true iff this iterator is basic
    static bool is_basic_iterator()
    { return false; }

  protected:
    /// Unmark all the marked darts during the iterator.
    void unmark_treated_darts()
    {
      CGAL_assertion(this->mmark_number != -1);
      if (this->mmap->is_whole_map_unmarked(this->mmark_number)) return;

      this->mmap->negate_mark(this->mmark_number);

      if (this->mmap->is_whole_map_unmarked(this->mmark_number)) return;

      Basic_iterator::rewind();
      while (this->mmap->number_of_unmarked_darts(this->mmark_number) > 0)
	this->operator++();
      this->mmap->negate_mark(this->mmark_number);
      CGAL_assertion(this->mmap->is_whole_map_unmarked(this->mmark_number));
    }
  };
  //****************************************************************************
  template <typename Map_,typename Basic_iterator>
  class CMap_non_basic_iterator<Map_,Basic_iterator,Tag_false>: 
    public Basic_iterator
  {
  public:
    typedef CMap_non_basic_iterator<Map_,Basic_iterator,Tag_false> Self;
    typedef typename Basic_iterator::Map Map;
    typedef typename Basic_iterator::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_non_basic_iterator(Map& amap, Dart_handle adart):
      Basic_iterator(amap, adart,-1)
    {}
    /// Return true iff this iterator is basic
    static bool is_basic_iterator()
    { return false; }
  };
  //****************************************************************************
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template<typename Map_,bool Const,int...Beta>
  class CMap_dart_iterator_of_orbit_generic: 
    public  CMap_non_basic_iterator<Map_, 
				    CMap_dart_iterator_basic_of_orbit_generic
				    <Map_,Const,Beta...> >
  {
  public:
    typedef CMap_dart_iterator_of_orbit_generic<Map_,Const,Beta...> Self;
    typedef CMap_non_basic_iterator<Map_,
				    CMap_dart_iterator_basic_of_orbit_generic
				    <Map_,Const,Beta...> > Base;

    typedef typename Base::Map Map;
    typedef typename Base::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_dart_iterator_of_orbit_generic(Map& amap, Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  template<typename Map_,unsigned int...Beta>
  class CMap_dart_iterator_of_orbit: 
    public CMap_dart_iterator_of_orbit_generic<Map_,false,Beta...>
  {
  public:
    typedef CMap_dart_iterator_of_orbit<Map_,Beta...> Self;
    typedef CMap_dart_iterator_of_orbit_generic<Map_,false,Beta...> Base;

    typedef typename Base::Dart_handle Dart_handle;
    
    /// Main constructor.
    CMap_dart_iterator_of_orbit(Map_& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  #else
  //****************************************************************************
  template<typename Map_,bool Const,int B1=-1,int B2=-1,int B3=-1,int B4=-1,
	   int B5=-1,int B6=-1,int B7=-1,int B8=-1,int B9=-1>
  class CMap_dart_iterator_of_orbit_generic:
    public CMap_non_basic_iterator<Map_,
				   typename 
				   Get_CMap_dart_iterator_basic_of_orbit
				   <Map_,Const,B1,B2,B3,B4,B5,
				    B6,B7,B8,B9>::type>
  {
  public:
    typedef CMap_dart_iterator_of_orbit_generic<Map_,Const,B1,B2,B3,B4,B5,
					       B6,B7,B8,B9> Self;
    typedef CMap_non_basic_iterator<Map_,
                                    typename 
				    Get_CMap_dart_iterator_basic_of_orbit
				    <Map_,Const,B1,B2,B3,B4,B5,
				     B6,B7,B8,B9>::type> Base;

    typedef typename Base::Map Map;
    typedef typename Base::Dart_handle Dart_handle;

    /// Main constructor.
    CMap_dart_iterator_of_orbit_generic(Map& amap, Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  template<typename Map,int B1=-1,int B2=-1,int B3=-1,int B4=-1,
	   int B5=-1,int B6=-1,int B7=-1,int B8=-1,int B9=-1>
  class CMap_dart_iterator_of_orbit: 
    public CMap_dart_iterator_of_orbit_generic<Map,false,
					      B1,B2,B3,B4,B5,B6,B7,B8,B9>
  {
  public:
    typedef CMap_dart_iterator_of_orbit<Map,B1,B2,B3,B4,B5,B6,B7,B8,B9> Self;
    typedef CMap_dart_iterator_of_orbit_generic<Map,false,
					      B1,B2,B3,B4,B5,B6,B7,B8,B9> Base;

    typedef typename Base::Dart_handle Dart_handle;
    
    /// Main constructor.
    CMap_dart_iterator_of_orbit(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  #endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_of_cell: 
    public CMap_non_basic_iterator<Map_,CMap_dart_iterator_basic_of_cell
				   <Map_,i,d,Const> >
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,i,d,Const> Self;
    typedef CMap_non_basic_iterator<Map_,
				    CMap_dart_iterator_basic_of_cell
				    <Map_,i,d,Const> > Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    
    /// Main constructor.
    CMap_dart_iterator_of_cell(Map& amap, Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  //********************ITERATOR*INVOLUTION*************************************
  //****************************************************************************
  // i-involution iterator in combinatorial map of dimension d, 
  // 2<i<=Map::dimension. Iterate by using all beta between 0 and d, 
  // except beta(i-1), betai and beta(i+1)
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_basic_of_involution;

  template<typename Map_,int i,int d,bool Const>
  class CMap_dart_iterator_basic_of_involution: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,i,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      CGAL_assertion( d>=3 && d<=Map::dimension );
      CGAL_assertion( i>=3 && i<=Map::dimension );
      if (adart!=NULL) this->mmap->mark(adart, mmark_number); 
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( unsigned int k=0; k<2; ++k )
	{
	  if ( this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      for ( unsigned int k=2; k<=d; ++k )
	{
	  if ( k!=i-1 && k!=i && k!=i+1 && 
	       this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // i-involution iterator in combinatorial map of dimension d, 
  // 2<i<=Map::dimension. Iterate by using all beta between 0 and d, 
  // except beta(i-1), betai and beta(i+1), by inversing order between
  // beta0 and beta1
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_basic_of_involution_inv: 
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,i,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      CGAL_assertion( i>=3 && i<=Map::dimension );
      if (adart!=NULL) this->mmap->mark(adart, mmark_number); 
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( int k=1; k>=0; --k )
	{
	  if ( this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}
      for ( unsigned int k=2; k<=d; ++k )
	{
	  if ( k!=i-1 && k!=i && k!=i+1 && 
	       this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension d.
  // Iterate by using all beta between 3 and d.
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,1,d,Const>:  
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,1,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int amark):
      Base(amap, adart),
      mmark_number(amark)
    { if (adart!=NULL) this->mmap->mark(adart, mmark_number); }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( unsigned int k=3; k<=d; ++k )
	{
	  if ( this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension d.
  // Iterate by using all beta between 3 and d.
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,1,d,Const>:  
    public CMap_dart_iterator_basic_of_involution<Map_,1,d,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,1,d,Const> Self;
    typedef CMap_dart_iterator_basic_of_involution<Map_,1,d,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int amark):
      Base(amap, adart,amark)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension d.
  // Iterate by using all beta between 4 and d.
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,2,d,Const>:  
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,2,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int amark):
      Base(amap, adart),
      mmark_number(amark)
    { this->mmap->mark(adart, mmark_number); }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = NULL;

      for ( unsigned int k=4; k<=d; ++k )
	{
	  if ( this->is_unmarked((*this), k, mmark_number) )
	    {
	      if (nd == NULL)
		{
		  nd = (*this)->beta(k);
		  CGAL_assertion(nd!=Map::null_dart_handle);
		  this->mmap->mark(nd, mmark_number);
		  this->mprev_op = OP_BETAI;
		}
	      else
		{
		  mto_treat.push((*this)->beta(k));
		  this->mmap->mark((*this)->beta(k), mmark_number);
		}
	    }
	}

      if (nd == NULL)
	{
	  if (!mto_treat.empty())
	    {
	      nd = mto_treat.front();
	      mto_treat.pop();
	      this->mprev_op = OP_POP;
	    }
	  else
	    {
	      this->mprev_op = OP_END;
	      this->set_current_dart(NULL);
	    }
	}
      
      this->set_current_dart(nd);
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  protected:
    /// Queue of darts to process.
    std::queue<Dart_handle> mto_treat;

    /// Index of the used mark.
    int mmark_number;
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension d.
  // Iterate by using all beta between 4 and d.
  template<typename Map_,int d,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,2,d,Const>:  
    public CMap_dart_iterator_basic_of_involution<Map_,2,d,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,2,d,Const> Self;
    typedef CMap_dart_iterator_basic_of_involution<Map_,2,d,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int amark):
      Base(amap, adart,amark)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 2.
  // Empty iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,1,2,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,1,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 2.
  // self iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,1,2,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,1,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 2.
  // self iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,2,2,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,2,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 2.
  // self iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,2,2,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,2,2,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 3.
  // Beta3 iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,1,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,1,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    
    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart):
      Base(amap, adart)
    {}    
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 3.
  // Beta3 iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,1,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,1,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 3.
  // Self iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,2,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,2,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int /* amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 3.
  // Self iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,2,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,2,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 3.
  // Beta1 iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution<Map_,3,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,3,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart,
					   int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
					   Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 3.
  // Beta0 iterator.
  template<typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_involution_inv<Map_,3,3,Const>:  
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,3,3,Const> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart,
					       int /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
					       Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_of_involution:
    public CMap_non_basic_iterator<Map_,
				   CMap_dart_iterator_basic_of_involution
				   <Map_,i,d,Const> >
  {
  public:
    typedef CMap_dart_iterator_of_involution<Map_,i,d,Const> Self;
    typedef CMap_non_basic_iterator<Map_,
				    CMap_dart_iterator_basic_of_involution
				    <Map_,i,d,Const> >  Base;

    /// Main constructor.
    CMap_dart_iterator_of_involution(typename Base::Map& amap, 
				     typename Base::Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class CMap_dart_iterator_of_involution_inv: 
    public CMap_non_basic_iterator<Map_,
				   CMap_dart_iterator_basic_of_involution_inv
				   <Map_,i,d,Const> >
  {
  public:
    typedef CMap_dart_iterator_of_involution_inv<Map_,i,d,Const> Self;
    typedef CMap_non_basic_iterator<Map_,
				    CMap_dart_iterator_basic_of_involution_inv
				    <Map_,i,d,Const> >  Base;

    /// Main constructor.
    CMap_dart_iterator_of_involution_inv(typename Base::Map& amap,
					 typename Base::Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_DART_ITERATORS_HH
//******************************************************************************
