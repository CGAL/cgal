// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GMAP_DART_ITERATORS_HH
#define CGAL_GMAP_DART_ITERATORS_HH 1

#include <CGAL/Generalized_map_iterators_base.h>

namespace CGAL {

  /** @file GMap_dart_iterators.h
   * Definition of gmap dart iterators.
   * There are 7 iterators:
   *  - GMap_dart_iterator_basic_of_orbit<Map,Alpha...>
   *  - GMap_dart_iterator_basic_of_cell<Map,i,d>
   *  - GMap_dart_iterator_basic_of_all<Map>
   *  - GMap_dart_iterator_basic_of_involution<Map,i,d>
   *  - GMap_dart_iterator_of_orbit<Map,Alpha...>
   *  - GMap_dart_iterator_of_cell<Map,i,d>
   *  - GMap_dart_iterator_of_involution<Map,i,d>
   * but many specializations to optimize specific cases.
   *
   */
  //****************************************************************************
  //**********************BASIC ITERATORS***************************************
  //****************************************************************************
  /* Class GMap_dart_iterator_basic_of_orbit<Map, Alpha...>: to iterate
   * on the darts of the orbit <Alpha...>
   */
  template<typename Map,bool Const,int... Alpha>
  class GMap_dart_iterator_basic_of_orbit_generic;
  //****************************************************************************
  // Case when Alpha... is empty: iterator of self
  template <typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_orbit_generic<Map_,Const>:
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;      ///< True iff this iterator uses mark

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              size_type /*amark*/):
      Base(amap, adart)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      this->set_current_dart(this->mmap->null_handle);
      this->mprev_op = OP_END;
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  /* Class GMap_dart_iterator_basic_of_orbit<Ai>: to iterate
   * on the darts of the orbit <Ai> (0<=Ai<=dimension)
   */
  template <typename Map_,bool Const,int Ai>
  class GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai>:
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;      ///< True iff this iterator uses mark

    CGAL_static_assertion( Ai>=0 && Ai<=Map::dimension );

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              size_type /*amark*/):
      Base(amap, adart)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if ((*this)!=this->mfirst_dart || this->mmap->is_free(*this, Ai))
      {
        this->set_current_dart(this->mmap->null_handle);
        this->mprev_op = OP_END;
      }
      else
      {
        this->set_current_dart(this->mmap->template alpha<Ai>(*this));
        this->mprev_op = OP_ALPHAI;
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_alpha<Ai,delta>: to iterate
   * on the darts of the orbit <Ai,Ai+delta>: Ai<Ai+delta<=dimension.
   * This general case if for delta>1 (ie at most 4 darts).
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_alpha.
   */
  template <typename Map_,bool Const,int Ai,unsigned int delta>
  class GMap_dart_iterator_basic_of_two_alpha :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_two_alpha<Map_,Const,Ai,delta> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;      ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    CGAL_static_assertion( (0<=Ai && Ai+delta<=Map::dimension && delta>1) );

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_two_alpha(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mcurdart(0)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_two_alpha(Map& amap, Dart_handle adart,
                                          size_type /*amark*/):
      Base(amap, adart),
      mcurdart(0)
    {}

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mcurdart=0;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mcurdart==0)
      {
        if (!this->mmap->template is_free<Ai>(*this))
        {
          this->set_current_dart(this->mmap->template alpha<Ai>(*this));
          this->mprev_op = OP_ALPHAI;
          mcurdart=1;
        }
        else
        {
          if (!this->mmap->template is_free<Ai+delta>(*this))
          {
            this->set_current_dart(this->mmap->template alpha<Ai+delta>(*this));
            this->mprev_op = OP_ALPHAJ;
            mcurdart=3;
          }
          else
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
          }
        }
      }
      else if (mcurdart==1)
      {
        if (!this->mmap->template is_free<Ai+delta>(*this))
        {
          this->set_current_dart(this->mmap->template alpha<Ai+delta>(*this));
          this->mprev_op = OP_ALPHAJ;
          mcurdart=2;
        }
        else
        {
          this->mprev_op = OP_END;
          this->set_current_dart(this->mmap->null_handle);
        }
      }
      else if (mcurdart==2)
      {
        CGAL_assertion(!this->mmap->template is_free<Ai>(*this));
        this->set_current_dart(this->mmap->template alpha<Ai>(*this));
        this->mprev_op = OP_ALPHAI;
        mcurdart=3;
      }
      else
      {
        CGAL_assertion (mcurdart==3);
        this->mprev_op = OP_END;
        this->set_current_dart(this->mmap->null_handle);
      }

      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  private:
    /// mcurdart: number of the current dart (0,1,2 or 3).
    char mcurdart;
  };
  //****************************************************************************
  template <typename Map_,bool Const,int Ai, unsigned int delta>
  class GMap_dart_iterator_basic_of_two_alpha;
  /* Class CMap_dart_iterator_basic_of_two_alpha<Ai,1>: to iterate
   * on the darts of the orbit <Ai,Ai+1>: Ai<Ai+1<=dimension.
   * specialisation because here Aio(Ai+1) is not an involution.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_alpha.
   */
  template <typename Map_,bool Const,int Ai>
  class GMap_dart_iterator_basic_of_two_alpha<Map_,Const,Ai,1> :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_two_alpha<Map_,Const,Ai,1> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;      ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    CGAL_static_assertion(0<=Ai && Ai+1<=Map_::dimension);

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_two_alpha(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_first_alpha(true)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_two_alpha(Map& amap, Dart_handle adart,
                                          size_type /*amark*/):
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_first_alpha(true)
    {}

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
      mnext_try_first_alpha = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      if (mfirst_dir)
      {
        if (mnext_try_first_alpha)
        {
          if (this->mmap->template is_free<Ai>(*this))
          {
            mfirst_dir = false;
            if (this->mmap->template is_free<Ai+1>(this->mfirst_dart))
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              this->set_current_dart
                (this->mmap->template alpha<Ai+1>(this->mfirst_dart));
              this->mprev_op = OP_JUMP;
            }
          }
          else
          {
            this->set_current_dart(this->mmap->template alpha<Ai>(*this));
            mnext_try_first_alpha = false;
            this->mprev_op = OP_ALPHAI;
          }
        }
        else
        {
          if (this->mmap->template is_free<Ai+1>(*this))
          {
            mfirst_dir = false;
            if (this->mmap->template is_free<Ai+1>(this->mfirst_dart))
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              this->set_current_dart
                (this->mmap->template alpha<Ai+1>(this->mfirst_dart));
              mnext_try_first_alpha = true;
              this->mprev_op = OP_JUMP;
            }
          }
          else
          {
            this->set_current_dart(this->mmap->template alpha<Ai+1>(*this));
            if ((*this)==this->mfirst_dart)
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              mnext_try_first_alpha = true;
              this->mprev_op = OP_ALPHAJ;
            }
          }
        }
      }
      else
      {
        if (mnext_try_first_alpha)
        {
          if (this->mmap->is_free(*this, Ai))
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
          }
          else
          {
            this->set_current_dart(this->mmap->template alpha<Ai>(*this));
            mnext_try_first_alpha = false;
            this->mprev_op = OP_ALPHAI;
          }
        }
        else
        {
          if (this->mmap->template is_free<Ai+1>(*this))
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
          }
          else
          {
            this->set_current_dart(this->mmap->template alpha<Ai+1>(*this));
            mnext_try_first_alpha = true;
            this->mprev_op = OP_ALPHAJ;
          }
        }
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

  private:
    /// Boolean: true iff we turn in the first direction (i.e. using Ai).
    bool mfirst_dir;

    /// Boolean: true iff the next ++ must use Ai.
    bool mnext_try_first_alpha;
  };
  //****************************************************************************
  /* Class GMap_dart_iterator_basic_of_orbit<Ai,Aj>: to iterate
   * on the darts of the orbit <Ai,Aj>: Ai<dimension.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use GMap_dart_iterator_basic_of_orbit.
   */
  template <typename Map_,bool Const, int Ai, int Aj>
  class GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj>:
    public GMap_dart_iterator_basic_of_two_alpha<Map_,Const,Ai,Aj-Ai>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj> Self;
    typedef GMap_dart_iterator_basic_of_two_alpha<Map_,Const,Ai,Aj-Ai> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart) :
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              size_type amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /// TODO The template specialization with 3 alpha

  //
  /* Class CMap_dart_iterator_basic_of_three_alpha<Ai,delta1,delta2>: to iterate
   * on the darts of the orbit <Ai,Ai+delta1,Ai+delta2>:
   * Ai<Ai+delta1<Ai+delta2<=dimension.
   * This general case if for delta1>1 and delta2>1 (ie at most 8 darts).
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_alpha.
   */
  /*
    template <typename Map_,bool Const,int Ai,unsigned int delta1,
            unsigned int delta2>
  class GMap_dart_iterator_basic_of_three_alpha :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_three_alpha<Map_,Const,Ai,
                                                    delta1,delta2> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;      ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    CGAL_static_assertion( (0<=Ai && delta1<delta2 &&
                            Ai+delta2<=Map::dimension &&
                            delta1>1) );

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_three_alpha(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_three_alpha(Map& amap, Dart_handle adart,
                                            size_type):
      Base(amap, adart)
    {}
  }; */
  //****************************************************************************
  /* Class GMap_dart_iterator_basic_of_orbit<Map,Ai,Aj,Ak>: to iterate onto
   * the darts of the orbit <Ai,Aj,Ak>, Ai<Aj<Ak<=dimension
   * Basic classes do not guaranty correct marks (i.e. do not unmark
   * darts in the destructor, possible problem with the rewind).
   * not for <B0,A3,Ak>, <A1,A3,Ak>, <A2,A3,Ak> which are specific cases.
   */

  /*  template <typename Map_,bool Const,int Ai,int Aj,int Ak>
  class GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj,Ak>:
    public GMap_extend_iterator<Map_,GMap_dart_iterator_basic_of_orbit_generic
                                <Map_,Const,Ai,Aj>, Ak>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj,Ak> Self;
    typedef GMap_extend_iterator<Map_,GMap_dart_iterator_basic_of_orbit_generic
                                 <Map_,Const,Ai,Aj>, Ak> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map         Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    CGAL_static_assertion( Ai<Aj && Aj<Ak && Ak<=Map::dimension );

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              size_type amark):
      Base(amap, adart, amark)
    {}
    };
  */
  //**************************************************************************
  /// Generic nD version.
  template <typename Map_,bool Const,int Ai, int Aj, int... Alpha>
  class GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj,Alpha...>:
    public GMap_extend_iterator<Map_,
                                GMap_dart_iterator_basic_of_orbit_generic
                                <Map_,Const,Aj,Alpha...>,
                                Ai>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Ai,Aj,Alpha...> Self;
    typedef GMap_extend_iterator<Map_,
                                 GMap_dart_iterator_basic_of_orbit_generic
                                 <Map_,Const,Aj,Alpha...>,
                                 Ai> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              size_type amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /// Non const basic of orbit iterator
  template<typename Map,int Ai, int Aj, int...Alpha>
  class GMap_dart_iterator_basic_of_orbit:
    public GMap_dart_iterator_basic_of_orbit_generic<Map,false,Ai,Aj,Alpha...>
  {
  public:
    typedef GMap_dart_iterator_basic_of_orbit<Map,Ai,Aj,Alpha...> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map,false,Ai,Aj,Alpha...> Base;

    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Base::Use_mark Use_mark; ///< True iff this iterator uses mark
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart):
      Base(amap,adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_orbit(Map& amap,Dart_handle adart,size_type amark):
      Base(amap,adart,amark)
    {}
  };
  //****************************************************************************
  // Generic i-Cell iterator in generalized map of dimension d,
  // i<=Map::dimension+1 (for i==Map::dimension+1, iterate on the connected
  // component)
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class GMap_dart_iterator_basic_of_cell: public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,i,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Use_mark; ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

    CGAL_static_assertion( i>=0 && i<=Map::dimension+1 );

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      if (adart!=this->mmap->null_handle)
        this->mmap->mark(adart, mmark_number);
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark(*this, mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      CGAL_assertion(this->cont());
      Dart_handle nd = this->mmap->null_handle;

      for ( unsigned int k=0; k<=d; ++k )
      {
        if ( k!=i &&
             !this->mmap->is_free(*this, k) &&
             !this->mmap->is_marked(this->mmap->alpha(*this, k),
                                    this->mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            CGAL_assertion(!this->mmap->is_free(*this, k));
            nd = this->mmap->alpha(*this, k);
            this->mmap->mark(nd, mmark_number);
            this->mprev_op = OP_ALPHAI;
          }
          else
          {
            mto_treat.push(this->mmap->alpha(*this, k));
          }
         this->mmap->mark(this->mmap->alpha(*this, k), mmark_number);
        }
      }

      if (nd == this->mmap->null_handle)
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
    size_type mmark_number;
  };
  //****************************************************************************
  // Specialization for vertex in 2D
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_cell<Map_,0,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,0,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type /*amark*/): Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for edge in 2D
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_cell<Map_,1,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,2>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,1,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic
    typedef typename Map::size_type size_type;

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type /*amark*/): Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for facet in 2D
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_cell<Map_,2,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,2,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic
    typedef typename Map::size_type size_type;

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart):
      Base(amap, adart)
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type /*amark*/):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // Specialization for edge in 3D
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_cell<Map_,1,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,2,3>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,1,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,2,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic
    typedef typename Map::size_type size_type;

    /// Main constructor.
    /// @TODO specialization of iterator
    /// GMap_dart_iterator_basic_of_orbit_generic
    /// with 3 non consecutive args
    /*GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart):
      Base(amap, adart)
    {}*/

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type amark): Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  // Specialization for facet in 3D
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_cell<Map_,2,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1,3>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,2,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic
    typedef typename Map::size_type size_type;

    /// Main constructor.
    /// @TODO specialization of iterator
    /// GMap_dart_iterator_basic_of_orbit_generic
    /// with 3 non consecutive args
    /*GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart):
      Base(amap, adart)
    {}*/

    /// Main constructor.
    GMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     size_type amark): Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /* Class GMap_dart_iterator_basic_of_all: to iterate onto all the
   * darts of the map.
   */
  template <typename Map_,bool Const=false>
  class GMap_dart_iterator_basic_of_all: public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_all Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark; ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_all(Map& amap):
      Base(amap, amap.darts().begin())
    {}

    /// Main constructor.
    GMap_dart_iterator_basic_of_all(Map& amap, size_type /*amark*/):
      Base(amap, amap.darts().begin())
    {}

    /// Constructor with a dart in parameter (for end iterator).
    GMap_dart_iterator_basic_of_all(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    GMap_dart_iterator_basic_of_all(Map& amap, Dart_handle adart,
                                    size_type /*amark*/):
      Base(amap, adart)
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
        this->set_current_dart(this->mmap->null_handle);
        this->mprev_op = OP_END;
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  //*************************ITERATORS*NON*BASIC*********************************
  //****************************************************************************
  template<typename Map_,bool Const,int...Alpha>
  class GMap_dart_iterator_of_orbit_generic:
    public  CMap_non_basic_iterator<Map_,
                                    GMap_dart_iterator_basic_of_orbit_generic
                                    <Map_,Const,Alpha...> >
  {
  public:
    typedef GMap_dart_iterator_of_orbit_generic<Map_,Const,Alpha...> Self;
    typedef CMap_non_basic_iterator<Map_,
                                    GMap_dart_iterator_basic_of_orbit_generic
                                    <Map_,Const,Alpha...> > Base;

    typedef typename Base::Map Map;
    typedef typename Base::Dart_handle Dart_handle;
    typedef Tag_false Basic_iterator; ///< True iff this iterator is basic

    /// Main constructor.
    GMap_dart_iterator_of_orbit_generic(Map& amap, Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  template<typename Map_,unsigned int...Alpha>
  class GMap_dart_iterator_of_orbit:
    public GMap_dart_iterator_of_orbit_generic<Map_,false,Alpha...>
  {
  public:
    typedef GMap_dart_iterator_of_orbit<Map_,Alpha...> Self;
    typedef GMap_dart_iterator_of_orbit_generic<Map_,false,Alpha...> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef Tag_false Basic_iterator; ///< True iff this iterator is basic

    /// Main constructor.
    GMap_dart_iterator_of_orbit(Map_& amap, Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class GMap_dart_iterator_of_cell:
    public CMap_non_basic_iterator<Map_,GMap_dart_iterator_basic_of_cell
                                   <Map_,i,d,Const> >
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,i,d,Const> Self;
    typedef CMap_non_basic_iterator<Map_,
                                    GMap_dart_iterator_basic_of_cell
                                    <Map_,i,d,Const> > Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef Tag_false Basic_iterator; ///< True iff this iterator is basic

    /// Main constructor.
    GMap_dart_iterator_of_cell(Map& amap, Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
  //********************ITERATOR*INVOLUTION*************************************
  //****************************************************************************
  // i-involution iterator in combinatorial map of dimension d,
  // 0<i<=Map::dimension. Iterate by using all alpha between 0 and d,
  // except alpha(i-1), alphai and alpha(i+1)
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class GMap_dart_iterator_basic_of_involution;

  template<typename Map_,int i,int d,bool Const>
  class GMap_dart_iterator_basic_of_involution:
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,i,d,Const> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Use_mark; ///< True iff this iterator uses mark
    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      CGAL_assertion( d<=Map::dimension );
      CGAL_assertion( i<=Map::dimension );
      if (adart!=this->mmap->null_handle)
        this->mmap->mark(adart, mmark_number);
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != Map::INVALID_MARK);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( int k=0; k<i-1; ++k )
      {
        if ( !this->mmap->is_free(*this, k) &&
             !this->mmap->is_marked(this->mmap->alpha(*this, k),
                                    this->mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            CGAL_assertion(!this->mmap->is_free(*this, k));
            nd = this->mmap->alpha(*this, k);
            this->mmap->mark(nd, mmark_number);
            this->mprev_op = OP_ALPHAI;
          }
          else
          {
            mto_treat.push(this->mmap->alpha(*this, k));
          }
         this->mmap->mark(this->mmap->alpha(*this, k), mmark_number);
        }
      }

      for ( unsigned int k=i+2; k<=d; ++k )
      {
        if ( !this->mmap->is_free(*this, k) &&
             !this->mmap->is_marked(this->mmap->alpha(*this, k),
                                    this->mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            CGAL_assertion(!this->mmap->is_free(*this, k));
            nd = this->mmap->alpha(*this, k);
            this->mmap->mark(nd, mmark_number);
            this->mprev_op = OP_ALPHAI;
          }
          else
          {
            mto_treat.push(this->mmap->alpha(*this, k));
          }
         this->mmap->mark(this->mmap->alpha(*this, k), mmark_number);
        }
      }

      if (nd == this->mmap->null_handle)
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
    size_type mmark_number;
  };
  //****************************************************************************
  // 0-involution iterator in combinatorial map of dimension 2.
  // Alpha2 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,0,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,1,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 2.
  // Self iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,1,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,1,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 2.
  // Alpha0 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,2,2,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,2,2,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 0-involution iterator in combinatorial map of dimension 3.
  // Alpha2, Alpha3 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,0,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,0,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,2,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 1-involution iterator in combinatorial map of dimension 3.
  // Alpha3 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,1,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,1,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 2-involution iterator in combinatorial map of dimension 3.
  // Alpha0 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,2,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,2,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /* amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  // 3-involution iterator in combinatorial map of dimension 3.
  // Alpha0, Alpha1 iterator.
  template<typename Map_,bool Const>
  class GMap_dart_iterator_basic_of_involution<Map_,3,3,Const>:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,3,3,Const> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0,1> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_true Basic_iterator; ///< True iff this iterator is basic

  public:
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           size_type /*amark*/):
      Base(amap, adart)
    {}
    /// Main constructor.
    GMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart):
      Base(amap, adart)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension,bool Const=false>
  class GMap_dart_iterator_of_involution:
    public CMap_non_basic_iterator<Map_,
                                   GMap_dart_iterator_basic_of_involution
                                   <Map_,i,d,Const> >
  {
  public:
    typedef GMap_dart_iterator_of_involution<Map_,i,d,Const> Self;
    typedef CMap_non_basic_iterator<Map_,
                                    GMap_dart_iterator_basic_of_involution
                                    <Map_,i,d,Const> >  Base;

    /// Main constructor.
    GMap_dart_iterator_of_involution(typename Base::Map& amap,
                                     typename Base::Dart_handle adart1):
      Base(amap, adart1)
    {}
  };
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_GMAP_DART_ITERATORS_HH
//******************************************************************************
