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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_DART_ITERATORS_HH
#define CGAL_DART_ITERATORS_HH 1

#include <CGAL/Combinatorial_map_iterators_base.h>

namespace CGAL {

  /** @file Dart_iterators.h
   * Definition of dart iterators. There are 9 iterators:
   * - CMap_dart_iterator_basic_of_orbit<Map,Beta...>
   * - CMap_dart_iterator_basic_of_cell<Map,i,d>
   * - CMap_dart_iterator_basic_of_all
   * - CMap_dart_iterator_basic_of_involution<Map,i,d>
   * - CMap_dart_iterator_basic_of_involution_inv<Map,i,d>
   * - CMap_dart_iterator_of_orbit<Map,Beta...>
   * - CMap_dart_iterator_of_cell<Map,i,d>
   * - CMap_dart_iterator_of_involution<Map,i,d>
   * - CMap_dart_iterator_of_involution_inv<Map,i,d>
   * but many specializations to optimize specific cases.
   *
   */
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
                                                      B1,B2,B3,B4,B5,
                                                      B6,B7> type;
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
      this->set_current_dart(this->mmap->null_handle);
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

      if (mfirst_dir && this->mmap->is_free(*this, 0))
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
        CGAL_assertion(!this->mmap->is_free(*this, 0));
        this->set_current_dart(this->mmap->beta(*this, 0));

        if ((*this)==this->mfirst_dart)
        {
          this->set_current_dart(this->mmap->null_handle);
          this->mprev_op = OP_END;
        }
      }
      else
      {
        if (this->mmap->is_free(*this, 1))
        {
          this->set_current_dart(this->mmap->null_handle);
          this->mprev_op = OP_END;
        }
        else
        {
          this->set_current_dart(this->mmap->beta(*this, 1));
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

      if (mfirst_dir && this->mmap->is_free(*this, 1))
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
        CGAL_assertion(!this->mmap->is_free(*this, 1));
        this->set_current_dart(this->mmap->beta(*this, 1));

        if ((*this)==this->mfirst_dart)
        {
          this->set_current_dart(this->mmap->null_handle);
          this->mprev_op = OP_END;
        }
      }
      else
      {
        if (this->mmap->is_free(*this, 0))
        {
          this->set_current_dart(this->mmap->null_handle);
          this->mprev_op = OP_END;
        }
        else
        {
          this->set_current_dart(this->mmap->beta(*this, 0));
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
    { CGAL_static_assertion( Bi>=2 && Bi<=Map::dimension ); }

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart,
                                              int /*amark*/):
      Base(amap, adart)
    { CGAL_static_assertion( Bi>=2 && Bi<=Map::dimension ); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if ((*this)!=this->mfirst_dart || this->mmap->is_free(*this, Bi))
      {
        this->set_current_dart(this->mmap->null_handle);
        this->mprev_op = OP_END;
      }
      else
      {
        this->set_current_dart(this->mmap->beta(*this, Bi));
        this->mprev_op = OP_BETAI;
      }
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * This general case if for Bi>1 and delta>1.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const,int Bi,unsigned int delta>
  class CMap_dart_iterator_basic_of_two_beta :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,Bi,delta> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

    CGAL_static_assertion( Bi>1 && delta>1 && Bi+delta<=Map::dimension );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mcurdart(0)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart,
                                         int /*amark*/):
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
        if (!this->mmap->is_free(*this, Bi))
        {
          this->set_current_dart(this->mmap->beta(*this, Bi));
          this->mprev_op = OP_BETAI;
          mcurdart=1;
        }
        else
        {
          if (!this->mmap->is_free(*this, Bi+delta))
          {
            this->set_current_dart(this->mmap->beta(*this, Bi+delta));
            this->mprev_op = OP_BETAJ;
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
        if (!this->mmap->is_free(*this, Bi+delta))
        {
          this->set_current_dart(this->mmap->beta(*this, Bi+delta));
          this->mprev_op = OP_BETAJ;
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
        CGAL_assertion(!this->mmap->is_free(*this, Bi));
        this->set_current_dart(this->mmap->beta(*this, Bi));
        this->mprev_op = OP_BETAI;
        mcurdart=1;
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
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * Special case for Bi==0 and delta==2.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_two_beta<Map_,Const,0,2> :
    public CMap_extend_iterator
  <Map_, CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>, 2>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,0,2> Self;
    typedef CMap_extend_iterator
    <Map_, CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0>, 2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

    CGAL_static_assertion( 2<=Map::dimension );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart,
                                         int amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * Special case for Bi==1 and delta==1.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const>
  class CMap_dart_iterator_basic_of_two_beta<Map_,Const,1,1> :
    public CMap_extend_iterator
  <Map_, CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1>, 2>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,1,1> Self;
    typedef CMap_extend_iterator
    <Map_, CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,1>, 2> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_true Use_mark;

    CGAL_static_assertion( 2<=Map::dimension );
    
  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart, 
                                         int amark):
      Base(amap, adart, amark)
    {}
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * Special case for Bi==0 and delta>2.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const, unsigned int delta>
  class CMap_dart_iterator_basic_of_two_beta<Map_,Const,0,delta> :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,0,delta> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

    CGAL_static_assertion( delta>1 && delta<=Map::dimension );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mit(amap, adart),
      mexist_betaj(false),
      mprev_betaj(false),
      mfirst_border(true)
    { if (adart!=this->mmap->null_handle)
        mexist_betaj=!this->mmap->is_free(adart, delta); }
    
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart, 
                                         int /*amark*/):
      Base(amap, adart),
      mit(amap, adart),
      mexist_betaj(false),
      mprev_betaj(false),
      mfirst_border(true)
    { if (adart!=this->mmap->null_handle)
        mexist_betaj=!this->mmap->is_free(adart, delta); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if (mexist_betaj && !mprev_betaj)
      {
        mprev_betaj = true;
        mfirst_border = ! mfirst_border;
        this->set_current_dart(this->mmap->beta(*this, delta));
        this->mprev_op = OP_BETAJ;
      }
      else
      {
        mprev_betaj = false;
        ++mit;
        this->mprev_op = mit.prev_operation();
        if ( !mit.cont() ) 
          this->set_current_dart(this->mmap->null_handle);
        else
        {          
          if ( !mfirst_border ) 
            this->set_current_dart(mit->beta(delta));
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
      mprev_betaj   = false;
      mfirst_border = true;
    }

  private:
    /// Iterator on beta0
    CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,0> mit;

    /// Boolean: true iff there are two half facets.
    bool mexist_betaj;

    /// Boolean: true iff the last ++ used betaj.
    bool mprev_betaj;

    /// Boolean: true iff the current dart is on the first border.
    bool mfirst_border;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * Special case for Bi==1 and delta>1.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const, unsigned int delta>
  class CMap_dart_iterator_basic_of_two_beta<Map_,Const,1,delta> :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,1,delta> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

    CGAL_static_assertion( delta>1 && delta+1<=Map::dimension );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart):
      Base(amap, adart),
      mit(amap, adart),
      mexist_betaj(false),
      mprev_betaj(false),
      mfirst_border(true)
    { if (adart!=this->mmap->null_handle)
        mexist_betaj=!this->mmap->is_free(adart, 1+delta); }
    
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart, 
                                         int /*amark*/):
      Base(amap, adart),
      mit(amap, adart),
      mexist_betaj(false),
      mprev_betaj(false),
      mfirst_border(true)
    { if (adart!=this->mmap->null_handle)
        mexist_betaj=!this->mmap->is_free(adart, 1+delta); }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      if (mexist_betaj && !mprev_betaj)
      {
        mprev_betaj = true;
        mfirst_border = ! mfirst_border;
        this->set_current_dart(this->mmap->beta(*this, 1+delta));
        this->mprev_op = OP_BETAJ;
      }
      else
      {
        mprev_betaj = false;
        ++mit;
        this->mprev_op = mit.prev_operation();
        if ( !mit.cont() ) 
          this->set_current_dart(this->mmap->null_handle);
        else
        {          
          if ( !mfirst_border ) 
            this->set_current_dart(this->mmap->beta(mit, 1+delta));
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
      mprev_betaj   = false;
      mfirst_border = true;
    }

  private:
    /// Iterator on beta1
    CMap_dart_iterator_basic_of_orbit_generic<Map_,Const, 1> mit;

    /// Boolean: true iff there are two half facets.
    bool mexist_betaj;

    /// Boolean: true iff the last ++ used betaj.
    bool mprev_betaj;

    /// Boolean: true iff the current dart is on the first border.
    bool mfirst_border;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_two_beta<Bi,delta>: to iterate
   * on the darts of the orbit <Bi,Bi+delta>: Bi<Bi+delta<=dimension.
   * Special case for Bi>1 and delta==1.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_two_beta.
   */
  template <typename Map_,bool Const, int Bi>
  class CMap_dart_iterator_basic_of_two_beta<Map_,Const,Bi,1> :
    public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,Bi,1> Self;
    typedef CMap_dart_iterator<Map_,Const> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef Tag_false Use_mark;

    CGAL_static_assertion( Bi>1 && Bi+1<=Map::dimension );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart):    
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_betai(true)
    {}

    /// Main constructor.
    CMap_dart_iterator_basic_of_two_beta(Map& amap, Dart_handle adart,
                                         int /*amark*/):
      Base(amap, adart),
      mfirst_dir(true),
      mnext_try_betai(true)
    {}
    
    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      Base::rewind();
      mfirst_dir = true;
      mnext_try_betai   = true;
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());
      
      if (mfirst_dir)
      {
        if (mnext_try_betai)
        {
          if (this->mmap->is_free(*this, Bi))
          {
            mfirst_dir = false;
            if (this->mmap->is_free(this->mfirst_dart, Bi+1))
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              this->set_current_dart(this->mmap->beta(this->mfirst_dart, Bi+1));
              this->mprev_op = OP_JUMP;
            }
          }
          else
          {
            this->set_current_dart(this->mmap->beta(*this, Bi));
            mnext_try_betai = false;
            this->mprev_op = OP_BETAI;
          }
        }
        else
        {
          if (this->mmap->is_free(*this, Bi+1))
          {
            mfirst_dir = false;
            if (this->mmap->is_free(this->mfirst_dart, Bi+1))
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              this->set_current_dart(this->mmap->beta(this->mfirst_dart, Bi+1));
              mnext_try_betai = true;
              this->mprev_op = OP_JUMP;
            }
          }
          else
          {
            this->set_current_dart(this->mmap->beta(*this, Bi+1));
            if ((*this)==this->mfirst_dart)
            {
              this->mprev_op = OP_END;
              this->set_current_dart(this->mmap->null_handle);
            }
            else
            {
              mnext_try_betai = true;
              this->mprev_op = OP_BETAJ;
            }
          }
        }
      }
      else
      {
        if (mnext_try_betai)
        {
          if (this->mmap->is_free(*this, Bi))
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
          }
          else
          {
            this->set_current_dart(this->mmap->beta(*this, Bi));
            mnext_try_betai = false;
            this->mprev_op = OP_BETAI;
          }
        }
        else
        {
          if (this->mmap->is_free(*this, Bi+1))
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
          }
          else
          {
            this->set_current_dart(this->mmap->beta(*this, Bi+1));
            mnext_try_betai = true;
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
    /// Boolean: true iff we turn in the first direction (i.e. using betai).
    bool mfirst_dir;

    /// Boolean: true iff the next ++ must use betai.
    bool mnext_try_betai;
  };
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_orbit<Bi,Bj>: to iterate
   * on the darts of the orbit <Bi,Bj>: Bi<Bj<=dimension.
   * Basic classes do not guaranty correct marks (i.e. do not unmark darts in
   * the destructor, possible problem with the rewind). If you are not sure,
   * use CMap_dart_iterator_basic_of_orbit.
   */
  template <typename Map_,bool Const,int Bi,int Bj>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj>: 
    public CMap_dart_iterator_basic_of_two_beta<Map_,Const,Bi,Bj-Bi>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj> Self;
    typedef CMap_dart_iterator_basic_of_two_beta<Map_,Const,Bi,Bj-Bi> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    typedef typename Base::Use_mark Use_mark;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart) :
      Base(amap, adart)
    {}    

    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
                                              int amark):
      Base(amap, adart, amark)
    {}    
  };
  //****************************************************************************
  /* Generic nD version. 
   */
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename Map_,bool Const,int Bi,int Bj, int Bk, int... Beta>
  class CMap_dart_iterator_basic_of_orbit_generic<Map_,Const,Bi,Bj,Bk,Beta...>: 
    public CMap_extend_iterator<Map_,
                                CMap_dart_iterator_basic_of_orbit_generic
                                <Map_,Const,Bi,Bj,Beta...>,
                                Bk>
  {
  public:
    typedef CMap_dart_iterator_basic_of_orbit_generic
    <Map_,Const,Bi,Bj,Bk,Beta...> Self;
    typedef CMap_extend_iterator<Map_,
                                 CMap_dart_iterator_basic_of_orbit_generic
                                 <Map_,Const,Bi,Bj,Beta...>,
                                 Bk> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
                                              int amark):
      Base(amap, adart, amark)
    {}
  };
#else //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename Map_,bool Const,int B1,int B2,int B3,int B4,int B5,int B6,
            int B7,int B8,int B9>
  class CMap_dart_iterator_basic_of_orbit_generic: 
    public CMap_extend_iterator
  <Map_,typename Get_CMap_dart_iterator_basic_of_orbit
   <Map_,Const,B1,B2,B4,B5,B6,B7,B8,B9>::type, B3>
  {
  public:
    typedef typename Get_CMap_dart_iterator_basic_of_orbit
    <Map_,Const,B1,B2,B3,B4,B5,B6,B7,B8,B9>::type Self;

    typedef CMap_extend_iterator
    <Map_,typename Get_CMap_dart_iterator_basic_of_orbit
     <Map_,Const,B1,B2,B4,B5,B6,B7,B8,B9>::type, B3> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_orbit_generic(Map& amap, Dart_handle adart, 
                                              int amark):
      Base(amap, adart, amark)
    {}
  };
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  //****************************************************************************
  // TODO? we can optimize the iterators<Bi,Bj,Bk> when
  // 1<Bi and Bi+2<=Bj and Bj+2<=Bk but there is no real interest...
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
    typedef typename Get_CMap_dart_iterator_basic_of_orbit
    <Map,false,B1,B2,B3,B4,B5,B6,B7,B8,B9>::type Base;

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

    /// Constructor with a dart in parameter (for end iterator).
    CMap_dart_iterator_basic_of_all(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    CMap_dart_iterator_basic_of_all(Map& amap, Dart_handle adart,
                                    int /*amark*/):
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
  //***************************CELL*ITERATORS***********************************
  //****************************************************************************
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

    CGAL_static_assertion( i>1 && i<=Map::dimension+1 );

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_cell(Map& amap,
                                     Dart_handle adart,
                                     int amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark_null_dart(mmark_number);
        this->mmap->mark(adart, mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark(*this, mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());
      Dart_handle nd = this->mmap->null_handle;
      
      for ( unsigned int k=0; k<i; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
        }
      }
      for ( unsigned int k=i+1; k<=d; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
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
    int mmark_number;
  };
  //****************************************************************************
  // i-Cell iterator in combinatorial map of dimension d, i==1.
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
    {
      if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;
      
      for ( unsigned int k=2; k<=d; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
        }
      }

      if (nd == this->mmap->null_handle)
      {
        if (!mto_treat.empty())
        {
          nd = mto_treat.front();
          CGAL_assertion(nd!=this->mmap->null_dart_handle);
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
    { if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }      
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( unsigned int k=2; k<=d; ++k )
      {
        if ( this->is_unmarked2((*this), 0, k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, 0, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETA0I;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, 0, k));
          }
          this->mmap->mark(this->mmap->beta(*this, 0, k), mmark_number);
        }
        if ( this->is_unmarked2((*this), k, 1, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k, 1);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI1;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k, 1));
          }
          this->mmap->mark(this->mmap->beta(*this, k, 1), mmark_number);
        }
        for ( unsigned int l=k+1; l<=d; ++l )
        {
          if ( this->is_unmarked2((*this), k, l, mmark_number) )
          {
            if (nd == this->mmap->null_handle)
            {
              nd = this->mmap->beta(*this, k, l);
              CGAL_assertion(nd!=this->mmap->null_dart_handle);
              this->mprev_op = OP_BETAIJ;
            }
            else
            {
              mto_treat.push(this->mmap->beta(*this, k, l));
            }
            this->mmap->mark(this->mmap->beta(*this, k, l), mmark_number);
          }
          if ( this->is_unmarked2((*this), l, k, mmark_number) )
          {
            if (nd == this->mmap->null_handle)
            {
              nd = this->mmap->beta(*this, l, k);
              CGAL_assertion(nd!=this->mmap->null_dart_handle);
              this->mprev_op = OP_BETAJI;
            }
            else
            {
              mto_treat.push(this->mmap->beta(*this, l, k));
            }
            this->mmap->mark(this->mmap->beta(*this, l, k), mmark_number);
          }
        }
      }

      if (nd == this->mmap->null_handle)
      {
        if (!mto_treat.empty())
        {
          nd = mto_treat.front();
          CGAL_assertion(nd!=this->mmap->null_dart_handle);
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
        this->set_current_dart(this->mmap->beta(*this, 0, 2));
        if ((*this)==this->mmap->null_dart_handle)
        {
          mfirst_dir = false;
          this->set_current_dart(this->mmap->beta(this->mfirst_dart, 2, 1));
          if ((*this)==this->mmap->null_dart_handle)
          {
            this->mprev_op = OP_END;
            this->set_current_dart(this->mmap->null_handle);
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
            this->set_current_dart(this->mmap->null_handle);
          }
          else
            this->mprev_op = OP_BETA0I;
        }
      }
      else
      {
        this->set_current_dart(this->mmap->beta(*this, 2, 1));
        if ((*this) == this->mmap->null_dart_handle)
        {
          this->mprev_op = OP_END;
          this->set_current_dart(this->mmap->null_handle);
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
  //*************************ITERATORS*NON*BASIC********************************
  //****************************************************************************
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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;
    
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
      if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( int k=0; k<2; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
        }
      }

      for ( int k=2; k<=d; ++k )
      {
        if ( k!=i-1 && k!=i && k!=i+1 && 
             this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution_inv(Map& amap,
                                               Dart_handle adart,
                                               int amark):
      Base(amap, adart),
      mmark_number(amark)
    {
      CGAL_assertion( i>=3 && i<=Map::dimension );
      if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( int k=1; k>=0; --k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
        }
      }
      for ( int k=2; k<=d; ++k )
      {
        if ( k!=i-1 && k!=i && k!=i+1 && 
             this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           int amark):
      Base(amap, adart),
      mmark_number(amark)
    { if (adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark_null_dart(mmark_number);
      this->mmap->mark((*this), mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( unsigned int k=3; k<=d; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

  public:
    /// Main constructor.
    CMap_dart_iterator_basic_of_involution(Map& amap,
                                           Dart_handle adart,
                                           int amark):
      Base(amap, adart),
      mmark_number(amark)
    { if ( adart!=this->mmap->null_handle)
      {
        this->mmap->mark(adart, mmark_number);
        this->mmap->mark_null_dart(mmark_number);
      }
    }

    /// Rewind of the iterator to its beginning.
    void rewind()
    {
      CGAL_assertion(mmark_number != -1);
      Base::rewind();
      mto_treat = std::queue<Dart_handle>();
      this->mmap->mark((*this), mmark_number);
      this->mmap->mark_null_dart(mmark_number);
    }

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(mmark_number != -1);
      CGAL_assertion(this->cont());

      Dart_handle nd = this->mmap->null_handle;

      for ( unsigned int k=4; k<=d; ++k )
      {
        if ( this->is_unmarked((*this), k, mmark_number) )
        {
          if (nd == this->mmap->null_handle)
          {
            nd = this->mmap->beta(*this, k);
            CGAL_assertion(nd!=this->mmap->null_dart_handle);
            this->mprev_op = OP_BETAI;
          }
          else
          {
            mto_treat.push(this->mmap->beta(*this, k));
          }
          this->mmap->mark(this->mmap->beta(*this, k), mmark_number);
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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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

    /// True iff this iterator is basic
    typedef Tag_true Basic_iterator;

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
