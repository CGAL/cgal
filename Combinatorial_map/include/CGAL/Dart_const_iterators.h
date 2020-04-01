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
#ifndef CGAL_DART_CONST_ITERATORS_HH
#define CGAL_DART_CONST_ITERATORS_HH 1

#include <CGAL/Dart_iterators.h>

namespace CGAL {

  /** @file Dart_const_iterators.h
   * Definition of dart const iterators. There are 9 iterators:
   *  - CMap_dart_const_iterator_basic_of_orbit<Map,Beta...>
   *  - CMap_dart_const_iterator_basic_of_cell<Map,i,d>
   *  - CMap_dart_const_iterator_basic_of_all<Map>
   *  - CMap_dart_const_iterator_basic_of_involution<Map,i,d>
   *  - CMap_dart_const_iterator_of_involution_inv<Map,i,d>
   *  - CMap_dart_const_iterator_of_orbit<Map,Beta...>
   *  - CMap_dart_const_iterator_of_cell<Map,i,d>
   *  - CMap_dart_const_iterator_of_involution<Map,i,d>
   *  - CMap_dart_const_iterator_basic_of_involution_inv<Map,i,d>
   */
  //****************************************************************************

  template<typename Map_,unsigned int...Beta>
  class CMap_dart_const_iterator_basic_of_orbit:
    public CMap_dart_iterator_basic_of_orbit_generic<Map_,true,Beta...>
  {
  public:
    typedef CMap_dart_const_iterator_basic_of_orbit<Map_,Beta...> Self;
    typedef CMap_dart_iterator_basic_of_orbit_generic<Map_,true,Beta...> Base;

    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /// Main constructor.
    CMap_dart_const_iterator_basic_of_orbit(const Map_& amap,
                                            Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Main constructor.
    CMap_dart_const_iterator_basic_of_orbit(const Map_& amap,
                                            Dart_const_handle adart,
                                            size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_basic_of_orbit
    (const CMap_dart_iterator_basic_of_orbit<Map_,Beta...>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart(),
           it.mmark_number)
    {}
  };
  //****************************************************************************
  template<typename Map_,unsigned int...Beta>
  class CMap_dart_const_iterator_of_orbit:
    public CMap_dart_iterator_of_orbit_generic<Map_,true,Beta...>
  {
  public:
    typedef CMap_dart_const_iterator_of_orbit<Map_,Beta...> Self;
    typedef CMap_dart_iterator_of_orbit_generic<Map_,true,Beta...> Base;

    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /// Main constructor.
    CMap_dart_const_iterator_of_orbit(const Map_& amap,
                                      Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_of_orbit
    (const CMap_dart_iterator_of_orbit<Map_,Beta...>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_>
  class CMap_dart_const_iterator_basic_of_all:
    public CMap_dart_iterator_basic_of_all<Map_,true>
  {
  public:
    typedef CMap_dart_iterator_basic_of_all<Map_,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_all(const Map_& amap,
                                          Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_all(const Map_& amap,
                                          Dart_const_handle adart,
                                          size_type /*amark*/):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_basic_of_all
    (const CMap_dart_iterator_basic_of_all<Map_,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class CMap_dart_const_iterator_basic_of_cell:
    public CMap_dart_iterator_basic_of_cell<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_basic_of_cell<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_cell(const Map_& amap,
                                           Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_cell(const Map_& amap,
                                           Dart_const_handle adart,
                                           size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_basic_of_cell
    (const CMap_dart_iterator_basic_of_cell<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_, int i, int d=Map_::dimension>
  class CMap_dart_const_iterator_of_cell:
    public CMap_dart_iterator_of_cell<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_of_cell<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /* Main constructor. */
    CMap_dart_const_iterator_of_cell(const Map_& amap,
                                     Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_of_cell
    (const CMap_dart_iterator_of_cell<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class CMap_dart_const_iterator_basic_of_involution:
    public CMap_dart_iterator_basic_of_involution<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_involution(const Map_& amap,
                                                 Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_involution(const Map_& amap,
                                                 Dart_const_handle adart,
                                                 size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_basic_of_involution
    (const CMap_dart_iterator_basic_of_involution<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart(), it.mmark_number)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class CMap_dart_const_iterator_of_involution:
    public CMap_dart_iterator_of_involution<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_of_involution<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /* Main constructor. */
    CMap_dart_const_iterator_of_involution(const Map_& amap,
                                           Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_of_involution
    (const CMap_dart_iterator_of_involution<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class CMap_dart_const_iterator_basic_of_involution_inv:
    public CMap_dart_iterator_basic_of_involution_inv<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_basic_of_involution_inv<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_involution_inv(const Map_& amap,
                                                     Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    CMap_dart_const_iterator_basic_of_involution_inv(const Map_& amap,
                                                     Dart_const_handle adart,
                                                     size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_basic_of_involution_inv
    (const CMap_dart_iterator_basic_of_involution_inv<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart(), it.mmark_number)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class CMap_dart_const_iterator_of_involution_inv:
    public CMap_dart_iterator_of_involution_inv<Map_,i,d,true>
  {
  public:
    typedef CMap_dart_iterator_of_involution_inv<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /* Main constructor. */
    CMap_dart_const_iterator_of_involution_inv(const Map_& amap,
                                               Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    CMap_dart_const_iterator_of_involution_inv
    (const CMap_dart_iterator_of_involution_inv<Map_,i,d>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_DART_CONST_ITERATORS_HH
//******************************************************************************
