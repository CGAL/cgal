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
#ifndef CGAL_GMAP_DART_CONST_ITERATORS_HH
#define CGAL_GMAP_DART_CONST_ITERATORS_HH 1

#include <CGAL/GMap_dart_iterators.h>

namespace CGAL {

  /** @file GMap_dart_const_iterators.h
   * Definition of gmap dart const iterators.
   * There are 7 iterators:
   *  - GMap_dart_const_iterator_basic_of_orbit<Map,Alpha...>
   *  - GMap_dart_const_iterator_basic_of_cell<Map,i,d>
   *  - GMap_dart_const_iterator_basic_of_involution<Map,i,d>
   *  - GMap_dart_const_iterator_basic_of_all<Map>
   *  - GMap_dart_const_iterator_of_orbit<Map,Alpha...>
   *  - GMap_dart_const_iterator_of_cell<Map,i,d>
   *  - GMap_dart_const_iterator_of_involution<Map,i,d>
   */
  //****************************************************************************
  template<typename Map_,unsigned int...Alpha>
  class GMap_dart_const_iterator_basic_of_orbit:
    public GMap_dart_iterator_basic_of_orbit_generic<Map_,true,Alpha...>
  {
  public:
    typedef GMap_dart_const_iterator_basic_of_orbit<Map_,Alpha...> Self;
    typedef GMap_dart_iterator_basic_of_orbit_generic<Map_,true,Alpha...> Base;

    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /// Main constructor.
    GMap_dart_const_iterator_basic_of_orbit(const Map_& amap,
                                            Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Main constructor.
    GMap_dart_const_iterator_basic_of_orbit(const Map_& amap,
                                            Dart_const_handle adart,
                                            size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_basic_of_orbit
    (const GMap_dart_const_iterator_basic_of_orbit<Map_,Alpha...>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart(),
           it.mmark_number)
    {}
  };
  //****************************************************************************
  template<typename Map_,unsigned int...Alpha>
  class GMap_dart_const_iterator_of_orbit:
    public GMap_dart_iterator_of_orbit_generic<Map_,true,Alpha...>
  {
  public:
    typedef GMap_dart_const_iterator_of_orbit<Map_,Alpha...> Self;
    typedef GMap_dart_iterator_of_orbit_generic<Map_,true,Alpha...> Base;

    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /// Main constructor.
    GMap_dart_const_iterator_of_orbit(const Map_& amap,
                                      Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_of_orbit
    (const GMap_dart_iterator_of_orbit<Map_,Alpha...>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_>
  class GMap_dart_const_iterator_basic_of_all:
    public GMap_dart_iterator_basic_of_all<Map_,true>
  {
  public:
    typedef GMap_dart_iterator_basic_of_all<Map_,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_all(const Map_& amap,
                                          Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_all(const Map_& amap,
                                          Dart_const_handle adart,
                                          size_type /*amark*/):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_basic_of_all
    (const GMap_dart_iterator_basic_of_all<Map_,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class GMap_dart_const_iterator_basic_of_cell:
    public GMap_dart_iterator_basic_of_cell<Map_,i,d,true>
  {
  public:
    typedef GMap_dart_iterator_basic_of_cell<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_cell(const Map_& amap,
                                           Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_cell(const Map_& amap,
                                           Dart_const_handle adart,
                                           size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_basic_of_cell
    (const GMap_dart_iterator_basic_of_cell<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_, int i, int d=Map_::dimension>
  class GMap_dart_const_iterator_of_cell:
    public GMap_dart_iterator_of_cell<Map_,i,d,true>
  {
  public:
    typedef GMap_dart_iterator_of_cell<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;

    /* Main constructor. */
    GMap_dart_const_iterator_of_cell(const Map_& amap,
                                     Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_of_cell
    (const GMap_dart_iterator_of_cell<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class GMap_dart_const_iterator_basic_of_involution:
    public GMap_dart_iterator_basic_of_involution<Map_,i,d,true>
  {
  public:
    typedef GMap_dart_iterator_basic_of_involution<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_involution(const Map_& amap,
                                                 Dart_const_handle adart):
      Base(amap,adart)
    {}
    /* Main constructor. */
    GMap_dart_const_iterator_basic_of_involution(const Map_& amap,
                                                 Dart_const_handle adart,
                                                 size_type amark):
      Base(amap,adart,amark)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_basic_of_involution
    (const GMap_dart_iterator_basic_of_involution<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart(), it.mmark_number)
    {}
  };
  //****************************************************************************
  template<typename Map_,int i,int d=Map_::dimension>
  class GMap_dart_const_iterator_of_involution:
    public GMap_dart_iterator_of_involution<Map_,i,d,true>
  {
  public:
    typedef GMap_dart_iterator_of_involution<Map_,i,d,true> Base;
    typedef typename Map_::Dart_const_handle Dart_const_handle;
    typedef typename Map_::size_type size_type;

    /* Main constructor. */
    GMap_dart_const_iterator_of_involution(const Map_& amap,
                                           Dart_const_handle adart):
      Base(amap,adart)
    {}
    /// Constructor from non const version.
    GMap_dart_const_iterator_of_involution
    (const GMap_dart_iterator_of_involution<Map_,i,d,false>& it):
      Base(*const_cast<const Map_*>(it.get_combinatorial_map()),
           it.get_first_dart())
    {}
  };
  //****************************************************************************
} // namespace CGAL
//******************************************************************************
#endif // CGAL_GMAP_DART_CONST_ITERATORS_HH
//******************************************************************************
