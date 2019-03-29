// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_ITERATORS_FOR_FACE_GRAPH_WRAPPER_H
#define CGAL_ITERATORS_FOR_FACE_GRAPH_WRAPPER_H 1

#include <CGAL/Combinatorial_map_iterators_base.h>

namespace CGAL
{
  
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_all: to iterate onto all the
   * darts of the face graph.
   */
  template <typename Map_,bool Const=false>
  class FGW_dart_iterator_basic_of_all: public CMap_dart_iterator<Map_,Const>
  {
  public:
    typedef FGW_dart_iterator_basic_of_all Self;
    typedef CMap_dart_iterator<Map_,Const> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Map Map;
    typedef typename Map::size_type size_type;

    typedef Tag_false Use_mark;

  public:
    /// Main constructor.
    FGW_dart_iterator_basic_of_all(Map& amap):
      Base(amap, *(CGAL::halfedges(amap.get_fg()).begin())),
      m_it(CGAL::halfedges(amap.get_fg()).begin())
    {
      if (m_it!=CGAL::halfedges(amap.get_fg()).end() &&
          CGAL::is_border(*m_it, amap.get_fg()))
      { operator++(0); }
    }
    /// Main constructor.
    FGW_dart_iterator_basic_of_all(Map& amap, size_type /*amark*/):
      Base(amap, CGAL::halfedges(amap.get_fg()).begin())
    {}

    /// Constructor with a dart in parameter (for end iterator).
    FGW_dart_iterator_basic_of_all(Map& amap, Dart_handle adart):
      Base(amap, adart)
    {}
    /// Constructor with a dart in parameter (for end iterator).
    FGW_dart_iterator_basic_of_all(Map& amap, Dart_handle adart,
                                    size_type /*amark*/):
      Base(amap, adart)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(this->cont());

      do
      {
        ++m_it;
        this->set_current_dart(*m_it);
      }
      while(m_it!=CGAL::halfedges(this->mmap->get_fg()).end() &&
            CGAL::is_border(*m_it, this->mmap->get_fg()));
      
      if (m_it!=CGAL::halfedges(this->mmap->get_fg()).end())
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

  protected:
    typename boost::graph_traits<typename Map_::HEG>::halfedge_iterator m_it;
  };
}

// ////////////////////////////////////////////////////////////////////////////////
// template<typename HEG, unsigned int i>
// class Cell_iterator
// {};
// ////////////////////////////////////////////////////////////////////////////////
// /// Vertex iterator
// template<typename HEG>
// class Cell_iterator<HEG, 0> // Vertex
// {
// public:
//   typedef Cell_iterator Self;
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   Cell_iterator(Map& amap, Dart_const_handle adart) : mmap(amap), mfirst_dart(adart),
//                                                       mcurrent_dart(adart),
//                                                       m_second_way(false)
//   {}

//   /// == operator.
//   bool operator==(const Self& other) const
//   { return (mfirst_dart==other.mfirst_dart &&
//             mcurrent_dart==other.mcurrent_dart); }

//   /// != operator.
//   bool operator!=(const Self& other) const
//   { return !operator==(other); }

//   /// Prefix ++ operator.
//   Self& operator++()
//   {
//     return *this;
//   }

//   /// Postfix ++ operator.
//   Self operator++(int)
//   { Self res=*this; operator ++(); return res; }

// protected:
//   /// The map containing the darts to iterate on.
//   const HEG& mmap;
//   /// The initial dart of the iterator.
//   Dart_const_handle mfirst_dart;
//   /// The current dart of the iterator.
//   Dart_const_handle mcurrent_dart;
//   /// True if we already found a border dart, and thus turn in the second way
//   bool m_second_way;
// };
// template<typename HEG>
// class Cell_iterator<HEG, 1> // Edge
// {
// public:
//   typedef Cell_iterator Self;
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   Cell_iterator(Map& amap, Dart_const_handle adart) : mmap(amap), mfirst_dart(adart)
//   {}

//   /// == operator.
//   bool operator==(const Self& other) const
//   { return (mfirst_dart==other.mfirst_dart &&
//             mcurrent_dart==other.mcurrent_dart); }

//   /// != operator.
//   bool operator!=(const Self& other) const
//   { return !operator==(other); }

//   /// Prefix ++ operator.
//   Self& operator++()
//   {
//     return *this;
//   }

//   /// Postfix ++ operator.
//   Self operator++(int)
//   { Self res=*this; operator ++(); return res; }

// protected:
//   /// The map containing the darts to iterate on.
//   const HEG& mmap;
//   /// The initial dart of the iterator.
//   Dart_const_handle mfirst_dart;
//   /// The current dart of the iterator.
//   Dart_const_handle mcurrent_dart;
// };
// template<typename HEG>
// class Cell_iterator<HEG, 2> // Face
// {
// public:
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

//   Cell_iterator(Map& amap, Dart_const_handle adart) : mmap(amap), mfirst_dart(adart)
//   {}
  
//   /// == operator.
//   bool operator==(const Self& other) const
//   { return (mfirst_dart==other.mfirst_dart &&
//             mcurrent_dart==other.mcurrent_dart); }

//   /// != operator.
//   bool operator!=(const Self& other) const
//   { return !operator==(other); }

//   /// Prefix ++ operator.
//   Self& operator++()
//   {
//     return *this;
//   }

//   /// Postfix ++ operator.
//   Self operator++(int)
//   { Self res=*this; operator ++(); return res; }

// protected:
//   /// The map containing the darts to iterate on.
//   const HEG& mmap;
//   /// The initial dart of the iterator.
//   Dart_const_handle mfirst_dart;
//   /// The current dart of the iterator.
//   Dart_const_handle mcurrent_dart;
// };
// template<typename HEG>
// class Cell_iterator<HEG, 3> // Connected component
// {
// public:
//   typedef Cell_iterator Self;
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  
//   Cell_iterator(Map& amap, Dart_const_handle adart) : mmap(amap), mfirst_dart(adart)
//   {}
  
// protected:
//   /// The map containing the darts to iterate on.
//   const HEG& mmap;
//   /// The initial dart of the iterator.
//   Dart_const_handle mfirst_dart;
//   /// The current dart of the iterator.
//   Dart_const_handle mcurrent_dart;
// };
// template<typename HEG>
// class Cell_iterator<HEG, 4> // All darts
// {
// public:
//   typedef Cell_iterator Self;
//   typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;
  
//   Cell_iterator(Map& amap, Dart_const_handle adart) : mmap(amap), mfirst_dart(adart)
//   {}
  
//   /// == operator.
//   bool operator==(const Self& other) const
//   { return (mfirst_dart==other.mfirst_dart &&
//             mcurrent_dart==other.mcurrent_dart); }

//   /// != operator.
//   bool operator!=(const Self& other) const
//   { return !operator==(other); }

//   /// Prefix ++ operator.
//   Self& operator++()
//   {
//     return *this;
//   }

//   /// Postfix ++ operator.
//   Self operator++(int)
//   { Self res=*this; operator ++(); return res; }

// protected:
//   /// The map containing the darts to iterate on.
//   const HEG& mmap;
//   /// The initial dart of the iterator.
//   Dart_const_handle mfirst_dart;
//   /// The current dart of the iterator.
//   Dart_const_handle mcurrent_dart;
// };
// ////////////////////////////////////////////////////////////////////////////////

#endif // CGAL_ITERATORS_FOR_FACE_GRAPH_WRAPPER_H //
// EOF //
