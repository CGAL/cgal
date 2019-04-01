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

namespace CGAL
{
  
  //****************************************************************************
  /* Class CMap_dart_iterator_basic_of_all: to iterate onto all the
   * darts of the face graph.
   */
  template <typename Map_,bool Const=false>
  class FGW_dart_iterator_basic_of_all
  {
  public:
    typedef FGW_dart_iterator_basic_of_all Self;

    typedef Map_ Map;
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::size_type size_type;

  public:
    /// Main constructor.
    FGW_dart_iterator_basic_of_all(const Map& amap):
      mmap(amap),
      m_it(CGAL::halfedges(amap.get_fg()).begin())
    {
      if (m_it!=CGAL::halfedges(amap.get_fg()).end() &&
          CGAL::is_border(*m_it, amap.get_fg()))
      { operator++(0); }
    }

    /// Constructor with a dart in parameter (for end iterator).
    FGW_dart_iterator_basic_of_all(const Map& amap, Dart_handle /*adart*/):
      mmap(amap),
      m_it(CGAL::halfedges(amap.get_fg()).end())
    {}

    FGW_dart_iterator_basic_of_all(const FGW_dart_iterator_basic_of_all& other):
      mmap(other.mmap),
      m_it(other.m_it)      
    {}
    
    bool operator==(const Self& other) const
    { return &mmap==&(other.mmap) && m_it==other.m_it; }

    bool operator!=(const Self& other) const
    { return !(operator==(other)); }
    
    /// Prefix ++ operator.
    Self& operator++()
    {
      CGAL_assertion(m_it!=CGAL::halfedges(this->mmap.get_fg()).end());

      do
      {
        ++m_it;
      }
      while(m_it!=CGAL::halfedges(this->mmap.get_fg()).end() &&
            CGAL::is_border(*m_it, this->mmap.get_fg()));
      
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }

    Dart_handle operator*()
    {
      CGAL_assertion(m_it!=CGAL::halfedges(this->mmap.get_fg()).end());
      return *m_it;
    }
    
  protected:
    const Map& mmap;
    typename boost::graph_traits<typename Map_::HEG>::halfedge_iterator m_it;
  };

////////////////////////////////////////////////////////////////////////////////
  template <typename Map_>
  class FGW_basis_for_cell_iterator
  {
  public:
    typedef FGW_basis_for_cell_iterator Self;
    typedef Map_ Map;    
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::size_type size_type;
    
    /// Main constructor.
    FGW_basis_for_cell_iterator(const Map& amap, Dart_handle adart):
      mmap(amap),
      m_firstdart(adart),
      m_curdart(adart)
    {}

    /// Constructor with two darts in parameter (for end iterator).
    FGW_basis_for_cell_iterator(const Map& amap, Dart_handle adart,
                                Dart_handle /* d2 */):
      mmap(amap),
      m_firstdart(adart),
      m_curdart(Dart_handle())
    {}

    bool operator==(const Self& other) const
    { return &mmap==&(other.mmap) && m_firstdart==other.m_firstdart &&
        m_curdart==other.m_curdart; }

    bool operator!=(const Self& other) const
    { return !(this->operator==(other)); }
    
    Dart_handle operator*()
    {
      CGAL_assertion(m_curdart!=Dart_handle());
      return m_curdart;
    }
    
  protected:
    const Map& mmap;
    Dart_handle m_firstdart, m_curdart;
  };
////////////////////////////////////////////////////////////////////////////////
  template<typename HEG, unsigned int i>
  class FGW_cell_iterator
  {};
////////////////////////////////////////////////////////////////////////////////
/// Vertex iterator
  template<typename Map_>
  class FGW_cell_iterator<Map_, 0>: public FGW_basis_for_cell_iterator<Map_> // Vertex
  {
  public:
    typedef FGW_cell_iterator Self;
    typedef FGW_basis_for_cell_iterator<Map_> Base;
    typedef Map_ Map;    
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::size_type size_type;
    
    FGW_cell_iterator(const Map& amap, Dart_handle adart) : Base(amap, adart),
                                                            m_second_way(false)
    {}

    /// Constructor with two darts in parameter (for end iterator).
    FGW_cell_iterator(const Map& amap, Dart_handle adart,
                      Dart_handle d2): Base(amap, adart, d2)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      if (!m_second_way)
      {
        if (this->mmap.template is_free<2>(this->m_curdart))
        {
          m_second_way=true;
          this->m_curdart=this->mmap.template beta<0>(this->m_firstdart);
          if (this->mmap.template is_free<2>(this->m_curdart))
          { this->m_curdart=Dart_handle(); }
          else { this->m_curdart=this->mmap.template beta<2>(this->m_curdart); }
        }
        else
        {
          this->m_curdart=this->mmap.template beta<2, 1>(this->m_curdart);
          if (this->m_curdart==this->m_firstdart) 
          { this->m_curdart=Dart_handle(); }
        }
      }
      else
      {
        this->m_curdart=this->mmap.template beta<0>(this->m_curdart);
        if (this->mmap.template is_free<2>(this->m_curdart))
        { this->m_curdart=Dart_handle(); }
        else { this->m_curdart=this->mmap.template beta<2>(this->m_curdart); }
      }
      
      return *this;
    }

    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }
    
protected:
  /// True if we already found a border dart, and thus turn in the second way
  bool m_second_way;
};
template<typename Map_>
class FGW_cell_iterator<Map_, 1>: public FGW_basis_for_cell_iterator<Map_>  // Edge
{
  public:
    typedef FGW_cell_iterator Self;
    typedef FGW_basis_for_cell_iterator<Map_> Base;
    typedef Map_ Map;    
    typedef typename Map::Dart_handle Dart_handle;
    typedef typename Map::size_type size_type;
    
    FGW_cell_iterator(const Map& amap, Dart_handle adart) : Base(amap, adart)
    {}

    /// Constructor with two darts in parameter (for end iterator).
    FGW_cell_iterator(const Map& amap, Dart_handle adart,
                      Dart_handle d2): Base(amap, adart, d2)
    {}

    /// Prefix ++ operator.
    Self& operator++()
    {
      if (this->m_curdart==this->m_firstdart)
      {
        if (this->mmap.template is_free<2>(this->m_curdart))
        { this->m_curdart=Dart_handle(); }
        else { this->m_curdart=this->mmap.template beta<2>(this->m_curdart); }
      }
      else
      { this->m_curdart=Dart_handle(); }
      
      return *this;
    }
    
    /// Postfix ++ operator.
    Self operator++(int)
    { Self res=*this; operator ++(); return res; }    
};
template<typename Map_>
class FGW_cell_iterator<Map_, 2>: public FGW_basis_for_cell_iterator<Map_> // Face
{
public:
  typedef FGW_cell_iterator Self;
  typedef FGW_basis_for_cell_iterator<Map_> Base;
  typedef Map_ Map;    
  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::size_type size_type;
  
  FGW_cell_iterator(const Map& amap, Dart_handle adart) : Base(amap, adart)
  {}
  
  /// Constructor with two darts in parameter (for end iterator).
  FGW_cell_iterator(const Map& amap, Dart_handle adart,
                    Dart_handle d2): Base(amap, adart, d2)
  {}

  /// Prefix ++ operator.
  Self& operator++()
  {
    this->m_curdart=this->mmap.template beta<1>(this->m_curdart);
    if (this->m_curdart==this->m_firstdart)
    { this->m_curdart=Dart_handle(); }
    
    return *this;
  }
  
  /// Postfix ++ operator.
  Self operator++(int)
  { Self res=*this; operator ++(); return res; }
};
////////////////////////////////////////////////////////////////////////////////

} // namespace CGAL

#endif // CGAL_ITERATORS_FOR_FACE_GRAPH_WRAPPER_H //
// EOF //
