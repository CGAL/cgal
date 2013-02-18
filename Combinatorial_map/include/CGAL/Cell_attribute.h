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
#ifndef CGAL_CELL_ATTRIBUTE_H
#define CGAL_CELL_ATTRIBUTE_H 1

#include <CGAL/Compact_container.h>

namespace CGAL {

  template < unsigned int d_, class Refs,
             class Items_, class Alloc_ >
  class Combinatorial_map_base;

  /** @file Cell_attribute.h
   * Definition of cell attribute, with or without info.
   */

  /// Info associated with a cell_attribute.
  template <typename Info>
  class Info_for_cell_attribute
  {
  public:
    /// Contructor without parameter.
    Info_for_cell_attribute()
    {}

    /// Contructor with an info in parameter.
    Info_for_cell_attribute(const Info& ainfo) : minfo(ainfo)
    {}

    /// Get the info associated with the cell_attribute.
    Info& info()
    { return minfo; }

    /// Get the info associated with the cell_attribute.
    const Info& info() const
    { return minfo; }

  protected:
    /// The info associated with the cell_attribute.
    Info minfo;
  };

  /// Cell_attribute_without_info
  template <class Refs, class Tag=Tag_true,
            class OnMerge=Null_functor, 
            class OnSplit=Null_functor>
  class Cell_attribute_without_info;

  // Cell_attribute_without_info without dart support.
  template <class Refs, class OnMerge, class OnSplit>
  class Cell_attribute_without_info<Refs, Tag_false, 
                                    OnMerge, OnSplit>
  {
    template < unsigned int d_, class Refs_,
               class Items_, class Alloc_ >
    friend class Combinatorial_map_base;
    
    template <int d, typename Refs_>
    friend struct Dart;

    template < unsigned int d_, class Refs_,
	       class Items_, class Alloc_ >
    friend class Generalized_map_base;
    
    template <int d, typename Refs_>
    friend struct GMap_dart;

    template <class T, class Alloc_>
    friend class Compact_container;

    template<typename CMap, unsigned int i, typename T>
    friend struct internal::Decrease_attribute_functor_run;

  public:
    typedef Tag_false                            Supports_cell_dart;

    typedef typename Refs::Dart_handle           Dart_handle;
    typedef typename Refs::Dart_const_handle     Dart_const_handle;
    typedef typename Refs::Alloc                 Alloc;

    typedef OnMerge On_merge;
    typedef OnSplit On_split;

    /// operator =
    /// We do nothing since we must not copy mrefcounting.
    Cell_attribute_without_info& 
    operator=(const Cell_attribute_without_info& /*acell*/)
    { return *this; }

    /// Get the dart associated with the cell.
    Dart_handle dart() { return NULL; }

    /// Get the dart associated with the cell.
    Dart_const_handle dart() const { return NULL; }

    /// Set the dart associated with the cell.
    void set_dart(Dart_handle) {}

    /// Test if the cell is valid.
    /// For cell without dart, return always true.
    bool is_valid() const
    { return true; }

  protected:
    /// Contructor without parameter.
    Cell_attribute_without_info(): mrefcounting(0)
    {}

    /// Copy contructor.
    Cell_attribute_without_info(const Cell_attribute_without_info&):
      mrefcounting(0)
    {}

    /// Increment the reference counting.
    void inc_nb_refs()
    { ++mrefcounting; }

    /// Decrement the reference counting.
    void dec_nb_refs()
    {
      CGAL_assertion( mrefcounting>0 );
      --mrefcounting; 
    }

    /// Get the reference counting.
    unsigned int get_nb_refs() const
    { return mrefcounting; }

  public:
    void * for_compact_container() const 
    { return vp; }
    void * & for_compact_container()       
    { return vp; }

  private:
    /// Reference counting: the number of darts linked to this cell.
    union
    {
      unsigned int mrefcounting;
      void        *vp;
    };
  };

  /** Definition of cell attribute.
   * Cell_attribute defines what is a a cell. This is an object allowing to
   * link to a dart of the cell (when T is true).
   * The refs class must provide the type of Combinatorial_map used.
   */
  template <class Refs, class OnMerge, class OnSplit>
  class Cell_attribute_without_info<Refs, Tag_true, 
                                    OnMerge, OnSplit>
  {
    template < unsigned int d_, class Refs_,
               class Items_, class Alloc_ >
    friend class Combinatorial_map_base;

    template <int d, typename Refs_>
    friend struct Dart;

    template < unsigned int d_, class Refs_,
	       class Items_, class Alloc_ >
    friend class Generalized_map_base;

    template <int d, typename Refs_>
    friend struct GMap_dart;

    template <class T, class Alloc_>
    friend class Compact_container;

    template<typename CMap, unsigned int i, typename T>
    friend struct internal::Decrease_attribute_functor_run;

  public:
    typedef Tag_true                             Supports_cell_dart;

    typedef typename Refs::Dart_handle           Dart_handle;
    typedef typename Refs::Dart_const_handle     Dart_const_handle;
    typedef typename Refs::Alloc                 Alloc;
    
    typedef OnMerge On_merge;
    typedef OnSplit On_split;

    /// operator =
    /// We must not copy mrefcounting.
    Cell_attribute_without_info& 
    operator=(const Cell_attribute_without_info& acell)
    {
      mdart = acell.mdart;
      return *this; 
    }

    /// Get the dart associated with the cell.
    Dart_handle dart() { return mdart; }

    /// Get the dart associated with the cell.
    Dart_const_handle dart() const { return mdart; }

    /// Set the dart associated with the cell.
    void set_dart(Dart_handle adart) { mdart = adart; }

    /// Test if the cell is valid.
    /// A cell is valid if its dart is not NULL.
    bool is_valid() const
    { return mdart!=NULL; }

  protected:
    /// Contructor without parameter.
    Cell_attribute_without_info() : mdart(NULL),
                                    mrefcounting(0)
    {}

    /// Copy contructor.
    Cell_attribute_without_info(const Cell_attribute_without_info& acell):
      mdart(acell.mdart),
      mrefcounting(0)
    {}

    /// Increment the reference counting.
    void inc_nb_refs()
    { ++mrefcounting; }

    /// Decrement the reference counting.
    void dec_nb_refs()
    {
      CGAL_assertion( mrefcounting>0 );
      --mrefcounting; 
    }

  public:
    /// Get the reference counting.
    unsigned int get_nb_refs() const
    { return mrefcounting; }

    void * for_compact_container() const
    { return mdart.for_compact_container(); }
    void * & for_compact_container()       
    { return mdart.for_compact_container(); }

  private:
    /// The dart handle associated with the cell.
    Dart_handle mdart;
   
    /// Reference counting: the number of darts linked to this cell.
    unsigned int mrefcounting;
  };

  /// Cell associated with an attribute, with or without info depending 
  /// if Info==void.
  template <class Refs, class Info_=void, class Tag_=Tag_true, 
            class OnMerge=Null_functor, 
            class OnSplit=Null_functor>
  class Cell_attribute;
  
  /// Specialization when Info==void.
  template <class Refs, class Tag_, 
            class OnMerge, class OnSplit>
  class Cell_attribute<Refs, void, Tag_, 
                       OnSplit, OnMerge> : 
    public Cell_attribute_without_info<Refs, Tag_, 
                                       OnSplit, OnMerge>
  {
    template < unsigned int d_, class Refs_,
               class Items_, class Alloc_ >
    friend class Combinatorial_map_base;
    
    template <class T, class Alloc_>
    friend class Compact_container;

  public:
    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
  };
  

  /// Specialization when Info!=void.
  template <class Refs, class Info_, class Tag_,
            class OnMerge, class OnSplit>
  class Cell_attribute :
    public Cell_attribute_without_info<Refs, Tag_,
                                       OnMerge, OnSplit>,
    public Info_for_cell_attribute<Info_>
  {
    template < unsigned int d_, class Refs_,
               class Items_, class Alloc_ >
    friend class Combinatorial_map_base;
    
    template <class T, class Alloc_>
    friend class Compact_container;

  public:
    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
    typedef Info_                            Info;

  protected:
    /// Default contructor.
    Cell_attribute()
    {}

    /// Contructor with an attribute in parameter.
    Cell_attribute(const Info_& ainfo) : Info_for_cell_attribute<Info_>(ainfo)
    {}
  };

} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_H //
// EOF //
