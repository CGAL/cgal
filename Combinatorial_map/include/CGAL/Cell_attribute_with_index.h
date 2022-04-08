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
#ifndef CGAL_CELL_ATTRIBUTE_WITH_INDEX_H
#define CGAL_CELL_ATTRIBUTE_WITH_INDEX_H 1

#include <CGAL/tags.h>
#include <CGAL/assertions.h>
#include <CGAL/Info_for_cell_attribute.h>
#include <cstddef>

namespace CGAL {

  template <class, class, class, class>
  class Compact_container_with_index_2;

  template<unsigned int, class, class, class>
  class Combinatorial_map_storage_2;

  namespace Index
  {
  // Versions to use with containers using index
  /// Cell_attribute_without_info
  template <class Refs, class Tag=Tag_true, class OnMerge=Null_functor,
            class OnSplit=Null_functor>
  class Cell_attribute_without_info;

  // Cell_attribute_without_info without dart support.
  template <class Refs, class OnMerge, class OnSplit>
  class Cell_attribute_without_info<Refs, Tag_false,
                                    OnMerge, OnSplit>
  {
    template<unsigned int, class, class, class>
    friend class CGAL::Combinatorial_map_storage_2;

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index_2;

  public:
    typedef Tag_false                        Supports_cell_dart;

    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef CGAL::Tag_false                  Has_id;

    typedef OnMerge On_merge;
    typedef OnSplit On_split;

    /// operator =
    Cell_attribute_without_info&
    operator=(const Cell_attribute_without_info& acell)
    {
      mrefcounting = acell.mrefcounting;
      return *this;
    }

    /// Get the dart associated with the cell.
    Dart_handle dart() { return Refs::null_handle; }

    /// Get the dart associated with the cell.
    Dart_const_handle dart() const { return Refs::null_handle; }

    /// Set the dart associated with the cell.
    void set_dart(Dart_handle) {}

    /// Test if the cell is valid.
    /// For cell without dart, return always true.
    bool is_valid() const
    { return true; }

    bool operator==(const Cell_attribute_without_info&) const
    { return true; }

    bool operator!=(const Cell_attribute_without_info& other) const
    { return !operator==(other); }

    // protected:
    /// Contructor without parameter.
    Cell_attribute_without_info(): mrefcounting(0)
    {}

    /// Copy contructor.
    Cell_attribute_without_info
    (const Cell_attribute_without_info& acell):
      mrefcounting(acell.mrefcounting)
    {}

  protected:
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
    typename Refs::size_type get_nb_refs() const
    { return mrefcounting; }

    typename Refs::size_type for_compact_container_with_index() const
    { return mrefcounting; }
    typename Refs::size_type & for_compact_container_with_index()
    { return mrefcounting; }

  private:
    /// Reference counting: the number of darts linked to this cell.
    typename Refs::size_type mrefcounting;
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
    template<unsigned int, class, class, class>
    friend class CGAL::Combinatorial_map_storage_2;

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index;

  public:
    typedef Tag_true                         Supports_cell_dart;

    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef CGAL::Tag_false                  Has_id;

    typedef OnMerge On_merge;
    typedef OnSplit On_split;

    /// operator =
    Cell_attribute_without_info&
    operator=(const Cell_attribute_without_info& acell)
    {
      mdart = acell.mdart;
      mrefcounting = acell.mrefcounting;
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
    { return mdart!=Refs::null_handle; }

    bool operator==(const Cell_attribute_without_info&) const
    { return true; }

    bool operator!=(const Cell_attribute_without_info& other) const
    { return !operator==(other); }

  protected:
    /// Contructor without parameter.
    Cell_attribute_without_info() : mdart(Refs::null_handle),
                                    mrefcounting(0)
    {}

    /// Copy contructor.
    Cell_attribute_without_info
    (const Cell_attribute_without_info& acell):
      mdart(acell.mdart),
      mrefcounting(acell.mrefcounting)
    {}

  protected:
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
    typename Refs::size_type get_nb_refs() const
    { return mrefcounting; }

    typename Refs::size_type for_compact_container_with_index() const
    { return mdart.for_compact_container_with_index(); }
    typename Refs::size_type & for_compact_container_with_index()
    { return mdart.for_compact_container_with_index(); }

  private:
    /// Reference counting: the number of darts linked to this cell.
    std::size_t mrefcounting;

    /// The dart handle associated with the cell.
    Dart_handle mdart;
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
    template < unsigned int, class, class, class, class >
    friend class CGAL::Combinatorial_map_base;

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index_2;

  public:
    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
    typedef void                             Info;

    //  protected:
    /// Default contructor.
    Cell_attribute()
    {}
  };

  /// Specialization when Info!=void.
  template <class Refs, class Info_, class Tag_,
            class OnMerge, class OnSplit>
  class Cell_attribute :
    public Cell_attribute_without_info<Refs, Tag_,
                                       OnMerge, OnSplit>,
    public CGAL::Info_for_cell_attribute<Info_>
  {
    template < unsigned int, class, class, class, class >
    friend class CGAL::Combinatorial_map_base;

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index_2;

  public:
    typedef Cell_attribute<Refs, Info_, Tag_, OnMerge, OnSplit> Self;

    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
    typedef Info_                            Info;

    bool operator==(const Self& other) const
    { return this->info()==other.info(); }

    bool operator!=(const Self& other) const
    { return !operator==(other); }

  protected:
    /// Default contructor.
    Cell_attribute()
    {}

    /// Contructor with an info in parameter.
    Cell_attribute(const Info_& ainfo) :
      Info_for_cell_attribute<Info_>(ainfo)
    {}
  };
  } // namespace Index

} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_WITH_INDEX_H //
// EOF //
