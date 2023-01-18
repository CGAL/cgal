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
#ifndef CGAL_CELL_ATTRIBUTE_H
#define CGAL_CELL_ATTRIBUTE_H 1

#include <CGAL/tags.h>
#include <CGAL/assertions.h>
#include <cstddef>
#include <CGAL/Info_for_cell_attribute.h>

namespace CGAL {

template <class, class, class, class>
class Compact_container;

template <class, class>
class Concurrent_compact_container;

template <class, class, class, class>
class Compact_container_with_index;

template<unsigned int, class, class>
class Combinatorial_map_storage_1;

template<unsigned int, class, class>
class Combinatorial_map_storage_with_index;

template<unsigned int, class, class>
class Generalized_map_storage_1;

template<unsigned int, class, class>
class Generalized_map_storage_with_index;

template<unsigned int, unsigned int, class, class, class>
class CMap_linear_cell_complex_storage_1;

template<unsigned int, unsigned int, class, class, class>
class CMap_linear_cell_complex_storage_with_index;

template<unsigned int, unsigned int, class, class, class>
class GMap_linear_cell_complex_storage_1;

template<unsigned int, unsigned int, class, class, class>
class GMap_linear_cell_complex_storage_with_index;

namespace internal {

template<class, class>
struct Init_id;

} // end namespace internal

  /** @file Cell_attribute.h
   * Definition of cell attribute, with or without info.
   */

  /// Id associated with a cell attribute
  template <class WithId>
  class Add_id
  {
  public:
    // Required to have "internal" property maps.
    std::size_t& id()
    { return m_id; }
    const std::size_t& id() const
    { return m_id; }

  protected:
    void set_id(std::size_t id)
    { m_id=id; }

  protected:
    std::size_t m_id; ///< id of the cell
  };

  /// If the tag WithId is false, we do not add id to cells.
  template <>
  class Add_id<Tag_false>
  {};

  /// Cell_attribute_without_info
  template <class Refs, class Tag=Tag_true, class OnMerge=Null_functor,
            class OnSplit=Null_functor, class WithID=Tag_false>
  class Cell_attribute_without_info;

  // Cell_attribute_without_info without dart support.
  template <class Refs, class OnMerge, class OnSplit, class WithID>
  class Cell_attribute_without_info<Refs, Tag_false,
                                    OnMerge, OnSplit, WithID>:
      public Add_id<WithID>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_with_index;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_with_index;

    template<class, class>
    friend struct internal::Init_id;

  public:
    typedef Tag_false                            Supports_cell_dart;

    typedef typename Refs::Dart_descriptor           Dart_descriptor;
    typedef typename Refs::Dart_const_descriptor     Dart_const_descriptor;
    typedef typename Refs::Alloc                 Alloc;

    typedef OnMerge On_merge;
    typedef OnSplit On_split;
    typedef WithID Has_id;
    using Type_for_compact_container=typename Refs::Type_for_compact_container;

    /// operator =
    Cell_attribute_without_info&
    operator=(const Cell_attribute_without_info& acell)
    {
      mrefcounting=acell.mrefcounting;
      m_for_cc=acell.m_for_cc;
      return *this;
    }

    /// Get the dart associated with the cell.
    Dart_descriptor dart() { return Refs::null_descriptor; }

    /// Get the dart associated with the cell.
    Dart_const_descriptor dart() const { return Refs::null_descriptor; }

    /// Set the dart associated with the cell.
    void set_dart(Dart_descriptor) {}

    /// Test if the cell is valid.
    /// For cell without dart, return always true.
    bool is_valid() const
    { return true; }

    bool operator==(const Cell_attribute_without_info&) const
    { return true; }

    bool operator!=(const Cell_attribute_without_info& other) const
    { return !operator==(other); }

  protected:
    /// Constructor without parameter.
    Cell_attribute_without_info(): mrefcounting(0), m_for_cc(Refs::null_descriptor)
    {}

    /// Copy constructor.
    Cell_attribute_without_info(const Cell_attribute_without_info& acell):
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
    std::size_t get_nb_refs() const
    { return mrefcounting; }

    Type_for_compact_container for_compact_container() const
    { return m_for_cc; }
    void for_compact_container(Type_for_compact_container p)
    { m_for_cc=p; }

  private:
    /// Reference counting: the number of darts linked to this cell.
    std::size_t                mrefcounting;
    Type_for_compact_container m_for_cc; // TODO better: this is memory consuming
    // TODO: or keep like that and never use an attribute without info and without dart !
  };

  /** Definition of cell attribute.
   * Cell_attribute defines what is a a cell. This is an object allowing to
   * link to a dart of the cell (when T is true).
   * The refs class must provide the type of Combinatorial_map used.
   */
  template <class Refs, class OnMerge, class OnSplit, class WithID>
  class Cell_attribute_without_info<Refs, Tag_true,
                                    OnMerge, OnSplit, WithID>: public Add_id<WithID>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_with_index;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_with_index;

    template<class, class>
    friend struct internal::Init_id;

  public:
    typedef Tag_true                             Supports_cell_dart;

    typedef typename Refs::Dart_descriptor           Dart_descriptor;
    typedef typename Refs::Dart_const_descriptor     Dart_const_descriptor;
    typedef typename Refs::Alloc                 Alloc;

    typedef OnMerge On_merge;
    typedef OnSplit On_split;
    typedef WithID Has_id;
    using Type_for_compact_container=typename Refs::Type_for_compact_container;

    /// operator =
    Cell_attribute_without_info&
    operator=(const Cell_attribute_without_info& acell)
    {
      mdart = acell.mdart;
      mrefcounting = acell.mrefcounting;
      return *this;
    }

    /// Get the dart associated with the cell.
    Dart_descriptor dart() { return mdart; }

    /// Get the dart associated with the cell.
    Dart_const_descriptor dart() const { return mdart; }

    /// Set the dart associated with the cell.
    void set_dart(Dart_descriptor adart) { mdart = adart; }

    /// Test if the cell is valid.
    /// A cell is valid if its dart is not null_descriptor.
    bool is_valid() const
    { return mdart!=Refs::null_descriptor; }

    bool operator==(const Cell_attribute_without_info&) const
    { return true; }

    bool operator!=(const Cell_attribute_without_info& other) const
    { return !operator==(other); }

  protected:
    /// Constructor without parameter.
    Cell_attribute_without_info() : mdart(Refs::null_descriptor),
                                    mrefcounting(0)
    {}

    /// Copy constructor.
    Cell_attribute_without_info(const Cell_attribute_without_info& acell):
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
    std::size_t get_nb_refs() const
    { return mrefcounting; }

    Type_for_compact_container for_compact_container() const
    { return mdart.for_compact_container(); }
    void for_compact_container(Type_for_compact_container p)
    { mdart.for_compact_container(p); }

  private:
    /// The dart descriptor associated with the cell.
    Dart_descriptor mdart;

    /// Reference counting: the number of darts linked to this cell.
    std::size_t mrefcounting;
  };

  /// Cell associated with an attribute, with or without info depending
  /// if Info==void.
  template <class Refs, class Info_=void, class Tag_=Tag_true,
            class OnMerge=Null_functor,
            class OnSplit=Null_functor,
            class WithID=Tag_false>
  class Cell_attribute;

  /// Specialization when Info==void.
  template <class Refs, class Tag_, class OnMerge, class OnSplit, class WithID>
  class Cell_attribute<Refs, void, Tag_, OnMerge, OnSplit, WithID> :
    public Cell_attribute_without_info<Refs, Tag_, OnMerge, OnSplit, WithID>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_with_index;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_with_index;

  public:
    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_descriptor       Dart_descriptor;
    typedef typename Refs::Dart_const_descriptor Dart_const_descriptor;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
    typedef void                             Info;

  protected:
    /// Default constructor.
    Cell_attribute()
    {}
  };

  /// Specialization when Info!=void.
  template <class Refs, class Info_, class Tag_,
            class OnMerge, class OnSplit, class WithID>
  class Cell_attribute :
    public Cell_attribute_without_info<Refs, Tag_, OnMerge, OnSplit, WithID>,
    public Info_for_cell_attribute<Info_>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_with_index;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, class, class>
    friend class Generalized_map_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class CMap_linear_cell_complex_storage_with_index;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class>
    friend class GMap_linear_cell_complex_storage_with_index;

  public:
    typedef Cell_attribute<Refs, Info_, Tag_, OnMerge, OnSplit, WithID> Self;

    typedef Tag_                             Supports_cell_dart;
    typedef typename Refs::Dart_descriptor       Dart_descriptor;
    typedef typename Refs::Dart_const_descriptor Dart_const_descriptor;
    typedef typename Refs::Alloc             Alloc;
    typedef OnMerge                          On_merge;
    typedef OnSplit                          On_split;
    typedef Info_                            Info;

    bool operator==(const Self& other) const
    { return this->info()==other.info(); }

    bool operator!=(const Self& other) const
    { return !operator==(other); }

  protected:
    /// Default constructor.
    Cell_attribute()
    {}

    /// Constructor with an info in parameter.
    Cell_attribute(const Info_& ainfo) :
      Info_for_cell_attribute<Info_>(ainfo)
    {}
  };

} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_H //
// EOF //
