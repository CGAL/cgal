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
#ifndef CGAL_DART_H
#define CGAL_DART_H 1

#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/tuple.h>
#include <bitset>
#include <CGAL/Cell_attribute.h>

namespace CGAL {

  template <class, class, class, class>
  class Compact_container;

  template <class, class>
  class Concurrent_compact_container;

  template<unsigned int, class, class, class>
  class Combinatorial_map_storage_1;

  template<unsigned int, class, class, class>
  class Generalized_map_storage_1;

  template<unsigned int, unsigned int, class, class, class, class>
  class CMap_linear_cell_complex_storage_1;

  template<unsigned int, unsigned int, class, class, class, class>
  class GMap_linear_cell_complex_storage_1;

  namespace internal {
    template<class, class>
    struct Init_id;
  } // end namespace internal

  /** @file Dart.h
   * Definition of nD dart.
   */

  /** Definition of nD dart without information.
   * The_dart class describes an nD dart (basic element of a combinatorial or generalized map).
   * A dart is composed with handle towards its neighbors,
   * a bitset containing Boolean marks, and handle towards enabled attributes.
   * n is the dimension of the space (2 for 2D, 3 for 3D...)
   * Refs the ref class
   */
  template <unsigned int d, typename Refs, class WithId>
  struct Dart_without_info: public Add_id<WithId>
  {
  public:
    template<unsigned int, class, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, class, class, class>
    friend class Combinatorial_map_storage_2;
    
    template<unsigned int, unsigned int, class, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class, class>
    friend class Linear_cell_complex_storage_2;

    template<unsigned int, unsigned int, class, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template <class, class, class, class>
    friend class Compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template <class, class, class, class>
    friend class Compact_container_with_index_2;

    template <class, class>
    friend class Concurrent_compact_container;

    template<class, class>
    friend struct internal::Init_id;

    typedef Dart_without_info<d,Refs, WithId> Self;
    typedef typename Refs::Dart_handle        Dart_handle;
    typedef typename Refs::size_type          size_type;
    typedef typename Refs::Dart_const_handle  Dart_const_handle;
    typedef typename Refs::Helper             Helper;
    typedef WithId                            Has_id;

    /// Typedef for attributes
    template<int i>
    struct Attribute_handle: public Refs::template Attribute_handle<i>
    {};
    template<int i>
    struct Attribute_const_handle:
      public Refs::template Attribute_const_handle<i>
    {};

    /// The number of used marks.
    static const size_type NB_MARKS = Refs::NB_MARKS;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d;

    void * for_compact_container() const
    { return mf[0].for_compact_container(); }
    void for_compact_container(void *p)
    { mf[0].for_compact_container(p); }

    Dart_handle get_f(unsigned int i) const
    {
      assert(i<=dimension);
      return mf[i];
    }

  protected:
    /** Default constructor: no real initialisation,
     *  because this is done in the combinatorial map class.
     */
    Dart_without_info()
    {}

    /** Copy constructor:
     * @param adart a dart.
     */
    Dart_without_info(const Dart_without_info& adart) : mmarks(adart.mmarks),
    mattribute_handles(adart.mattribute_handles)
    {
      for (unsigned int i = 0; i <= dimension; ++i)
        mf[i] = adart.mf[i];
    }

    /** Return the mark value of a given mark number.
     * @param amark the mark number.
     * @return the value for this number.
     */
    bool get_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      return mmarks[amark];
    }

    /** Set the mark of a given mark number to a given value.
     * @param amark the mark number.
     * @param AValue the value.
     */
    void set_mark(size_type amark, bool avalue) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      mmarks.set(amark, avalue);
    }
    /** Flip the mark of a given mark number.
     * @param amark the mark number.
     */
    void flip_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      mmarks.flip(amark);
    }

    /** Return all the marks of this dart.
     * @return the marks.
     */
     std::bitset<NB_MARKS> get_marks() const
    { return mmarks; }

    /** Set simultaneously all the marks of this dart to a given value.
     * @param amarks the value of the marks.
     */
     void set_marks(const std::bitset<NB_MARKS>& amarks) const
    { mmarks = amarks; }

    /// @return a handle on the i-attribute
    template<int i>
    typename Attribute_handle<i>::type attribute()
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_handles);
    }
    template<int i>
    typename Attribute_const_handle<i>::type attribute() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_handles);
    }

  protected:
    /// Neighboors for each dimension +1 (from 0 to dimension).
    Dart_handle mf[dimension+1];

    /// Values of Boolean marks.
    mutable std::bitset<NB_MARKS> mmarks;

    /// Attributes enabled
    typename Helper::Attribute_handles mattribute_handles;
  };

  // Dart definition with an info;
  //  (there is a specialization below when Info_==void)
  template <unsigned int d, typename Refs, typename Info_=void,
            class WithID=Tag_false>
  struct Dart : public Dart_without_info<d, Refs, WithID>
  {
  public:
    template<unsigned int, class, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class, class>
    friend class Generalized_map_storage_1;

    template<unsigned int, unsigned int, class, class, class, class>
    friend class CMap_linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class, class>
    friend class GMap_linear_cell_complex_storage_1;

    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

    typedef Dart<d, Refs, Info_, WithID> Self;
    typedef Info_                        Info;

  protected:
    /** Default constructor: no real initialisation,
     *  because this is done in the combinatorial or generalized map class.
     */
    Dart()
    {}

    Dart(const Info_& info) : minfo(info)
    {}

    Info_& info()
    { return minfo; }
    const Info_& info() const
    { return minfo; }

  protected:
    Info minfo;
  };

  // Specialization of Dart class when info==void
  template <unsigned int d, typename Refs, class WithID>
  struct Dart<d, Refs, void, WithID> : public Dart_without_info<d, Refs, WithID>
  {
  public:
    typedef CGAL::Void Info;
  };
  namespace Index
  {
  template <unsigned int d, typename Refs>
  struct Dart
  {
    template < unsigned int, class, class, class, class >
    friend class Combinatorial_map_base;

    template<unsigned int, class, class>
    friend class Combinatorial_map_storage_1;

    template<unsigned int, class, class, class>
    friend class Combinatorial_map_storage_2;

    template<unsigned int, unsigned int, class, class, class>
    friend class Linear_cell_complex_storage_1;

    template<unsigned int, unsigned int, class, class, class, class>
    friend class Linear_cell_complex_storage_2;

    template <class, class, class, class>
    friend class Compact_container;

    template <class, class, class, class>
    friend class Compact_container_with_index;

    template <class, class, class, class>
    friend class Compact_container_with_index_2;

    template<class, unsigned int, unsigned int>
    friend struct Remove_cell_functor;

    template<class, unsigned int>
    friend struct Contract_cell_functor;

    template <typename,unsigned int>
    friend struct internal::link_beta_functor;

    template <typename, typename>
    friend struct internal::Reverse_orientation_of_map_functor;

    template <typename, typename>
    friend struct internal::Reverse_orientation_of_connected_component_functor;

  public:
    typedef Dart<d,Refs>                     Self;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::size_type         size_type;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Helper            Helper;
    /// Typedef for attributes
    template<int i>
    struct Attribute_handle: public Refs::template Attribute_handle<i>
    {};
    template<int i>
    struct Attribute_const_handle:
      public Refs::template Attribute_const_handle<i>
    {};

    /// The number of used marks.
    static const size_type NB_MARKS = Refs::NB_MARKS;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d;

    /** Return the beta of this dart for a given dimension.
     * @param i the dimension.
     * @return beta(\em i).
     */
    template<unsigned int i>
    Dart_handle beta()
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i];
    }
    Dart_handle beta(unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i];
    }
    template<unsigned int i>
    Dart_const_handle beta() const
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i];
    }
    Dart_const_handle beta(unsigned int i) const
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i];
    }

    /** Return the beta inverse of this dart for a given dimension.
     * @param i the dimension.
     * @return beta^{-1}(\em i).
     */
    template<unsigned int i>
    Dart_handle beta_inv()
    { return beta<CGAL_BETAINV(i)>(); }
    Dart_handle beta_inv(unsigned int i)
    { return beta(CGAL_BETAINV(i)); }
    template<unsigned int i>
    Dart_const_handle beta_inv() const
    { return beta<CGAL_BETAINV(i)>(); }
    Dart_const_handle beta_inv(unsigned int i) const
    { return beta(CGAL_BETAINV(i)); }

    /// @return a handle on the i-attribute
    template<int i>
    typename Attribute_handle<i>::type attribute()
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_handles);
    }
    template<int i>
    typename Attribute_const_handle<i>::type attribute() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_handles);
    }

    /// Set the handle on the i th attribute
    template<int i>
    void set_attribute( typename Attribute_handle<i>::type ahandle )
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "set_attribute<i> called but i-attributes are disabled.");
      CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_handles) = ahandle;
      if (ahandle!=NULL) ahandle->inc_nb_refs();
    }

    /** Return the mark value of a given mark number.
     * @param amark the mark number.
     * @return the value for this number.
     */
    bool get_mark(int amark) const
    {
      CGAL_assertion(amark>=0 && (size_type)amark<NB_MARKS);
      return mmarks[(size_type)amark];
    }

    /** Set the mark of a given mark number to a given value.
     * @param amark the mark number.
     * @param AValue the value.
     */
    void set_mark(int amark, bool avalue) const
    {
      CGAL_assertion(amark>=0 && (size_type)amark<NB_MARKS);
      mmarks.set((size_type)amark, avalue);
    }

    /** Flip the mark of a given mark number.
     * @param amark the mark number.
     */
    void flip_mark(int amark) const
    {
      CGAL_assertion(amark>=0 && (size_type)amark<NB_MARKS);
      mmarks.flip((size_type)amark);
    }

    /** Return all the marks of this dart.
     * @return the marks.
     */
     std::bitset<NB_MARKS> get_marks() const
    { return mmarks; }

    /** Set simultaneously all the marks of this dart to a given value.
     * @param amarks the value of the marks.
     */
    void set_marks(const std::bitset<NB_MARKS>& amarks) const
    { mmarks = amarks; }

  protected:
    /** Default constructor: no real initialisation,
     *  because this is done in the combinatorial map class.
     */
    Dart()
    {}

    /** Copy constructor:
     * @param adart a dart.
     */
    Dart(const Dart& adart) : mmarks(adart.mmarks),
      mattribute_handles(adart.mattribute_handles)
    {
      for (unsigned int i = 0; i <= dimension; ++i)
        mbeta[i] = adart.mbeta[i];
    }

  public:
    size_type for_compact_container_with_index() const
    { return mbeta[0].for_compact_container_with_index(); }
    size_type& for_compact_container_with_index()
    { return mbeta[0].for_compact_container_with_index(); }

  protected:
    /// Beta for each dimension +1 (from 0 to dimension).
    Dart_handle mbeta[dimension+1];

    /// Values of Boolean marks.
    mutable std::bitset<NB_MARKS> mmarks;

    /// Attributes enabled
    typename Helper::Attribute_handles mattribute_handles;
  };
  } // namespace Index

} // namespace CGAL

#endif // CGAL_DART_H //
// EOF //
