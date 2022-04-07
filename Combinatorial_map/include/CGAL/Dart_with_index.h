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
#ifndef CGAL_DART_WITH_INDEX_H
#define CGAL_DART_WITH_INDEX_H 1

#include <CGAL/assertions.h>
#include <CGAL/tags.h>
#include <CGAL/tuple.h>
#include <bitset>
#include <CGAL/Cell_attribute.h>

namespace CGAL {

  namespace Index
  {
  template <unsigned int d, typename Refs>
  struct Dart_without_info
  {
    template<unsigned int, class, class, class>
    friend class CGAL::Combinatorial_map_storage_2;

    // TODO template<unsigned int, unsigned int, class, class, class, class>
    // friend class CGAL::CMap_linear_cell_complex_storage_2;

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index_2;

    typedef Dart_without_info<d,Refs>        Self;
    typedef typename Refs::Dart_handle       Dart_handle;
    typedef typename Refs::size_type         size_type;
    typedef typename Refs::Dart_const_handle Dart_const_handle;
    typedef typename Refs::Helper            Helper;
    typedef CGAL::Tag_false                  Has_id;

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

    size_type for_compact_container_with_index() const
    { return mf[0].for_compact_container_with_index(); }
    size_type& for_compact_container_with_index()
    { return mf[0].for_compact_container_with_index(); }

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
  template <unsigned int d, typename Refs, typename Info_=void>
  struct Dart : public Dart_without_info<d, Refs>
  {
  public:
    template<unsigned int, class, class, class>
    friend class CGAL::Combinatorial_map_storage_2;

    /* TODO template<unsigned int, unsigned int, class, class, class, class>
    friend class CMap_linear_cell_complex_storage_2; */

    template <class, class, class, class>
    friend class CGAL::Compact_container_with_index_2;

    typedef Dart<d, Refs, Info_> Self;
    typedef Info_                Info;

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
  template <unsigned int d, typename Refs>
  struct Dart<d, Refs, void> : public Dart_without_info<d, Refs>
  {
  public:
    typedef CGAL::Void Info;
  };

  } // namespace Index

} // namespace CGAL

#endif // CGAL_DART_WITH_INDEX_H //
// EOF //
