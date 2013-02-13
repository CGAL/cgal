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
#ifndef CGAL_DART_H
#define CGAL_DART_H 1

#include <CGAL/Compact_container.h>
#include <CGAL/assertions.h>
#include <bitset>

namespace CGAL {

  /** @file Dart.h
   * Definition of nD dart.
   */

  template < unsigned int d_, class Refs,
             class Items_, class Alloc_ >
  class Combinatorial_map_base;

  template<class Map, unsigned int i, unsigned int nmi>
  struct Remove_cell_functor;

  namespace internal {
    template <typename Map,unsigned int i>
    struct basic_link_beta_functor;

    template <typename CMap,unsigned int i>
    struct link_beta_functor;
  }

#define CGAL_BETAINV(i) (i>1?i:(i==1?0:1))

  /** Definition of nD dart.
   * The Dart class describes an nD dart (basic element of a
   * combinatorial map). A dart is composed with handle towards its neighbors,
   * a bitset containing Boolean marks, and handle towards enabled attributes.
   * n is the dimension of the space (2 for 2D, 3 for 3D...)
   * Refs the ref class
   */
  template <int d, typename Refs>
  struct Dart
  {
    template < unsigned int d_, class Refs_,
               class Items_, class Alloc_ >
    friend class Combinatorial_map_base;

    template <class T, class Alloc_>
    friend class Compact_container;

    template<class Map, unsigned int i, unsigned int nmi>
    friend struct Remove_cell_functor;

    template<class Map, unsigned int i>
    friend struct Contract_cell_functor;

    template <typename Map,unsigned int i>
    friend struct internal::link_beta_functor;

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

    /** Return if this dart is free for adimension.
     * @param i the dimension.
     * @return true iff the dart is linked with NULL for \em adimension.
     */
    template<unsigned int i>
    bool is_free() const
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i] == Refs::null_dart_handle;
    }
    bool is_free(unsigned int i) const
    {
      CGAL_assertion(i <= dimension);
      return mbeta[i] == Refs::null_dart_handle;
    }

    /** Return the highest dimension for which the dart is not free.
     * @return the dimension d such that the dart is not d-free but k-free for
     *         all k>d. -1 if the dart is free for all d in {0..n}
     */
    int highest_nonfree_dimension() const
    {
      for (int i=(int)dimension; i>=0; --i)
      { if ( !is_free(i) ) return i; }
      return -1;
    }

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

    /** Return a dart belonging to the same edge and to the second vertex
     * of the current edge (NULL if such a dart does not exist).
     * @return An handle to the opposite dart.
     */
    Dart_handle opposite()
    {
      for (unsigned int i = 2; i <= dimension; ++i)
        if (!is_free(i)) return beta(i);
      return NULL;
    }
    Dart_const_handle opposite() const
    {
      for (unsigned int i = 2; i <= dimension; ++i)
        if (!is_free(i)) return beta(i);
      return NULL;
    }

    /** Return a dart incident to the other extremity of the current edge,
     *  but contrary to opposite, non necessary to the same edge
     *  (NULL if such a dart does not exist).
     * @return An handle to the opposite dart.
     */
    Dart_handle other_extremity()
    {
      for (unsigned int i = 1; i <= dimension; ++i)
        if (!is_free(i)) return beta(i);
      return NULL;
    }
    Dart_const_handle other_extremity() const
    {
      for (unsigned int i = 1; i <= dimension; ++i)
        if (!is_free(i)) return beta(i);
      return NULL;
    }

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

    /** Link this dart with a given dart for a given dimension.
     * @param adart the dart to link with.
     * @param i the dimension.
     */
    template<unsigned int i>
    void basic_link_beta(Dart_handle adart)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(this!=&*Refs::null_dart_handle);
      mbeta[i] = adart;
    }
    void basic_link_beta(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(this!=&*Refs::null_dart_handle);
      mbeta[i] = adart;
    }

    /** Unlink this dart for a given dimension.
     * @param i the dimension.
     */
    template<unsigned int i>
    void unlink_beta()
    {
      CGAL_assertion(i <= dimension);
      mbeta[i] = Refs::null_dart_handle;
    }
    void unlink_beta(unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      mbeta[i] = Refs::null_dart_handle;
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

  protected:
    /** Default constructor: initialise marks and beta of this dart.
     * @param amarks the marks.
     */
    Dart(const std::bitset<NB_MARKS>& amarks) : mmarks(amarks)
    {
      for (unsigned int i = 0; i <= dimension; ++i)
        mbeta[i] = Refs::null_dart_handle;

      Helper::template Foreach_enabled_attributes<Init_attribute_functor>::
        run(this);
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

    /// Functor used to initialize all attributes to NULL.
    struct Init_attribute_functor
    {
      template <int i>
      static void run(Self* adart)
      { adart->template set_attribute<i>(NULL); }
    };

  public:
    void * for_compact_container() const
    { return mbeta[0].for_compact_container(); }
    void * & for_compact_container()
    { return mbeta[0].for_compact_container(); }

  protected:
    /// Beta for each dimension +1 (from 0 to dimension).
    Dart_handle mbeta[dimension+1];

    /// Attributes enabled
    typename Helper::Attribute_handles mattribute_handles;

    /// Values of Boolean marks.
    mutable std::bitset<NB_MARKS> mmarks;
  };

} // namespace CGAL

#endif // CGAL_DART_H //
// EOF //
