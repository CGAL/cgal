// Copyright (c) 2014 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GMAP_DART_H
#define CGAL_GMAP_DART_H 1

#include <CGAL/Compact_container.h>
#include <CGAL/assertions.h>
#include <bitset>

namespace CGAL {

/** @file GMap_dart.h
 * Definition of nD dart for generalized maps.
 */

  namespace internal {
    template <typename Map,unsigned int i>
    struct basic_link_alpha_functor;

    template <typename CMap,unsigned int i>
    struct link_alpha_functor;
  }

/** Definition of nD dart.
 * The GMap_dart class describes an nD dart (basic element of a
 * generalized map). A dart is composed with handle towards its neighbors,
 * a bitset containing Boolean marks, and handle towards enabled attributes.
 * n is the dimension of the space (2 for 2D, 3 for 3D...)
 * Refs the ref class
 */
template <unsigned int d, typename Refs>
struct GMap_dart
{
  template < unsigned int, class, class, class, class >
  friend class Generalized_map_base;

  template < unsigned int, class, class >
  friend class Generalized_map_storage_1;

  template < unsigned int, class, class >
  friend class Generalized_map_storage_2;

  template < unsigned int, unsigned int, class, class, class >
  friend class GMap_linear_cell_complex_storage_1;

  template <class, class, class, class>
  friend class Compact_container;

  template<class, unsigned int, unsigned int>
  friend struct Remove_cell_functor;

  template<class, unsigned int>
  friend struct Contract_cell_functor;

  template <typename, unsigned int>
  friend struct internal::link_alpha_functor;

public:
  typedef GMap_dart<d,Refs>                Self;
  typedef typename Refs::Dart_handle       Dart_handle;
  typedef typename Refs::size_type         size_type;
  typedef typename Refs::Dart_const_handle Dart_const_handle;
  typedef typename Refs::Helper            Helper;
  /// Typedef for attributes
  template<int i>
  struct Attribute_handle: public Refs::template Attribute_handle<i>
  {};
  template<int i>
  struct Attribute_const_handle: public Refs::template Attribute_const_handle<i>
  {};

   /// The number of used marks.
   static const size_type NB_MARKS = Refs::NB_MARKS;

   /// The dimension of the generalized map.
   static const unsigned int dimension = d;

   /** Return the alpha of this dart for a given dimension.
    * @return alpha(\em i).
    */
  template<unsigned int i>
  Dart_handle alpha()
  {
    CGAL_assertion(i <= dimension);
    return malpha[i];
  }
  template<unsigned int i>
  Dart_const_handle alpha() const
  {
    CGAL_assertion(i <= dimension);
    return malpha[i];
  }
  Dart_handle alpha(unsigned int i)
  {
    CGAL_assertion(i <= dimension);
    return malpha[i];
  }
  Dart_const_handle alpha(unsigned int i) const
  {
    CGAL_assertion(i <= dimension);
    return malpha[i];
  }

  /// @return a handle on the i-attribute
  template<int i>
  typename Attribute_handle<i>::type attribute()
  {
    CGAL_static_assertion_msg( Helper::template Dimension_index<i>::value>=0,
                               "attribute<i> but i-attributes are disabled.");
    return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
      (mattribute_handles);
  }
  template<int i>
  typename Attribute_const_handle<i>::type attribute() const
  {
    CGAL_static_assertion_msg( Helper::template Dimension_index<i>::value>=0,
                               "attribute<i> but i-attributes are disabled.");
    return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
      (mattribute_handles);
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

protected:
    /** Default constructor: no real initialisation,
     *  because this is done in the generalized map class.
     */
    GMap_dart()
    {}

    /** Copy constructor:
     * @param adart a dart.
     */
    GMap_dart(const GMap_dart& adart) :
      mmarks(adart.mmarks),
      mattribute_handles(adart.mattribute_handles)
    {
      for (unsigned int i = 0; i <= dimension; ++i)
        malpha[i] = adart.malpha[i];
    }

public:
  void * for_compact_container() const
  { return malpha[0].for_compact_container(); }
  void * & for_compact_container()
  { return malpha[0].for_compact_container(); }

protected:
   /// Alpha for each dimension +1 (from 0 to dimension).
   Dart_handle malpha[dimension+1];

   /// Values of Boolean marks.
   mutable std::bitset<NB_MARKS> mmarks;

  /// Attributes enabled
   typename Helper::Attribute_handles mattribute_handles;
};

} // namespace CGAL

#endif // CGAL_GMAP_DART_H //
// EOF //
