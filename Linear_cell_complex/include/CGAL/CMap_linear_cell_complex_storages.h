// Copyright (c) 2013 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_CMAP_LINEAR_CELL_COMPLEX_STORAGES_H
#define CGAL_CMAP_LINEAR_CELL_COMPLEX_STORAGES_H 1

#include <CGAL/Compact_container.h>
#include <CGAL/Dart.h>
#include <CGAL/Handle_hash_function.h>
#include <bitset>

#include <boost/config.hpp>
#if  (BOOST_GCC >= 40900)
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Warray-bounds\"")
#endif

namespace CGAL {

  namespace internal {
    template <typename M>
    struct Combinatorial_map_helper;
  }

  /** @file CMap_linear_cell_complex_storages.h
   * Definition of storages for dD Linear cell complex for combinatorial maps.
   */

  // Storage of darts with compact container, beta with handles
  // Copy of Combinatorial_map_storage_1 and add new types related
  // to geometry (not possible to inherith because we use Self type
  // as template parameter of Dart_wrapper. If we inherit, Self is not
  // the correct type).
  template<unsigned int d_, unsigned int ambient_dim,
           class Traits_, class Items_, class Alloc_ >
  class CMap_linear_cell_complex_storage_1
  {
  public:
    typedef typename Traits_::Point  Point;
    typedef typename Traits_::Vector Vector;
    typedef typename Traits_::FT     FT;

    typedef CMap_linear_cell_complex_storage_1<d_, ambient_dim, Traits_,
    Items_, Alloc_> Self;
    typedef CGAL::Tag_false Use_index;

    typedef internal::Combinatorial_map_helper<Self>      Helper;

    typedef typename Items_::template Dart_wrapper<Self>  Dart_wrapper;

#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
    typedef typename Dart_wrapper::Dart                   Dart;
#else
    typedef typename internal::template Get_dart_info<Dart_wrapper>::type
                                                           Dart_info;
    typedef typename internal::template Get_darts_with_id<Dart_wrapper>::type
                                                           Darts_with_id;
    typedef CGAL::Dart<d_, Self, Dart_info, Darts_with_id> Dart;
#endif

#ifdef CGAL_CXX11
    typedef std::allocator_traits<Alloc_> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Dart> Dart_allocator;
#else
    typedef typename Alloc_::template rebind<Dart>::other  Dart_allocator;
#endif

    typedef Compact_container<Dart, Dart_allocator>        Dart_container;

    typedef typename Dart_container::iterator              Dart_handle;
    typedef typename Dart_container::const_iterator        Dart_const_handle;
    typedef typename Dart_container::size_type             size_type;

    typedef CGAL::Void* Null_handle_type;
    static const Null_handle_type null_handle;

    typedef Items_ Items;
    typedef Alloc_ Alloc;

#ifdef CGAL_CXX11
    template <typename T>
    struct Container_for_attributes :
      public Compact_container<T, typename std::allocator_traits<Alloc_>::template rebind_alloc<T> >
    {};
#else
    template <typename T>
    struct Container_for_attributes :
        public Compact_container<T, typename Alloc_::template rebind<T>::other>
    {};
#endif
    
    /// Typedef for attributes
    typedef typename internal::template Get_attributes_tuple<Dart_wrapper>::type
                                   Attributes;

    template<int i>
    struct Attribute_type: public Helper::template Attribute_type<i>
    {};
    template<int i>
    struct Attribute_handle: public Helper::template Attribute_handle<i>
    {};
    template<int i>
    struct Attribute_const_handle:
      public Helper::template Attribute_const_handle<i>
    {};
    template<int i>
    struct Attribute_range: public Helper::template Attribute_range<i>
    {};
    template<int i>
    struct Attribute_const_range:
      public Helper::template Attribute_const_range<i>
    {};

    typedef typename Attribute_type<0>::type Vertex_attribute;
    typedef typename Attribute_handle<0>::type Vertex_attribute_handle;
    typedef typename Attribute_const_handle<0>::type
    Vertex_attribute_const_handle;

    typedef typename Attribute_range<0>::type Vertex_attribute_range;
    typedef typename Attribute_const_range<0>::type
    Vertex_attribute_const_range;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d_;

    typedef Handle_hash_function Hash_function;

    // Init
    void init_storage()
    {
      // emplace null_dart; initialized in Combinatorial_map class
      null_dart_handle = mnull_dart_container.emplace();
    }

   /** Return if this dart is free for adimension.
     * @param dh a dart handle
     * @param i the dimension.
     * @return true iff dh is linked with NULL for \em adimension.
     */
    template<unsigned int i>
    bool is_free(Dart_const_handle dh) const
    {
      CGAL_assertion( dh!=NULL );
      CGAL_assertion(i <= dimension);
      return dh->mf[i]==null_dart_handle;
    }
    bool is_free(Dart_const_handle dh, unsigned int i) const
    {
      CGAL_assertion( dh!=NULL );
      CGAL_assertion(i <= dimension);
      return dh->mf[i]==null_dart_handle;
    }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(Dart_const_handle ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      CGAL_assertion( ADart!=NULL );
      ADart->set_marks(amarks);
    }
    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(Dart_const_handle ADart) const
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->get_marks();
    }
    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(Dart_const_handle ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_handle ADart, size_type amark, bool avalue) const
    {
      CGAL_assertion( ADart!=NULL );
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_handle ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=NULL );
      ADart->flip_mark(amark);
    }

    // Access to beta maps
    Dart_handle get_beta(Dart_handle ADart, int B1)
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return ADart->mf[B1];
    }
    Dart_const_handle get_beta(Dart_const_handle ADart, int B1) const
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }
    template<int B1>
    Dart_handle get_beta(Dart_handle ADart)
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }
    template<int B1>
    Dart_const_handle get_beta(Dart_const_handle ADart) const
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }

    // return a handle on the i-attribute
    template<unsigned int i>
    typename Attribute_handle<i>::type attribute(Dart_handle ADart)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (ADart->mattribute_handles);
    }
    template<unsigned int i>
    typename Attribute_const_handle<i>::type
    attribute(Dart_const_handle ADart) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (ADart->mattribute_handles);
    }

    // Copy a given attribute
    template<unsigned int i>
    typename Attribute_handle<i>::type copy_attribute
    (typename Attribute_const_handle<i>::type ah)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "copy_attribute<i> called but i-attributes are disabled.");
      typename Attribute_handle<i>::type res=
        CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(*ah);
      this->template init_attribute_ref_counting<i>(res);
      return res;
    }

    // Test if a given attribute is valid
    template<unsigned int i>
    bool is_valid_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=NULL );
      return ah->is_valid();
    }
    
    // accessors and modifiers to the attribute ref counting given its handle
    template<unsigned int i>
    std::size_t get_attribute_ref_counting
    (typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=NULL );
      return ah->get_nb_refs();
    }
    template<unsigned int i>
    void init_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      ah->mrefcounting=0;
    }
    template<unsigned int i>
    void inc_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      ah->inc_nb_refs();
    }
    template<unsigned int i>
    void dec_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      ah->dec_nb_refs();
    }

    // get the attribute given its handle
    template<unsigned int i>
    typename Attribute_type<i>::type&
    get_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      return *ah;
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=NULL );
      return *ah;
    }

    // Get the dart of the given attribute
    template<unsigned int i>
    Dart_handle dart_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      return ah->dart();
    }
    template<unsigned int i>
    Dart_const_handle
    dart_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=NULL );
      return ah->dart();
    }

    // Set the dart of the given attribute
    template<unsigned int i>
    void set_dart_of_attribute(typename Attribute_handle<i>::type ah,
                               Dart_handle adart)
    {
      CGAL_assertion( ah!=NULL );
      ah->set_dart(adart);
    }

#if !defined(CGAL_CMAP_DART_DEPRECATED) || defined(CGAL_NO_DEPRECATED_CODE)
    // Get the information associated with a given dart
    Dart_info& info(Dart_handle adart)
    { return adart->info(); }
    const Dart_info& info(Dart_const_handle adart) const
    { return adart->info(); }
#endif

    // Get the info of the given attribute
    template<unsigned int i>
    typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=NULL );
      return ah->info();
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=NULL );
      return ah->info();
    }

    // Get the info of the i-cell attribute associated with the given dart
    template<unsigned int i>
    typename Attribute_type<i>::type::Info & info(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL );
      CGAL_assertion( attribute<i>(adart)!=NULL );
      return info_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info(Dart_const_handle adart) const
    {
      CGAL_assertion( adart!=NULL );
      CGAL_assertion( attribute<i>(adart)!=NULL );
      return info_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the i-cell attribute associated with the given dart
    template<unsigned int i>
    Dart_handle dart(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL );
      CGAL_assertion( attribute<i>(adart)!=NULL );
      return dart_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    Dart_const_handle dart(Dart_const_handle adart) const
    {
      CGAL_assertion( adart!=NULL );
      CGAL_assertion( attribute<i>(adart)!=NULL );
      return dart_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the given 0-attribute
    Point & point_of_vertex_attribute(typename Attribute_handle<0>::type vh)
    {
      CGAL_assertion( vh!=NULL );
      return get_attribute<0>(vh).point();
    }

    const Point & point_of_vertex_attribute
    (typename Attribute_const_handle<0>::type vh) const
    {
      CGAL_assertion( vh!=NULL );
      return get_attribute<0>(vh).point();
    }

    // Debug function
    void display_dart(Dart_const_handle ADart) const
    { std::cout<<mdarts.index(ADart); }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_handle<i>::type ah) const
    { std::cout<< CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).index(ah); }

  protected:
    // Set the handle on the i th attribute
    template<unsigned int i>
    void basic_set_dart_attribute(Dart_handle dh,
                                  typename Attribute_handle<i>::type ah)
    {
      CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (dh->mattribute_handles) = ah;
    }

    /** Link a dart with a given dart for a given dimension.
     * @param adart the dart to link.
     * @param adart2 the dart to link with.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_link_beta(Dart_handle adart, Dart_handle adart2)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=NULL && adart2!=NULL);
      CGAL_assertion(adart!=null_dart_handle);
      adart->mf[i] = adart2;
    }
    void dart_link_beta(Dart_handle adart, Dart_handle adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=NULL && adart2!=NULL);
      CGAL_assertion(adart!=null_dart_handle);
      adart->mf[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_beta(Dart_handle adart)
    {
      CGAL_assertion(adart!=NULL && i <= dimension);
      adart->mf[i] = null_dart_handle;
    }
    void dart_unlink_beta(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(adart!=NULL && i <= dimension);
      adart->mf[i] = null_dart_handle;
    }

  public:
    /// Void dart. A dart d is i-free if beta_i(d)=null_dart_handle.
    Dart_handle null_dart_handle; // Todo Dart_const_handle ??

  protected:
    /// Dart container.
    Dart_container mdarts;

    /// Container for the null_dart_handle.
    Dart_container mnull_dart_container;

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

  /// null_handle
  template<unsigned int d_, unsigned int ambient_dim,
           class Traits_, class Items_, class Alloc_ >
  const typename CMap_linear_cell_complex_storage_1<d_, ambient_dim, Traits_,
                                         Items_, Alloc_>::Null_handle_type
  CMap_linear_cell_complex_storage_1<d_, ambient_dim, Traits_,
                                Items_, Alloc_>::null_handle = NULL;

} // namespace CGAL

#if  (BOOST_GCC >= 40900)
 _Pragma("GCC diagnostic pop")
#endif

#endif // CGAL_CMAP_LINEAR_CELL_COMPLEX_STORAGES_H //
// EOF //
