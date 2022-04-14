// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_STORAGES_WITH_INDEX_H
#define CGAL_GENERALIZED_MAP_STORAGES_WITH_INDEX_H 1

#include <CGAL/Compact_container_with_index_2.h>
#include <CGAL/Dart.h>
#include <bitset>

#include <boost/config.hpp>
#if defined(BOOST_GCC)
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Warray-bounds\"")
#endif

namespace CGAL {

  namespace internal {
    template <typename M>
    struct Combinatorial_map_helper;

    template<typename Concurrent_tag, class T, class Alloc_>
    struct Container_type;
  }

  // Storage of darts with compact container, alpha with indices
  template<unsigned int d_, class Items_, class Alloc_>
  class Generalized_map_storage_2
  {
  public:
    using Self=Generalized_map_storage_2<d_, Items_, Alloc_>;
    using Use_index=CGAL::Tag_true;
    using Concurrent_tag=typename internal::Get_concurrent_tag<Items_>::type;

    typedef internal::Combinatorial_map_helper<Self> Helper;

    using Index_type=typename internal::Get_index_type<Items_>::type;

    typedef typename Items_::template Dart_wrapper<Self>  Dart_wrapper;

    typedef typename internal::template Get_dart_info<Dart_wrapper>::type
                                                           Dart_info;
    typedef CGAL::Dart<d_, Self, Dart_info> Dart;

    typedef std::allocator_traits<Alloc_> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Dart> Dart_allocator;

    typedef Compact_container_with_index_2<Dart,Dart_allocator,
    Multiply_by_two_policy_for_cc_with_size<64>, Index_type>
    Dart_container;

    typedef typename Dart_container::Index Dart_index;

    // Definition of old types, for backward compatibility.
    typedef Dart_index Dart_handle;
    typedef Dart_index Dart_const_handle;
    typedef typename Dart_container::size_type size_type;

    typedef Dart_index Null_handle_type;
    static Null_handle_type null_handle;

    using Type_for_compact_container=Index_type;

    typedef Items_ Items;
    typedef Alloc_ Alloc;
    template <typename T>
    struct Container_for_attributes : public
        Compact_container_with_index_2<T,
        typename Alloc_::template rebind<T>::other,
        Multiply_by_two_policy_for_cc_with_size<64>, size_type >
    {};
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

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /// The dimension of the generalized map.
    static const unsigned int dimension = d_;

    typedef internal::Index_hash_function Hash_function;

    typedef Dart_container       Dart_range;
    typedef const Dart_container Dart_const_range;
    /// @return a Dart_range (range through all the darts of the map).
    Dart_range& darts()             { return mdarts;}
    Dart_const_range& darts() const { return mdarts; }
    
    // Init
    void init_storage()
    {}

    void clear_storage()
    {}
    
    /** Test if the map is empty.
     *  @return true iff the map is empty.
     */
    bool is_empty() const
    { return mdarts.empty(); }

    /// @return the number of darts.
    size_type number_of_darts() const
    { return mdarts.size(); }

    /** Return if this dart is free for adimension.
     * @param dh a dart handle
     * @param i the dimension.
     * @return true iff dh is linked with itself for \em adimension.
     */
    template<unsigned int i>
    bool is_free(Dart_const_handle dh) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mf[i]==dh;
    }
    bool is_free(Dart_const_handle dh, unsigned int i) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mf[i]==dh;
    }
    bool is_perforated(Dart_const_handle /*dh*/) const
    { return false; }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(Dart_const_handle ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      mdarts[ADart].set_marks(amarks);
    }
    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(Dart_const_handle ADart) const
    {
      return mdarts[ADart].get_marks();
    }
    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(Dart_const_handle ADart, size_type amark) const
    {
      return mdarts[ADart].get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_handle ADart, size_type amark, bool avalue) const
    {
      mdarts[ADart].set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_handle ADart, size_type amark) const
    {
      mdarts[ADart].flip_mark(amark);
    }

    // Access to alpha maps
    Dart_handle get_alpha(Dart_handle ADart, int B1)
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return mdarts[ADart].mf[B1];
    }
    Dart_const_handle get_alpha(Dart_const_handle ADart, int B1) const
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }
    template<int B1>
    Dart_handle get_alpha(Dart_handle ADart)
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }
    template<int B1>
    Dart_const_handle get_alpha(Dart_const_handle ADart) const
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }

    // return a handle on the i-attribute
    template<unsigned int i>
    typename Attribute_handle<i>::type attribute(Dart_handle ADart)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_handles);
    }
    template<unsigned int i>
    typename Attribute_const_handle<i>::type
    attribute(Dart_const_handle ADart) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_handles);
    }

    // Copy a given attribute
    template<unsigned int i>
    typename Attribute_handle<i>::type copy_attribute
    (typename Attribute_const_handle<i>::type ah)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "copy_attribute<i> called but i-attributes are disabled.");
      typename Attribute_handle<i>::type res=
        std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(get_attribute<i>(ah));
      this->template init_attribute_ref_counting<i>(res);
      return res;
    }

    // Test if a given attribute is valid
    template<unsigned int i>
    bool is_valid_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      return get_attribute<i>(ah).is_valid();
    }

    // accessors and modifiers to the attribute ref counting given its handle
    template<unsigned int i>
    std::size_t get_attribute_ref_counting
    (typename Attribute_const_handle<i>::type ah) const
    {
      return get_attribute<i>(ah).get_nb_refs();
    }
    template<unsigned int i>
    void init_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      get_attribute<i>(ah).mrefcounting=0;
    }
    template<unsigned int i>
    void inc_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      get_attribute<i>(ah).inc_nb_refs();
    }
    template<unsigned int i>
    void dec_attribute_ref_counting(typename Attribute_handle<i>::type ah)
    {
      get_attribute<i>(ah).dec_nb_refs();
    }

    // get the attribute given its index
    template<unsigned int i>
    typename Attribute_type<i>::type&
    get_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }

    // Get the dart of the given attribute
    template<unsigned int i>
    Dart_handle dart_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return get_attribute<i>(ah).dart();
    }
    template<unsigned int i>
    Dart_const_handle
    dart_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return get_attribute<i>(ah).dart();
    }

    // Set the dart of the given attribute
    template<unsigned int i>
    void set_dart_of_attribute(typename Attribute_handle<i>::type ah,
                               Dart_handle adart)
    {
      CGAL_assertion( ah!=null_handle );
      get_attribute<i>(ah).set_dart(adart);
    }

    // Get the information associated with a given dart
    Dart_info& info(Dart_handle adart)
    { return mdarts[adart].info(); }
    const Dart_info& info(Dart_const_handle adart) const
    { return mdarts[adart].info(); }

    // Get the info of the given attribute
    template<unsigned int i>
    typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return get_attribute<i>(ah).info();
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return get_attribute<i>(ah).info();
    }

    // Get the info of the i-cell attribute associated with the given dart
    template<unsigned int i>
    typename Attribute_type<i>::type::Info & info(Dart_handle adart)
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( this->template attribute<i>(adart)!=null_handle );
      return info_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info(Dart_const_handle adart) const
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( attribute<i>(adart)!=null_handle );
      return info_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the i-cell attribute associated with the given dart
    template<unsigned int i>
    Dart_handle dart(Dart_handle adart)
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( attribute<i>(adart)!=null_handle );
      return dart_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    Dart_const_handle dart(Dart_const_handle adart) const
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( attribute<i>(adart)!=null_handle );
      return dart_of_attribute<i>(attribute<i>(adart));
    }

    // Debug function
    void display_dart(Dart_const_handle ADart) const
    { std::cout<<ADart; }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_handle<i>::type ah) const
    { std::cout<<ah; }

  protected:
    // Set the handle on the i th attribute
    template<unsigned int i>
    void basic_set_dart_attribute(Dart_handle dh,
                                  typename Attribute_handle<i>::type ah)
    {
      std::get<Helper::template Dimension_index<i>::value>
        (mdarts[dh].mattribute_handles) = ah;
    }

    /** Link a dart with a given dart for a given dimension.
     * @param adart the dart to link.
     * @param adart2 the dart to link with.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_link_alpha(Dart_handle adart, Dart_handle adart2)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=null_handle && adart2!=null_handle);
      mdarts[adart].mf[i] = adart2;
    }
    void dart_link_alpha(Dart_handle adart, Dart_handle adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=null_handle && adart2!=null_handle);
      mdarts[adart].mf[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_alpha(Dart_handle adart)
    {
      CGAL_assertion(adart!=null_handle && i <= dimension);
      mdarts[adart].mf[i] = adart;
    }
    void dart_unlink_alpha(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(adart!=null_handle && i <= dimension);
      mdarts[adart].mf[i] = adart;
    }

  protected:
    Dart_handle null_dart_handle; // To be compatible with combinatorial map

    /// Dart container.
    Dart_container mdarts;
    
    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

  /// null_handle
  template<unsigned int d_, class Items_, class Alloc_>
  typename Generalized_map_storage_2<d_, Items_, Alloc_>::Null_handle_type
      Generalized_map_storage_2<d_, Items_, Alloc_>::null_handle((std::numeric_limits<Index_type>::max)()/2);

} // namespace CGAL

#if defined(BOOST_GCC)
 _Pragma("GCC diagnostic pop")
#endif

#endif // CGAL_GENERALIZED_MAP_STORAGES_WITH_INDEX_H //
// EOF //
