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
#ifndef CGAL_GMAP_LINEAR_CELL_COMPLEX_STORAGES_H
#define CGAL_GMAP_LINEAR_CELL_COMPLEX_STORAGES_H 1

#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>
#include <CGAL/Dart.h>
#include <CGAL/Handle_hash_function.h>
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

  // Storage of darts with compact container, alpha with handles
  // Copy of Generalized_map_storage_1 and add new types related
  // to geometry (not possible to inherit because we use Self type
  // as template parameter of Dart_wrapper. If we inherit, Self is not
  // the correct type).
  template<unsigned int d_, unsigned int ambient_dim,
           class Traits_, class Items_, class Alloc_>
  class GMap_linear_cell_complex_storage_1
  {
  public:
    using Self=GMap_linear_cell_complex_storage_1<d_, ambient_dim, Traits_,
    Items_, Alloc_>;
    using Use_index=CGAL::Tag_false;
    using Concurrent_tag=typename internal::Get_concurrent_tag<Items_>::type;

    typedef typename Traits_::Point  Point;
    typedef typename Traits_::Vector Vector;
    typedef typename Traits_::FT     FT;

    typedef internal::Combinatorial_map_helper<Self>      Helper;

    typedef typename Items_::template Dart_wrapper<Self>  Dart_wrapper;

    typedef typename internal::template Get_dart_info<Dart_wrapper>::type
                                                           Dart_info;
    typedef typename internal::template Get_darts_with_id<Dart_wrapper>::type
                                                           Darts_with_id;
    typedef CGAL::Dart<d_, Self, Dart_info, Darts_with_id> Dart;

    typedef std::allocator_traits<Alloc_> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Dart> Dart_allocator;

    typedef typename internal::Container_type
                 <Concurrent_tag, Dart, Dart_allocator>::type Dart_container;

    typedef typename Dart_container::iterator              Dart_descriptor;
    typedef typename Dart_container::const_iterator        Dart_const_descriptor;
    typedef typename Dart_container::size_type             size_type;

    typedef std::nullptr_t Null_descriptor_type;
    inline static constexpr Null_descriptor_type null_descriptor=nullptr;

    using Type_for_compact_container=void*;

    typedef Items_ Items;
    typedef Alloc_ Alloc;
    template <typename T>
    struct Container_for_attributes :
      public internal::Container_type
                     <Concurrent_tag, T,
                      typename Allocator_traits::template rebind_alloc<T>>::type
    {};
    /// Typedef for attributes
    typedef typename internal::template Get_attributes_tuple<Dart_wrapper>::type
                                   Attributes;

    template<int i>
    struct Attribute_type: public Helper::template Attribute_type<i>
    {};
    template<int i>
    struct Attribute_descriptor: public Helper::template Attribute_descriptor<i>
    {};
    template<int i>
    struct Attribute_const_descriptor:
      public Helper::template Attribute_const_descriptor<i>
    {};
    template<int i>
    struct Attribute_range: public Helper::template Attribute_range<i>
    {};
    template<int i>
    struct Attribute_const_range:
      public Helper::template Attribute_const_range<i>
    {};

    typedef typename Attribute_type<0>::type Vertex_attribute;
    typedef typename Attribute_descriptor<0>::type Vertex_attribute_descriptor;
    typedef typename Attribute_const_descriptor<0>::type
    Vertex_attribute_const_descriptor;

    typedef typename Attribute_range<0>::type Vertex_attribute_range;
    typedef typename Attribute_const_range<0>::type
    Vertex_attribute_const_range;

    /// Deprecated types, keep for now for backward compatibility
    using Dart_handle=Dart_descriptor;
    using Dart_const_handle=Dart_const_descriptor;

    template<int i>
    using Attribute_handle=Attribute_descriptor<i>;
    template<int i>
    using Attribute_const_handle=Attribute_const_descriptor<i>;
    using Vertex_attribute_handle=Vertex_attribute_descriptor;
    using Vertex_attribute_const_handle=Vertex_attribute_const_descriptor;

    inline static constexpr Null_descriptor_type null_handle=null_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /// The dimension of the generalized map.
    static const unsigned int dimension = d_;

    typedef Handle_hash_function Hash_function;

    typedef Dart_container       Dart_range;
    typedef const Dart_container Dart_const_range;
    /// @return a Dart_range (range through all the darts of the map).
    Dart_range& darts()             { return mdarts;}
    Dart_const_range& darts() const { return mdarts; }

    // Init
    void init_storage()
    {
      null_dart_descriptor=nullptr;
      null_dart_handle=null_dart_descriptor;
    }

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

    size_type upper_bound_on_dart_ids() const
    { return 0; }

   /** Return if this dart is free for adimension.
     * @param dh a dart handle
     * @param i the dimension.
     * @return true iff dh is linked with itself for \em adimension.
     */
    template<unsigned int i>
    bool is_free(Dart_const_descriptor dh) const
    {
      CGAL_assertion( dh!=nullptr );
      CGAL_assertion(i <= dimension);
      return dh->mf[i]==dh;
    }
    bool is_free(Dart_const_descriptor dh, unsigned int i) const
    {
      CGAL_assertion( dh!=nullptr );
      CGAL_assertion(i <= dimension);
      return dh->mf[i]==dh;
    }
    bool is_perforated(Dart_const_descriptor /*dh*/) const
    { return false; }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_marks(amarks);
    }
    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(Dart_const_descriptor ADart) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_marks();
    }
    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_descriptor ADart, size_type amark, bool avalue) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_descriptor ADart, size_type amark) const
    {
      CGAL_assertion( ADart!=nullptr );
      ADart->flip_mark(amark);
    }

    // Access to alpha maps
    Dart_descriptor get_alpha(Dart_descriptor ADart, int B1)
    {
      CGAL_assertion(ADart!=nullptr && B1>=0 && B1<=(int)dimension);
      return ADart->mf[B1];
    }
    Dart_const_descriptor get_alpha(Dart_const_descriptor ADart, int B1) const
    {
      CGAL_assertion(ADart!=nullptr && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }
    template<int B1>
    Dart_descriptor get_alpha(Dart_descriptor ADart)
    {
      CGAL_assertion(ADart!=nullptr && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }
    template<int B1>
    Dart_const_descriptor get_alpha(Dart_const_descriptor ADart) const
    {
      CGAL_assertion(ADart!=nullptr && B1>=0 && B1<=(int)dimension);
      return  ADart->mf[B1];
    }

    // return a handle on the i-attribute
    template<unsigned int i>
    typename Attribute_descriptor<i>::type attribute(Dart_descriptor ADart)
    {
      static_assert(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (ADart->mattribute_descriptors);
    }
    template<unsigned int i>
    typename Attribute_const_descriptor<i>::type
    attribute(Dart_const_descriptor ADart) const
    {
      static_assert(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (ADart->mattribute_descriptors);
    }

    // Copy a given attribute
    template<unsigned int i>
    typename Attribute_descriptor<i>::type copy_attribute
    (typename Attribute_const_descriptor<i>::type ah)
    {
      static_assert(Helper::template Dimension_index<i>::value>=0,
                     "copy_attribute<i> called but i-attributes are disabled.");
      typename Attribute_descriptor<i>::type res=
        std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(*ah);
      this->template init_attribute_ref_counting<i>(res);
      return res;
    }

    // Test if a given attribute is valid
    template<unsigned int i>
    bool is_valid_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=nullptr );
      return ah->is_valid();
    }

    // accessors and modifiers to the attribute ref counting given its handle
    template<unsigned int i>
    std::size_t get_attribute_ref_counting
    (typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=nullptr );
      return ah->get_nb_refs();
    }
    template<unsigned int i>
    void init_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      ah->mrefcounting=0;
    }
    template<unsigned int i>
    void inc_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      ah->inc_nb_refs();
    }
    template<unsigned int i>
    void dec_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      ah->dec_nb_refs();
    }

    // get the attribute given its handle
    template<unsigned int i>
    typename Attribute_type<i>::type&
    get_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      return *ah;
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=nullptr );
      return *ah;
    }

    // Get the dart of the given attribute
    template<unsigned int i>
    Dart_descriptor dart_of_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      return ah->dart();
    }
    template<unsigned int i>
    Dart_const_descriptor
    dart_of_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=nullptr );
      return ah->dart();
    }

    // Set the dart of the given attribute
    template<unsigned int i>
    void set_dart_of_attribute(typename Attribute_descriptor<i>::type ah,
                               Dart_descriptor adart)
    {
      CGAL_assertion( ah!=nullptr );
      ah->set_dart(adart);
    }

    // Get the information associated with a given dart
    Dart_info& info(Dart_descriptor adart)
    { return adart->info(); }
    const Dart_info& info(Dart_const_descriptor adart) const
    { return adart->info(); }

    // Get the info of the given attribute
    template<unsigned int i>
    typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=nullptr );
      return ah->info();
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=nullptr );
      return ah->info();
    }

    // Get the info of the i-cell attribute associated with the given dart
    template<unsigned int i>
    typename Attribute_type<i>::type::Info & info(Dart_descriptor adart)
    {
      CGAL_assertion( adart!=nullptr );
      CGAL_assertion( attribute<i>(adart)!=nullptr );
      return info_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info(Dart_const_descriptor adart) const
    {
      CGAL_assertion( adart!=nullptr );
      CGAL_assertion( attribute<i>(adart)!=nullptr );
      return info_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the i-cell attribute associated with the given dart
    template<unsigned int i>
    Dart_descriptor dart(Dart_descriptor adart)
    {
      CGAL_assertion( adart!=nullptr );
      CGAL_assertion( attribute<i>(adart)!=nullptr );
      return dart_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    Dart_const_descriptor dart(Dart_const_descriptor adart) const
    {
      CGAL_assertion( adart!=nullptr );
      CGAL_assertion( attribute<i>(adart)!=nullptr );
      return dart_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the given 0-attribute
    Point & point_of_vertex_attribute(typename Attribute_descriptor<0>::type vh)
    {
      CGAL_assertion( vh!=nullptr );
      return get_attribute<0>(vh).point();
    }

    const Point & point_of_vertex_attribute
    (typename Attribute_const_descriptor<0>::type vh) const
    {
      CGAL_assertion( vh!=nullptr );
      return get_attribute<0>(vh).point();
    }

    // Debug function
    void display_dart(Dart_const_descriptor ADart) const
    { std::cout<<mdarts.index(ADart); }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_descriptor<i>::type ah) const
    { std::cout<< std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).index(ah); }

    template <unsigned int i>
    size_type upper_bound_on_attribute_ids() const
    { return 0; }

  protected:
    // Set the handle on the i th attribute
    template<unsigned int i>
    void basic_set_dart_attribute(Dart_descriptor dh,
                                  typename Attribute_descriptor<i>::type ah)
    {
      std::get<Helper::template Dimension_index<i>::value>
        (dh->mattribute_descriptors) = ah;
    }

    /** Link a dart with a given dart for a given dimension.
     * @param adart the dart to link.
     * @param adart2 the dart to link with.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_link_alpha(Dart_descriptor adart, Dart_descriptor adart2)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=nullptr && adart2!=nullptr);
      adart->mf[i] = adart2;
    }
    void dart_link_alpha(Dart_descriptor adart, Dart_descriptor adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=nullptr && adart2!=nullptr);
      adart->mf[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_alpha(Dart_descriptor adart)
    {
      CGAL_assertion(adart!=nullptr && i <= dimension);
      adart->mf[i] = adart;
    }
    void dart_unlink_alpha(Dart_descriptor adart, unsigned int i)
    {
      CGAL_assertion(adart!=nullptr && i <= dimension);
      adart->mf[i] = adart;
    }

  public:
    Dart_descriptor null_dart_descriptor; // To be compatible with combinatorial map
    Dart_descriptor null_dart_handle; // Deprecated: kept for backward compatibility

  protected:
    /// Dart container.
    Dart_container mdarts;

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

} // namespace CGAL

#if defined(BOOST_GCC)
 _Pragma("GCC diagnostic pop")
#endif
#endif // CGAL_GMAP_LINEAR_CELL_COMPLEX_STORAGES_H //
// EOF //
