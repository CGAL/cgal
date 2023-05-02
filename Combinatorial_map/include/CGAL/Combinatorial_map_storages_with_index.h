// Copyright (c) 2013 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_STORAGES_WITH_INDEX_H
#define CGAL_COMBINATORIAL_MAP_STORAGES_WITH_INDEX_H 1

#include <CGAL/Compact_container_with_index.h>
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

  // Storage with combinatorial maps using index
  template<unsigned int d_, class Items_, class Alloc_>
  class Combinatorial_map_storage_with_index
  {
  public:
    using Self=Combinatorial_map_storage_with_index<d_, Items_, Alloc_>;
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

    typedef Compact_container_with_index<Dart,Dart_allocator,
    Multiply_by_two_policy_for_cc_with_size<64>, Index_type>
    Dart_container;

    typedef typename Dart_container::Index Dart_index;

    typedef Dart_index Dart_descriptor;
    typedef Dart_index Dart_const_descriptor;
    typedef typename Dart_container::size_type size_type;

    typedef Dart_index Null_descriptor_type;
    static const Index_type null_descriptor=Dart_container::null_descriptor;

    using Type_for_compact_container=Index_type;

    typedef Items_ Items;
    typedef Alloc_ Alloc;
    template <typename T>
    struct Container_for_attributes : public
        Compact_container_with_index<T,
        typename Allocator_traits::template rebind_alloc<T>,
        Multiply_by_two_policy_for_cc_with_size<64>, size_type >
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

    /// Deprecated types, keep for now for backward compatibility
    using Dart_handle=Dart_descriptor;
    using Dart_const_handle=Dart_const_descriptor;

    template<int i>
    using Attribute_handle=Attribute_descriptor<i>;
    template<int i>
    using Attribute_const_handle=Attribute_const_descriptor<i>;

    static const Index_type null_handle=null_descriptor;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d_;

    typedef internal::Index_hash_function Hash_function;

    //**************************************************************************
    // Dart_range
    struct Dart_range
    {
      typedef typename Dart_container::iterator iterator;
      typedef typename Dart_container::const_iterator const_iterator;
      Dart_range(Self &amap) : mmap(amap)
      {}
      iterator begin()
      {
        CGAL_assertion(!mmap.mdarts.empty()); // The container mdarts is not empty
                                              // it has at least the null_dart
        iterator res=mmap.mdarts.begin(); ++res; return res;
      }
      iterator end() { return mmap.mdarts.end(); }
      const_iterator begin() const
      { const_iterator res=mmap.mdarts.begin(); ++res; return res; }
      const_iterator end() const { return mmap.mdarts.end(); }
      size_type size() const
      { return mmap.mdarts.size()-1; }
      bool empty() const
      { return size()==0; } // mmap.is_empty(); }
      size_type index(Dart_index idx) const
      { return idx; }
      size_type index(const_iterator cit) const
      { return cit; }
      bool is_used(size_type i) const
      { return mmap.mdarts.is_used(i); }
      bool owns(size_type i) const
      { return mmap.mdarts.owns(i); }
    private:
      Self & mmap;
    };
    typedef const Dart_range Dart_const_range;

    /// @return a Dart_range (range through all the darts of the map).
    Dart_range& darts()             { return mdarts_range;}
    Dart_const_range& darts() const { return mdarts_range; }
    //**************************************************************************

    Combinatorial_map_storage_with_index() : mdarts_range(*this)
    {}

    void init_storage()
    {
      // Allocate a dart for null_dart_descriptor
      assert(mdarts.empty()); // the compact container is empty
      Dart_index local_null_dart_descriptor = mdarts.emplace();
      if(local_null_dart_descriptor!=0)
      {
        std::cerr<<"[ERROR] fatal in Combinatorial_map_storage_with_index::init_storage"
                 <<std::endl;
        CGAL_assertion(false);
      }
    }

    void clear_storage()
    { mdarts.clear(); }

    /** Test if the map is empty.
     *  @return true iff the map is empty.
     */
    bool is_empty() const
    { return darts().empty(); }

    /// @return the number of darts.
    size_type number_of_darts() const
    { return darts().size(); }

    size_type upper_bound_on_dart_ids() const
    { return mdarts.upper_bound(); }

    /** Return if this dart is free for adimension.
     * @param dh a dart handle
     * @param i the dimension.
     * @return true iff dh is linked with null_dart_descriptor for \em adimension.
     */
    template<unsigned int i>
    bool is_free(Dart_const_descriptor dh) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mf[i]==null_dart_descriptor;
    }
    bool is_free(Dart_const_descriptor dh, unsigned int i) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mf[i]==null_dart_descriptor;
    }
    bool is_perforated(Dart_const_descriptor /*dh*/) const
    { return false; }

    /// Set simultaneously all the marks of this dart to a given value.
    void set_dart_marks(Dart_const_descriptor ADart,
                        const std::bitset<NB_MARKS>& amarks) const
    {
      mdarts[ADart].set_marks(amarks);
    }
    /// Return all the marks of a dart.
    std::bitset<NB_MARKS> get_dart_marks(Dart_const_descriptor ADart) const
    {
      return mdarts[ADart].get_marks();
    }
    /// Return the mark value of dart a given mark number.
    bool get_dart_mark(Dart_const_descriptor ADart, size_type amark) const
    {
      return mdarts[ADart].get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_descriptor ADart, size_type amark, bool avalue) const
    {
      mdarts[ADart].set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_descriptor ADart, size_type amark) const
    {
      mdarts[ADart].flip_mark(amark);
    }

    // Access to beta maps
    Dart_descriptor get_beta(Dart_descriptor ADart, int B1)
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return mdarts[ADart].mf[B1];
    }
    Dart_const_descriptor get_beta(Dart_const_descriptor ADart, int B1) const
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }
    template<int B1>
    Dart_descriptor get_beta(Dart_descriptor ADart)
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }
    template<int B1>
    Dart_const_descriptor get_beta(Dart_const_descriptor ADart) const
    {
      CGAL_assertion(B1>=0 && B1<=(int)dimension);
      return  mdarts[ADart].mf[B1];
    }

    // return a handle on the i-attribute
    template<unsigned int i>
    typename Attribute_descriptor<i>::type attribute(Dart_descriptor ADart)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_descriptors);
    }
    template<unsigned int i>
    typename Attribute_const_descriptor<i>::type
    attribute(Dart_const_descriptor ADart) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return std::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_descriptors);
    }

    // Copy a given attribute
    template<unsigned int i>
    typename Attribute_descriptor<i>::type copy_attribute
    (typename Attribute_const_descriptor<i>::type ah)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "copy_attribute<i> called but i-attributes are disabled.");
      // We need to do a reserve before the emplace in order to avoid a bug of
      // invalid reference when the container is reallocated.
      std::get<Helper::template Dimension_index<i>::value>
          (mattribute_containers).reserve
          (std::get<Helper::template Dimension_index<i>::value>
           (mattribute_containers).size()+1);

      typename Attribute_descriptor<i>::type res=
        std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(get_attribute<i>(ah));
      this->template init_attribute_ref_counting<i>(res);
      return res;
    }

    // Test if a given attribute is valid
    template<unsigned int i>
    bool is_valid_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      return get_attribute<i>(ah).is_valid();
    }

    // accessors and modifiers to the attribute ref counting given its handle
    template<unsigned int i>
    std::size_t get_attribute_ref_counting
    (typename Attribute_const_descriptor<i>::type ah) const
    {
      return get_attribute<i>(ah).get_nb_refs();
    }
    template<unsigned int i>
    void init_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      get_attribute<i>(ah).mrefcounting=0;
    }
    template<unsigned int i>
    void inc_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      get_attribute<i>(ah).inc_nb_refs();
    }
    template<unsigned int i>
    void dec_attribute_ref_counting(typename Attribute_descriptor<i>::type ah)
    {
      get_attribute<i>(ah).dec_nb_refs();
    }

    // get the attribute given its index
    template<unsigned int i>
    typename Attribute_type<i>::type&
    get_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=null_descriptor );
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=null_descriptor );
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }

    // Get the dart of the given attribute
    template<unsigned int i>
    Dart_descriptor dart_of_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=null_descriptor );
      return get_attribute<i>(ah).dart();
    }
    template<unsigned int i>
    Dart_const_descriptor
    dart_of_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=null_descriptor );
      return get_attribute<i>(ah).dart();
    }

    // Set the dart of the given attribute
    template<unsigned int i>
    void set_dart_of_attribute(typename Attribute_descriptor<i>::type ah,
                               Dart_descriptor adart)
    {
      CGAL_assertion( ah!=null_descriptor );
      get_attribute<i>(ah).set_dart(adart);
    }

    // Get the information associated with a given dart
    Dart_info& info(Dart_descriptor adart)
    { return mdarts[adart].info(); }
    const Dart_info& info(Dart_const_descriptor adart) const
    { return mdarts[adart].info(); }

    // Get the info of the given attribute
    template<unsigned int i>
    typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_descriptor<i>::type ah)
    {
      CGAL_assertion( ah!=null_descriptor );
      return get_attribute<i>(ah).info();
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_const_descriptor<i>::type ah) const
    {
      CGAL_assertion( ah!=null_descriptor );
      return get_attribute<i>(ah).info();
    }

    // Get the info of the i-cell attribute associated with the given dart
    template<unsigned int i>
    typename Attribute_type<i>::type::Info & info(Dart_descriptor adart)
    {
      CGAL_assertion( adart!=null_descriptor );
      CGAL_assertion( this->template attribute<i>(adart)!=null_descriptor );
      return info_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info(Dart_const_descriptor adart) const
    {
      CGAL_assertion( adart!=null_descriptor );
      CGAL_assertion( attribute<i>(adart)!=null_descriptor );
      return info_of_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the i-cell attribute associated with the given dart
    template<unsigned int i>
    Dart_descriptor dart(Dart_descriptor adart)
    {
      CGAL_assertion( adart!=null_descriptor );
      CGAL_assertion( attribute<i>(adart)!=null_descriptor );
      return dart_of_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    Dart_const_descriptor dart(Dart_const_descriptor adart) const
    {
      CGAL_assertion( adart!=null_descriptor );
      CGAL_assertion( attribute<i>(adart)!=null_descriptor );
      return dart_of_attribute<i>(attribute<i>(adart));
    }

    // Debug function
    void display_dart(Dart_const_descriptor ADart) const
    { std::cout<<ADart; }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_descriptor<i>::type ah) const
    { std::cout<<ah; }

    template <unsigned int i>
    size_type upper_bound_on_attribute_ids() const
    {
      return std::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).upper_bound();
    }

    protected:
    // Set the handle on the i th attribute
    template<unsigned int i>
    void basic_set_dart_attribute(Dart_descriptor dh,
                                  typename Attribute_descriptor<i>::type ah)
    {
      std::get<Helper::template Dimension_index<i>::value>
        (mdarts[dh].mattribute_descriptors) = ah;
    }

    /** Link a dart with a given dart for a given dimension.
     * @param adart the dart to link.
     * @param adart2 the dart to link with.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_link_beta(Dart_descriptor adart, Dart_descriptor adart2)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=null_dart_descriptor);
      mdarts[adart].mf[i] = adart2;
    }
    void dart_link_beta(Dart_descriptor adart, Dart_descriptor adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=null_dart_descriptor);
      mdarts[adart].mf[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_beta(Dart_descriptor adart)
    {
      CGAL_assertion(i <= dimension);
      mdarts[adart].mf[i] = null_dart_descriptor;
    }
    void dart_unlink_beta(Dart_descriptor adart, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      mdarts[adart].mf[i] = null_dart_descriptor;
    }

  public:
    /// Void dart. A dart d is i-free if beta_i(d)=null_dart_descriptor.
    Index_type null_dart_descriptor=0;

    // Deprecated: kept for backward compatibility
    Index_type null_dart_handle=null_dart_descriptor;

  protected:
    /// Dart container.
    Dart_container mdarts;
    Dart_range mdarts_range;

    /// Container for the null_dart_descriptor: unused; to be compatible with handle version
    Dart_container mnull_dart_container;

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

} // namespace CGAL

#if defined(BOOST_GCC)
 _Pragma("GCC diagnostic pop")
#endif

#endif // CGAL_COMBINATORIAL_MAP_STORAGES_WITH_INDEX_H //
// EOF //
