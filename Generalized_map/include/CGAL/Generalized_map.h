// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_H
#define CGAL_GENERALIZED_MAP_H 1

#include <CGAL/internal/Combinatorial_map_utility.h>
#include <CGAL/internal/Generalized_map_group_functors.h>
#include <CGAL/internal/Combinatorial_map_copy_functors.h>
#include <CGAL/internal/Generalized_map_sewable.h>

#include <CGAL/Generalized_map_storages.h>
#include <CGAL/Combinatorial_map_functors.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Generalized_map_operations.h>
#include <CGAL/Generic_map_min_items.h>
#include <CGAL/GMap_dart_const_iterators.h>
#include <CGAL/GMap_cell_const_iterators.h>
#include <CGAL/Generalized_map_save_load.h>

#include <CGAL/Unique_hash_map.h>
#include <bitset>
#include <vector>
#include <deque>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/config.h>

#if defined( __INTEL_COMPILER )
// Workarounf for warning in function basic_link_beta_0
#pragma warning disable 1017
#endif

#include <boost/config.hpp>
#if  (BOOST_GCC >= 40900)
_Pragma("GCC diagnostic push")
_Pragma("GCC diagnostic ignored \"-Warray-bounds\"")
#endif

namespace CGAL {

  /** @file Generalized_map.h
   * Definition of generic dD Generalized map.
   */

  struct Generalized_map_tag {};

  /** Generic definition of generalized map in dD.
   * The Generalized_map class describes an dD generalized map. It allows
   * mainly to create darts, to use marks onto these darts, to get and set
   * the alpha links, and to manage enabled attributes.
   */
  template < unsigned int d_, class Refs,
             class Items_=Generic_map_min_items,
             class Alloc_=CGAL_ALLOCATOR(int),
             class Storage_= Generalized_map_storage_1<d_, Items_, Alloc_> >
  class Generalized_map_base: public Storage_
  {
    template<typename Gmap,unsigned int i,typename Enabled>
    friend struct CGAL::internal::Call_merge_functor;
    template<typename Gmap,unsigned int i,typename Enabled>
    friend struct CGAL::internal::Call_split_functor;
    template<class Map, unsigned int i, unsigned int nmi>
    friend struct Remove_cell_functor;
    template<class Map, unsigned int i>
    friend struct Contract_cell_functor;
    template<typename Gmap>
    friend struct internal::Init_attribute_functor;
    template<typename CMap>
    friend struct Swap_attributes_functor;

  public:
    template < unsigned int A, class B, class I, class D, class S >
    friend class Generalized_map_base;

    typedef Generalized_map_tag Combinatorial_data_structure;

    /// Types definition
    typedef Storage_                                                  Storage;
    typedef Storage                                                   Base;
    typedef Generalized_map_base<d_, Refs, Items_, Alloc_, Storage_ > Self;

    typedef typename Base::Dart Dart;
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;
    typedef typename Base::Dart_container Dart_container;
    typedef typename Base::size_type size_type;
    typedef typename Base::Helper Helper;
    typedef typename Base::Attributes Attributes;
    typedef typename Base::Items Items;
    typedef typename Base::Alloc Alloc;
    typedef typename Base::Use_index Use_index;

    static const size_type NB_MARKS = Base::NB_MARKS;
    static const size_type INVALID_MARK = NB_MARKS;

    static const unsigned int dimension = Base::dimension;

    typedef typename Base::Null_handle_type Null_handle_type;
    using Base::null_handle;

    using Base::mdarts;
    using Base::get_alpha;
    using Base::is_free;
    using Base::set_dart_mark;
    using Base::get_dart_mark;
    using Base::flip_dart_mark;
    using Base::set_dart_marks;
    using Base::get_dart_marks;
    using Base::dart_link_alpha;
    using Base::dart_unlink_alpha;
    using Base::attribute;
    using Base::mattribute_containers;
    using Base::dart_of_attribute;
    using Base::set_dart_of_attribute;
    using Base::info_of_attribute;
    using Base::info;
    using Base::dart;

    /// Typedef for Dart_range, a range through all the darts of the map.
    typedef Dart_container       Dart_range;
    typedef const Dart_container Dart_const_range;

    /// Typedef for attributes
    template<int i>
    struct Attribute_type: public Base::template Attribute_type<i>
    {};
    template<int i>
    struct Attribute_handle: public Base::template Attribute_handle<i>
    {};
    template<int i>
    struct Attribute_const_handle:
      public Base::template Attribute_const_handle<i>
    {};
    template<int i>
    struct Attribute_range: public Base::template Attribute_range<i>
    {};
    template<int i>
    struct Attribute_const_range:
      public Base::template Attribute_const_range<i>
    {};

    class Exception_no_more_available_mark {};

  public:
    /** Default Generalized_map constructor.
     * The map is empty.
     */
    Generalized_map_base()
    {
      CGAL_static_assertion_msg(Helper::nb_attribs<=dimension+1,
                  "Too many attributes in the tuple Attributes_enabled");
      this->init_storage();

      this->mnb_used_marks = 0;
      this->mmask_marks.reset();

      for ( size_type i = 0; i < NB_MARKS; ++i)
      {
        this->mfree_marks_stack[i]        = i;
        this->mindex_marks[i]             = i;
        this->mnb_marked_darts[i]         = 0;
        this->mnb_times_reserved_marks[i] = 0;
      }

      this->automatic_attributes_management = true;

      CGAL_assertion(number_of_darts()==0);
    }

    /** Copy the given generalized map into *this.
     *  Note that both Gmap can have different dimensions and/or non void attributes.
     *  @param amap the generalized map to copy.
     *  @post *this is valid.
     */
    template <typename GMap2, typename Converters, typename DartInfoConverter,
              typename PointConverter>
    void copy(const GMap2& amap, const Converters& converters,
              const DartInfoConverter& dartinfoconverter,
              const PointConverter& pointconverter)
    {
      this->clear();

      this->mnb_used_marks = amap.mnb_used_marks;
      this->mmask_marks    = amap.mmask_marks;
      this->automatic_attributes_management =
          amap.automatic_attributes_management;

      for (size_type i = 0; i < NB_MARKS; ++i)
      {
        this->mfree_marks_stack[i]        = amap.mfree_marks_stack[i];
        this->mused_marks_stack[i]        = amap.mused_marks_stack[i];
        this->mindex_marks[i]             = amap.mindex_marks[i];
        this->mnb_marked_darts[i]         = amap.mnb_marked_darts[i];
        this->mnb_times_reserved_marks[i] = amap.mnb_times_reserved_marks[i];
      }

      // Create an mapping between darts of the two maps (originals->copies).
      // TODO: replace the std::map by a boost::unordered_map
      // (here we cannot use CGAL::Unique_hash_map because it does not provide
      // iterators...
      std::map<typename GMap2::Dart_const_handle, Dart_handle> dartmap;

      for (typename GMap2::Dart_const_range::const_iterator
             it=amap.darts().begin(), itend=amap.darts().end();
           it!=itend; ++it)
      {
        dartmap[it]=mdarts.emplace();
        init_dart(dartmap[it], amap.get_marks(it));
        internal::Copy_dart_info_functor<GMap2, Refs, DartInfoConverter>::
          run(amap, static_cast<Refs&>(*this), it, dartmap[it],
              dartinfoconverter);
      }

      unsigned int min_dim=(dimension<amap.dimension?dimension:amap.dimension);

      typename std::map<typename GMap2::Dart_const_handle,Dart_handle>
        ::iterator dartmap_iter, dartmap_iter_end=dartmap.end();
      for (dartmap_iter=dartmap.begin(); dartmap_iter!=dartmap_iter_end;
           ++dartmap_iter)
      {
        for (unsigned int i=0; i<=min_dim; i++)
        {
          if (!amap.is_free(dartmap_iter->first,i) &&
              (dartmap_iter->first)<(amap.alpha(dartmap_iter->first,i)))
          {
            basic_link_alpha(dartmap_iter->second,
                            dartmap[amap.alpha(dartmap_iter->first,i)], i);
          }
        }
      }

      /** Copy attributes */
      for (dartmap_iter=dartmap.begin(); dartmap_iter!=dartmap_iter_end;
           ++dartmap_iter)
      {
        Helper::template Foreach_enabled_attributes
          < internal::Copy_attributes_functor <GMap2, Refs, Converters,
            PointConverter> >::
          run(amap, static_cast<Refs&>(*this),
              dartmap_iter->first, dartmap_iter->second,
              converters, pointconverter);
      }

      CGAL_assertion (is_valid () == 1);
    }

    template <typename GMap2>
    void copy(const GMap2& amap)
    {
      CGAL::cpp11::tuple<> converters;
      Default_converter_dart_info<GMap2, Refs> dartinfoconverter;
      Default_converter_cmap_0attributes_with_point<GMap2, Refs> pointconverter;
      copy(amap, converters, dartinfoconverter, pointconverter);
    }

    template <typename GMap2, typename Converters>
    void copy(const GMap2& amap, const Converters& converters)
    {
      Default_converter_cmap_0attributes_with_point<GMap2, Refs> pointconverter;
      Default_converter_dart_info<GMap2, Refs> dartinfoconverter;
      copy(amap, converters, dartinfoconverter, pointconverter);
    }

    template <typename GMap2, typename Converters, typename DartInfoConverter>
    void copy(const GMap2& amap, const Converters& converters,
              const DartInfoConverter& dartinfoconverter)
    {
      Default_converter_cmap_0attributes_with_point<GMap2, Refs> pointconverter;
      copy(amap, converters, dartinfoconverter, pointconverter);
    }

    // Copy constructor from a map having exactly the same type.
    Generalized_map_base (const Self & amap)
    { copy<Self>(amap); }

    // "Copy constructor" from a map having different type.
    template <typename GMap2>
    Generalized_map_base(const GMap2& amap)
    { copy(amap); }

    // "Copy constructor" from a map having different type.
    template <typename GMap2, typename Converters>
    Generalized_map_base(const GMap2& amap, Converters& converters)
    { copy(amap, converters); }

    // "Copy constructor" from a map having different type.
    template <typename GMap2, typename Converters, typename DartInfoConverter>
    Generalized_map_base(const GMap2& amap, Converters& converters,
                         const DartInfoConverter& dartinfoconverter)
    { copy(amap, converters, dartinfoconverter); }

    // "Copy constructor" from a map having different type.
    template <typename GMap2, typename Converters, typename DartInfoConverter,
              typename PointConverter>
    Generalized_map_base(const GMap2& amap, Converters& converters,
                         const DartInfoConverter& dartinfoconverter,
                         const PointConverter& pointconverter)
    { copy(amap, converters, dartinfoconverter, pointconverter); }

    /** Affectation operation. Copies one map to the other.
     * @param amap a generalized map.
     * @return A copy of that generalized map.
     */
    Self & operator= (const Self & amap)
    {
      if (this!=&amap)
      {
        Self tmp(amap);
        this->swap(tmp);
      }
      return *this;
    }

    /** Swap this generalized map with amap, a second generalized map.
     * Note that the two maps have exactly the same type.
     * @param amap a generalized map.
     */
    void swap(Self & amap)
    {
      if (this!=&amap)
      {
        mdarts.swap(amap.mdarts);

        Helper::template Foreach_enabled_attributes
          < internal::Swap_attributes_functor <Self> >::run(*this, amap);

        std::swap_ranges(mnb_times_reserved_marks,
                         mnb_times_reserved_marks+NB_MARKS,
                         amap.mnb_times_reserved_marks);
        std::swap(mmask_marks,amap.mmask_marks);
        std::swap(mnb_used_marks, amap.mnb_used_marks);
        std::swap_ranges(mindex_marks,mindex_marks+NB_MARKS,
                         amap.mindex_marks);
        std::swap_ranges(mfree_marks_stack, mfree_marks_stack+NB_MARKS,
                         amap.mfree_marks_stack);
        std::swap_ranges(mused_marks_stack,mused_marks_stack+NB_MARKS,
                         amap.mused_marks_stack);
        std::swap_ranges(mnb_marked_darts,mnb_marked_darts+NB_MARKS,
                         amap.mnb_marked_darts);

        std::swap(automatic_attributes_management,
                  amap.automatic_attributes_management);
      }
    }

    /** Clear the generalized map. Remove all darts and all attributes.
     *  Note that reserved marks are not free.
     */
    void clear()
    {
      mdarts.clear();
      for ( size_type i = 0; i < NB_MARKS; ++i)
        this->mnb_marked_darts[i]  = 0;

      internal::Clear_all::run(mattribute_containers);
      this->init_storage();
    }

    /** Test if the map is empty.
     *  @return true iff the map is empty.
     */
    bool is_empty() const
    { return mdarts.empty(); }

    friend std::ostream& operator<< (std::ostream& os, const Self& amap)
    {
      save_generalized_map(amap, os);
      return os;
    }

    friend std::ifstream& operator>> (std::ifstream& is, Self& amap)
    {
      load_generalized_map(is, amap);
      return is;
    }

    /** Create a new dart and add it to the map.
     * The marks of the darts are initialised with mmask_marks, i.e. the dart
     * is unmarked for all the marks.
     * @return a Dart_handle on the new dart.
     */
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template < typename... Args >
    Dart_handle create_dart(const Args&... args)
    {
      Dart_handle res=mdarts.emplace(args...);
      init_dart(res);
      return res;
    }
#else
    Dart_handle create_dart()
    {
      Dart_handle res=mdarts.emplace();
      init_dart(res);
      return res;
    }
    template < typename T1 >
    Dart_handle create_dart(const T1 &t1)
    {
      Dart_handle res=mdarts.emplace(t1);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2)
    {
      Dart_handle res=mdarts.emplace(t1, t2);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4, typename T5 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4, const T5 &t5)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4, t5);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4, const T5 &t5, const T6 &t6)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4, t5, t6);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6, typename T7 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4, const T5 &t5, const T6 &t6,
                            const T7 &t7)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4, t5, t6, t7);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6, typename T7, typename T8 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4, const T5 &t5, const T6 &t6,
                            const T7 &t7, const T8 &t8)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4, t5, t6, t7, t8);
      init_dart(res);
      return res;
    }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6, typename T7, typename T8, typename T9 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3,
                            const T4 &t4, const T5 &t5, const T6 &t6,
                            const T7 &t7, const T8 &t8, const T9 &t9)
    {
      Dart_handle res=mdarts.emplace(t1, t2, t3, t4, t5, t6, t7, t8, t9);
      init_dart(res);
      return res;
    }
#endif

    /** Erase a dart from the list of darts.
     * @param adart the dart to erase.
     */
    void erase_dart(Dart_handle adart)
    {
      // 1) We update the number of marked darts.
      for ( size_type i = 0; i < mnb_used_marks; ++i)
      {
        if (is_marked(adart, mused_marks_stack[i]))
          --mnb_marked_darts[mused_marks_stack[i]];
      }

      // 2) We update the attribute_ref_counting.
      Helper::template Foreach_enabled_attributes
        <internal::Decrease_attribute_functor<Self> >::run(*this,adart);

      // 3) We erase the dart.
      mdarts.erase(adart);
    }

    /// @return true if dh points to a used dart (i.e. valid).
    bool is_dart_used(Dart_const_handle dh) const
    { return mdarts.is_used(dh); }

    /// @return a Dart_range (range through all the darts of the map).
    Dart_range& darts()             { return mdarts;}
    Dart_const_range& darts() const { return mdarts; }

    /** Get the first dart of this map.
     * @return the first dart.
     */
    Dart_handle first_dart()
    {
      if (darts().begin() == darts().end()) return null_handle;
      return mdarts.begin();
    }
    Dart_const_handle first_dart() const
    {
      if (darts().begin() == darts().end()) return null_handle;
      return mdarts.begin();
    }

    /// @return the Dart_handle corresponding to the given dart.
    Dart_handle dart_handle(Dart& adart)
    { return mdarts.iterator_to(adart); }
    Dart_const_handle dart_handle(const Dart& adart) const
    { return mdarts.iterator_to(adart); }

    /** Return the highest dimension for which dh is not free.
     * @param dh a dart handle
     * @return the dimension d such that dh is not d-free but k-free for
     *         all k>d. -1 if the dart is free for all d in {0..n}
     */
    int highest_nonfree_dimension(Dart_const_handle dh) const
    {
      for (int i=(int)dimension; i>=0; --i)
      { if ( !is_free(dh, i) ) return i; }
      return -1;
    }

    /** Return a dart belonging to the same edge and to the second vertex
     * of the current edge (NULL if such a dart does not exist).
     * @return An handle to a dart belonging to the other extremity.
     */
    Dart_handle other_extremity(Dart_handle dh)
    {
      if (!is_free(dh, 0)) return alpha(dh, 0);
      return null_handle;
    }
    Dart_const_handle other_extremity(Dart_const_handle dh) const
    {
      if (!is_free(dh, 0)) return alpha(dh, 0);
      return null_handle;
    }

    // Set the handle on the i th attribute
    template<unsigned int i>
    void set_dart_attribute(Dart_handle dh,
                            typename Attribute_handle<i>::type ah)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "set_dart_attribute<i> called but i-attributes are disabled.");

      if ( this->template attribute<i>(dh)==ah ) return;

      if ( this->template attribute<i>(dh)!=null_handle )
      {
        this->template dec_attribute_ref_counting<i>(this->template attribute<i>(dh));
        if ( this->are_attributes_automatically_managed() &&
             this->template get_attribute_ref_counting<i>
             (this->template attribute<i>(dh))==0 )
          this->template erase_attribute<i>(this->template attribute<i>(dh));
      }

      this->template basic_set_dart_attribute<i>(dh, ah);

      if ( ah!=null_handle )
      {
        this->template set_dart_of_attribute<i>(ah, dh);
        this->template inc_attribute_ref_counting<i>(ah);
      }
    }

  protected:
    /// Marks can be modified even for const handle; otherwise it is not
    /// possible to iterate through const generalized maps.

    // Initialize a given dart: all alpha to null_dart_handle and all
    // attributes to null, all marks unmarked.
    void init_dart(Dart_handle adart)
    {
      set_dart_marks(adart, mmask_marks);

      for (unsigned int i = 0; i <= dimension; ++i)
        dart_unlink_alpha(adart, i);

      Helper::template Foreach_enabled_attributes
          <internal::Init_attribute_functor<Self> >::run(*this, adart);
    }
    // Initialize a given dart: all alpha to this and all
    // attributes to null, marks are given.
    void init_dart(Dart_handle adart,
                   const std::bitset<NB_MARKS>& amarks)
    {
      set_marks(adart, amarks);

      for (unsigned int i = 0; i <= dimension; ++i)
        dart_unlink_alpha(adart, i);

      Helper::template Foreach_enabled_attributes
          <internal::Init_attribute_functor<Self> >::run(*this, adart);
    }

  public:

    /// @return the alphas of ADart (alpha are used in the same order than
    ///         they are given as parameters)

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template<typename ...Alphas>
    Dart_handle alpha(Dart_handle ADart, Alphas... alphas)
    { return CGAL::internal::Alpha_functor<Self, Dart_handle, Alphas...>::
        run(*this, ADart, alphas...); }
    template<typename ...Alphas>
    Dart_const_handle alpha(Dart_const_handle ADart, Alphas... alphas) const
    { return CGAL::internal::Alpha_functor<const Self, Dart_const_handle, Alphas...>::
        run(*this, ADart, alphas...); }
    template<int... Alphas>
    Dart_handle alpha(Dart_handle ADart)
    { return CGAL::internal::Alpha_functor_static<Self, Dart_handle, Alphas...>::
        run(*this, ADart); }
    template<int... Alphas>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return CGAL::internal::Alpha_functor_static<const Self, Dart_const_handle, Alphas...>::
        run(*this, ADart); }
#else
    Dart_handle alpha(Dart_handle ADart, int B1)
    { return this->get_alpha(ADart, B1); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2)
    { return alpha(alpha(ADart, B1), B2); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3)
    { return alpha(alpha(ADart, B1), B2, B3); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4)
    { return alpha(alpha(ADart, B1), B2, B3, B4); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4, int B5)
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6)
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7)
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7, int B8)
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7, B8); }
    Dart_handle alpha(Dart_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7, int B8, int B9)
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7, B8, B9); }

    template<int B1>
    Dart_handle alpha(Dart_handle ADart)
    { return this->template get_alpha<B1>(ADart); }
    template<int B1, int B2>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4, B5>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4, B5, B6>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
             int B7>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4, B5, B6, B7>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
             int B7, int B8>
    Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4, B5, B6, B7, B8>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
             int B7, int B8, int B9>
     Dart_handle alpha(Dart_handle ADart)
    { return alpha<B2, B3, B4, B5, B6, B7, B8, B9>(alpha<B1>(ADart)); }

    Dart_const_handle alpha(Dart_const_handle ADart, int B1) const
    { return this->get_alpha(ADart, B1); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2) const
    { return alpha(alpha(ADart, B1), B2); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3) const
    { return alpha(alpha(ADart, B1), B2, B3); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4) const
    { return alpha(alpha(ADart, B1), B2, B3, B4); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4, int B5) const
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6) const
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7) const
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7, int B8) const
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7, B8); }
    Dart_const_handle alpha(Dart_const_handle ADart, int B1, int B2, int B3,
                           int B4, int B5, int B6, int B7, int B8, int B9) const
    { return alpha(alpha(ADart, B1), B2, B3, B4, B5, B6, B7, B8, B9); }

    template<int B1>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return this->template get_alpha<B1>(ADart); }
    template<int B1, int B2>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4, B5>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4, B5, B6>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
              int B7>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4, B5, B6, B7>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
              int B7, int B8>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4, B5, B6, B7, B8>(alpha<B1>(ADart)); }
    template<int B1, int B2, int B3, int B4, int B5, int B6,
              int B7, int B8, int B9>
    Dart_const_handle alpha(Dart_const_handle ADart) const
    { return alpha<B2, B3, B4, B5, B6, B7, B8, B9>(alpha<B1>(ADart)); }
#endif

    // Generic function to iterate on CMap or GMap in a generic way
    bool is_previous_exist(Dart_const_handle ADart) const
    { return !this->template is_free<1>(ADart) &&
        !this->template is_free<0>(alpha<1>(ADart)); }
    bool is_next_exist(Dart_const_handle ADart) const
    { return !this->template is_free<0>(ADart) &&
        !this->template is_free<1>(this->template alpha<0>(ADart)); }
    template<unsigned int dim>
    bool is_opposite_exist(Dart_const_handle ADart) const
    { return !this->template is_free<dim>(ADart) &&
        !this->template is_free<0>(this->template alpha<dim>(ADart)); }

    Dart_handle previous(Dart_handle ADart)
    { return this->template alpha<1, 0>(ADart); }
    Dart_const_handle previous(Dart_const_handle ADart) const
    { return this->template alpha<1, 0>(ADart); }

    Dart_handle next(Dart_handle ADart)
    { return this->template alpha<0, 1>(ADart); }
    Dart_const_handle next(Dart_const_handle ADart) const
    { return this->template alpha<0, 1>(ADart); }

    template<unsigned int dim>
    Dart_handle opposite(Dart_handle ADart)
    { return this->template alpha<0, dim>(ADart); }
    template<unsigned int dim>
    Dart_const_handle opposite(Dart_const_handle ADart) const
    { return this->template alpha<0, dim>(ADart); }

    void set_next(Dart_handle dh1, Dart_handle dh2)
    {
      assert(!this->template is_free<0>(dh1));
      this->template link_alpha<1>(this->template alpha<0>(dh1), dh2); }

    template<unsigned int dim>
    void set_opposite(Dart_handle dh1, Dart_handle dh2)
    {
      assert(!this->template is_free<0>(dh1));
      assert(!this->template is_free<0>(dh2));
      this->template link_alpha<dim>(this->template alpha<0>(dh1), dh2);
      this->template link_alpha<dim>(dh1, this->template alpha<0>(dh2));
    }

    Dart_handle other_orientation(Dart_handle ADart)
    {
      CGAL_assertion(!this->template is_free<0>(ADart));
      return this->alpha<0>(ADart);
    }
    Dart_const_handle other_orientation(Dart_const_handle ADart) const
    {
      CGAL_assertion(!this->template is_free<0>(ADart));
      return this->alpha<0>(ADart);
    }

    size_type number_of_halfedges() const
    {
      assert(is_without_boundary(0));
      return number_of_darts()/2;
    }

    bool are_all_faces_closed() const
    {
      for ( typename Dart_const_range::const_iterator it(darts().begin()),
              itend(darts().end()); it!=itend; ++it)
      {
        if (this->template is_free<0>(it) ||
            this->template is_free<1>(it))
          return false;
      }

      return true;
    }

    /** Count the number of used marks.
     * @return the number of used marks.
     */
    size_type number_of_used_marks() const
    { return mnb_used_marks; }

    /** Test if a given mark is reserved.
     *  @return true iff the mark is reserved (ie in used).
     */
    bool is_reserved(size_type amark) const
    {
      CGAL_assertion(amark<NB_MARKS);
      return (mnb_times_reserved_marks[amark]!=0);
    }

    /**  Count the number of marked darts for a given mark.
     * @param amark the mark index.
     * @return the number of marked darts for amark.
     */
    size_type number_of_marked_darts(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      return mnb_marked_darts[amark];
    }

    /**  Count the number of unmarked darts for a given mark.
     * @param amark the mark index.
     * @return the number of unmarked darts for amark.
     */
    size_type number_of_unmarked_darts(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      return number_of_darts() - number_of_marked_darts(amark);
    }

    /** Test if all the darts are unmarked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are unmarked for amark.
     */
    bool is_whole_map_unmarked(size_type amark) const
    { return number_of_marked_darts(amark) == 0; }

    /** Test if all the darts are marked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are marked for amark.
     */
    bool is_whole_map_marked(size_type amark) const
    {  return number_of_marked_darts(amark) == number_of_darts(); }

    /** Reserve a new mark.
     * Get a new free mark and return its index.
     * All the darts are unmarked for this mark.
     * @return the index of the new mark.
     * @pre mnb_used_marks < NB_MARKS
     */
    size_type get_new_mark() const
    {
      if (mnb_used_marks == NB_MARKS)
      {
        std::cerr << "Not enough Boolean marks: "
          "increase NB_MARKS in item class." << std::endl;
        std::cerr << "  (exception launched)" << std::endl;
        throw Exception_no_more_available_mark();
      }

      size_type m = mfree_marks_stack[mnb_used_marks];
      mused_marks_stack[mnb_used_marks] = m;

      mindex_marks[m] = mnb_used_marks;
      mnb_times_reserved_marks[m]=1;

      ++mnb_used_marks;
      CGAL_assertion(is_whole_map_unmarked(m));

      return m;
    }

    /** Increase the number of times a mark is reserved.
     *  @param amark the mark to share.
     */
    void share_a_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      ++mnb_times_reserved_marks[amark];
    }

    /** @return the number of times a mark is reserved.
     *  @param amark the mark to share.
     */
    size_type get_number_of_times_mark_reserved(size_type amark) const
    {
      CGAL_assertion( amark<NB_MARKS );
      return mnb_times_reserved_marks[amark];
    }

    /** Negate the mark of all the darts for a given mark.
     * After this call, all the marked darts become unmarked and all the
     * unmarked darts become marked (in constant time operation).
     * @param amark the mark index
     */
    void negate_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      mnb_marked_darts[amark] = number_of_darts() - mnb_marked_darts[amark];

      mmask_marks.flip(amark);
    }

    /** Test if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(Dart_const_handle adart, size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      return get_dart_mark(adart, amark)!=mmask_marks[amark];
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(Dart_const_handle adart, size_type amark,
                     bool astate) const
    {
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark) != astate)
      {
        if (astate) ++mnb_marked_darts[amark];
        else --mnb_marked_darts[amark];

        flip_dart_mark(adart, amark);
      }
    }

    /** Mark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void mark(Dart_const_handle adart, size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark)) return;

      ++mnb_marked_darts[amark];
      flip_dart_mark(adart, amark);
    }

    /** Unmark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void unmark(Dart_const_handle adart, size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if (!is_marked(adart, amark)) return;

      --mnb_marked_darts[amark];
      flip_dart_mark(adart, amark);
    }

    /** Unmark all the darts of the map for a given mark.
     * If all the darts are marked or unmarked, this operation takes O(1)
     * operations, otherwise it traverses all the darts of the map.
     * @param amark the given mark.
     */
    void unmark_all(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if ( is_whole_map_marked(amark) )
      {
        negate_mark(amark);
      }
      else if ( !is_whole_map_unmarked(amark) )
      {
        for ( typename Dart_range::const_iterator it(darts().begin()),
               itend(darts().end()); it!=itend; ++it)
          unmark(it, amark);
      }
      CGAL_assertion(is_whole_map_unmarked(amark));
    }

    /** Free a given mark, previously calling unmark_all_darts.
     * @param amark the given mark.
     */
    void free_mark(size_type amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if ( mnb_times_reserved_marks[amark]>1 )
      {
        --mnb_times_reserved_marks[amark];
        return;
      }

      unmark_all(amark);

      // 1) We remove amark from the array mused_marks_stack by
      //    replacing it with the last mark in this array.
      mused_marks_stack[mindex_marks[amark]] =
        mused_marks_stack[--mnb_used_marks];
      mindex_marks[mused_marks_stack[mnb_used_marks]] =
        mindex_marks[amark];

      // 2) We add amark in the array mfree_marks_stack and update its index.
      mfree_marks_stack[ mnb_used_marks ] = amark;
      mindex_marks[amark] = mnb_used_marks;

      mnb_times_reserved_marks[amark]=0;
    }

    /** Test if this map is without boundary for a given dimension.
     * @param i the dimension.
     * @return true iff all the darts are not i-free.
     * @pre 0<=i<=n
     */
    bool is_without_boundary(unsigned int i) const
    {
      CGAL_assertion(i<=dimension);
      for ( typename Dart_const_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
        if (is_free(it, i)) return false;
      return true;
    }

    /** Test if this map is without boundary for all the dimensions.
     * @return true iff all the darts are non free.
     */
    bool is_without_boundary() const
    {
      for ( typename Dart_const_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
        for ( unsigned int i = 0; i<=dimension; ++i)
          if (is_free(it, i)) return false;
      return true;
    }

    /** Close the generalized map for a given dimension.
     *  @param i the dimension to close
     *  @return the number of new darts.
     *  @pre 0<=i<=n
     */
    template<unsigned int i>
    unsigned int close()
    {
      CGAL_assertion( i<=dimension );
      unsigned int res = 0;
      Dart_handle d, d2;

      for ( typename Dart_range::iterator it(darts().begin());
           it!=darts().end(); ++it)
      {
        if ( this->template is_free<i>(it) )
        {
          d = create_dart();
          ++res;

          link_alpha<i>(it, d);

          for ( unsigned int j=0; j<=dimension; ++j)
          {
            if ( j+1!=i && j!=i && j!=i+1 &&
                 !is_free(it, j) && !this->template is_free<i>(alpha(it, j)) )
            {
              basic_link_alpha(d, alpha(it, j, i), j);
            }
          }

          if ( i>0 )
          {
            d2 = alpha<i-1>(it);
            while ( !this->template is_free<i>(d2) &&
                    !this->template is_free<i-1>(alpha<i>(d2)) )
            { d2 = alpha<i, i-1>(d2); }
            if (!this->template is_free<i>(d2))
            {
              basic_link_alpha<i-1>(alpha<i>(d2), d);
            }
          }
        }
      }
      return res;
    }

    /** Test if the map is valid.
     * @return true iff the map is valid.
     */
    bool is_valid() const
    {
      bool valid = true;
      unsigned int i = 0, j = 0;
      std::vector<size_type> marks(dimension+1);
      for ( i=0; i<=dimension; ++i)
        marks[i] = INVALID_MARK;

      Helper::template
        Foreach_enabled_attributes<Reserve_mark_functor<Self> >::
          run(*this, marks);

      for ( typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        if ( !valid )
        { // We continue the traversal to mark all the darts.
          for ( i=0; i<=dimension; ++i)
            if (marks[i]!=INVALID_MARK) mark(it,marks[i]);
        }
        else
        {
          // Each alpha_i must be an involution
          for ( i = 0; i <= dimension; ++i)
            if (alpha(it, i, i)!=it)
            {
              std::cerr << "Map not valid: alpha(" << i
                        << ") is not an involution for dart "
                        <<darts().index(it)<< std::endl;
              valid = false;
            }

          // alpha_i o alpha_j must be an involution for j>=i+2
          for ( i = 0; i <= dimension-2; ++i)
          {
            for ( j = i + 2; j <= dimension; ++j)
              if (alpha(it, i, j)!=alpha(it, j, i))
              {
                std::cerr <<"Map not valid: alpha(" << i
                          <<") o alpha(" << j
                          <<") is not an involution for dart "
                          <<darts().index(it)<< std::endl;
                valid = false;
              }
          }
          Helper::template Foreach_enabled_attributes
            <internal::Test_is_valid_attribute_functor<Self> >::
            run(*this, it, marks, valid);
        }
      }
      for ( i=0; i<=dimension; ++i)
        if ( marks[i]!=INVALID_MARK )
        {
          CGAL_assertion( is_whole_map_marked(marks[i]) );
          free_mark(marks[i]);
        }

      return valid;
    }

    /// correct invalid attributes in the gmap
    void correct_invalid_attributes()
    {
      std::vector<size_type> marks(dimension+1);
      for ( unsigned int i=0; i<=dimension; ++i)
        marks[i] = INVALID_MARK;

      Helper::template
        Foreach_enabled_attributes<Reserve_mark_functor<Self> >::
          run(*this, marks);

      for ( typename Dart_range::iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        Helper::template Foreach_enabled_attributes
          <internal::Correct_invalid_attributes_functor<Self> >::
          run(*this, it, marks);
      }

      for ( unsigned int i=0; i<=dimension; ++i)
        if ( marks[i]!=INVALID_MARK )
        {
          CGAL_assertion( is_whole_map_marked(marks[i]) );
          free_mark(marks[i]);
        }

      Helper::template
        Foreach_enabled_attributes<internal::Cleanup_useless_attributes<Self> >::
          run(*this);
    }

    /// @return the number of darts.
    size_type number_of_darts() const
    { return mdarts.size(); }

    /// @return an estimation of the bytes used by the generalized map.
    size_type bytes() const
    {
      return mdarts.capacity() * sizeof(Dart) +
        internal::Count_bytes_all_attributes_functor<Self>::run(*this);
    }

    /** Write the content of the map: each dart and each alpha links.
     * @param os the ostream.
     * @return the ostream.
     */
    std::ostream& display_darts(std::ostream & os, bool attribs=false) const
    {
      unsigned int nb = 0;
      for ( typename Dart_range::const_iterator it=darts().begin();
           it!=darts().end(); ++it)
      {
        os << " dart " << darts().index(it)<<"; alpha[i]=";
        for ( unsigned int i=0; i<=dimension; ++i)
        {
          os << darts().index(alpha(it, i)) << ",\t";
        }
        if ( attribs )
        {
          Helper::template Foreach_enabled_attributes
              <Display_attribute_functor<Self> >::run(*this, it);
        }
        os << std::endl;
        ++nb;
      }
      os << "Number of darts: " << nb <<"(sizeofdarts="
         <<number_of_darts()<<")" << std::endl;
      return os;
    }

    /** Write the content of each given orbit of the map.
     * @param aos the ostream.
     * @return the ostream.
     */
    template < class Ite >
    std::ostream& display_orbits(std::ostream & aos) const
    {
      CGAL_static_assertion( (boost::is_same<typename Ite::Basic_iterator,
                              Tag_true>::value) );
      unsigned int nb = 0;
      size_type amark = get_new_mark();
      for ( typename Dart_range::const_iterator it1(darts().begin()),
             itend(darts().end()); it1!=itend; ++it1)
      {
        if ( !is_marked(it1, amark) )
        {
          ++nb;
          for ( Ite it2(*this, it1, amark); it2.cont(); ++it2 )
          {
            aos << darts().index(it2) << " - " << std::flush;
            mark(it2, amark);
          }
          aos << std::endl;
        }
      }
      CGAL_assertion( is_whole_map_marked(amark) );
      free_mark(amark);
      aos << "Number of orbits: " << nb << std::endl;
      return aos;
    }

    /** Write the content of each i-cell of the map.
     * @param aos the ostream.
     * @return the ostream.
     */
    template < unsigned int i >
    std::ostream& display_cells(std::ostream & aos) const
    {
      return display_orbits<GMap_dart_const_iterator_basic_of_cell<Self,i> >
        (aos);
    }

    /** Write the number of darts and cells of the map into a given ostream.
     * @param os the ostream.
     * @return the ostream.
     */
    std::ostream& display_characteristics(std::ostream & os) const
    {
      std::vector<unsigned int> cells(dimension+2);
      for ( unsigned int i=0; i<=dimension+1; ++i)
      { cells[i]=i; }

      std::vector<unsigned int> res = count_cells(cells);
      bool orientable=is_orientable();

      os << "#Darts=" << number_of_darts();
      for ( unsigned int i=0; i<=dimension; ++i)
        os<<", #"<<i<<"-cells="<<res[i];
      os<<", #ccs="<<res[dimension+1]
        <<", orientable="<<(orientable?"true":"false");

      return os;
    }

    /// Create a new attribute.
    /// @return a handle on the new attribute.
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template<unsigned int i, typename ...Args>
    typename Attribute_handle<i>::type create_attribute(const Args&... args)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
     typename Attribute_handle<i>::type res=
       CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(args...);
     // Reinitialize the ref counting of the new attribute. This is normally
     // not required except if create_attribute is used as "copy contructor".
     this->template init_attribute_ref_counting<i>(res);
     return res;
    }
#else
    template<unsigned int i>
    typename Attribute_handle<i>::type
    create_attribute()
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace();
    }
    template<unsigned int i, typename T1>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
     typename Attribute_handle<i>::type res=
       CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1);
     this->template init_attribute_ref_counting<i>(res);
     return res;
    }
    template<unsigned int i, typename T1, typename T2>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2);
    }
    template<unsigned int i, typename T1, typename T2, typename T3>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4,
             typename T5>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                     const T5 &t5)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                     const T5 &t5, const T6 &t6)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5, t6);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                     const T5 &t5, const T6 &t6, const T7 &t7)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5, t6, t7);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7, typename T8>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                     const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5, t6, t7, t8);
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4,
             typename T5, typename T6, typename T7, typename T8, typename T9>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                     const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8,
                     const T9 &t9)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5, t6, t7, t8, t9);
    }
#endif

    /// Erase an attribute.
    /// @param h a handle to the attribute to erase.
    template<unsigned int i>
    void erase_attribute(typename Attribute_handle<i>::type h)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "erase_attribute<i> but i-attributes are disabled");
      CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).erase(h);
    }

    /// @return true if ah points to a used i-attribute (i.e. valid).
    template<unsigned int i>
    bool is_attribute_used(typename Attribute_const_handle< i >::type ah) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "is_attribute_used<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).is_used(ah);
    }

    /// @return the number of attributes.
    template <unsigned int i>
    size_type number_of_attributes() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "number_of_attributes<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).size();
    }

    /** Set the i th attribute of all the darts of a given i-cell.
     * @param adart a dart of the i-cell.
     * @param ah the vertex to set.
     */
    template<unsigned int i>
    void set_attribute(Dart_handle dh,
                       typename Attribute_handle<i>::type ah)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "set_attribute<i> but i-attributes are disabled");
      for ( typename Dart_of_cell_range<i>::iterator it(*this, dh);
            it.cont(); ++it)
      {
        this->template set_dart_attribute<i>(it, ah);
      }
    }

    /// @return a Attributes_range<i> (range through all the
    /// attributes<i> of the map).
    template<unsigned int i>
    typename Attribute_range<i>::type & attributes()
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "attributes<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers);
    }

    template<unsigned int i>
    typename Attribute_const_range<i>::type & attributes() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "attributes<i> but i-attributes are disabled");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers);
    }

    // Get the ith dynamic onsplit functor (by reference so that we can
    // modify it directly).
    template<int i>
    boost::function<void(typename Attribute_type<i>::type&,
                         typename Attribute_type<i>::type&)>&
    onsplit_functor()
    {
      CGAL_static_assertion_msg
          (Helper::template Dimension_index<i>::value>=0,
           "onsplit_functor<i> but "
           "i-attributes are disabled");

      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (m_onsplit_functors);
    }

    // Get the ith dynamic onsplit functor (by reference so that we can
    // modify it directly).
    template<int i>
    const boost::function<void(typename Attribute_type<i>::type&,
                               typename Attribute_type<i>::type&)>&
    onsplit_functor() const
    {
      CGAL_static_assertion_msg
          (Helper::template Dimension_index<i>::value>=0,
           "onsplit_functor<i> but "
           "i-attributes are disabled");

      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (m_onsplit_functors);
    }

    // Get the ith dynamic onmerge functor (by reference so that we can
    // modify it directly).
    template<int i>
    boost::function<void(typename Attribute_type<i>::type&,
                               typename Attribute_type<i>::type&)>&
    onmerge_functor()
    {
      CGAL_static_assertion_msg
          (Helper::template Dimension_index<i>::value>=0,
           "onsplit_functor<i> but "
           "i-attributes are disabled");

      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (m_onmerge_functors);
    }
    // Get the ith dynamic onmerge functor (by reference so that we can
    // modify it directly).
    template<int i>
    const boost::function<void(typename Attribute_type<i>::type&,
                               typename Attribute_type<i>::type&)>&
    onmerge_functor() const
    {
      CGAL_static_assertion_msg
          (Helper::template Dimension_index<i>::value>=0,
           "onsplit_functor<i> but "
           "i-attributes are disabled");

      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (m_onmerge_functors);
    }

    /** Double link a dart with alpha i to a second dart.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i-linked
     * with \em adart1. Attributes are not updated, thus we can obtain
     * a non-valid map with darts belonging to a same orbit and having
     * different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param i the dimension of the alpha.
     */
    template<unsigned int i>
    void basic_link_alpha(Dart_handle dart1, Dart_handle dart2)
    {
      CGAL_assertion( i<=dimension );
      this->template dart_link_alpha<i>(dart1, dart2);
      this->template dart_link_alpha<i>(dart2, dart1);
    }
    void basic_link_alpha(Dart_handle dart1, Dart_handle dart2, unsigned int i)
    {
      CGAL_assertion( i<=dimension );
      dart_link_alpha(dart1, dart2, i);
      dart_link_alpha(dart2, dart1, i);
    }

    /** Double link two darts, and update the NULL attributes.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa.
     * If both darts have an attribute, the attribute of adart1 is
     * associated to adart2.
     * We can obtain a non-valid map with darts belonging to a same cell
     * and having different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param i the dimension of the alpha.
     * @pre i<=dimension.
     */
    template<unsigned int i>
    void link_alpha(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( i<=dimension );
      Helper::template Foreach_enabled_attributes_except
        <internal::GMap_group_attribute_functor_of_dart<Self, i>, i>::
        run(*this,adart1,adart2);
      this->template dart_link_alpha<i>(adart1, adart2);
      this->template dart_link_alpha<i>(adart2, adart1);
    }

    /** Double link a dart with alphai to a second dart.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i^-1-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa,
     * if both darts have an attribute, the attribute of adart1 is
     * associated to adart2 (only if update_attributes==true).
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param update_attributes a boolean to update the enabled attributes.
     *       (deprecated, now we use are_attributes_automatically_managed())
     */
    template<unsigned int i>
    void link_alpha(Dart_handle adart1, Dart_handle adart2,
                   bool update_attributes)
    {
      if ( update_attributes ) link_alpha<i>(adart1, adart2);
      else basic_link_alpha<i>(adart1, adart2);
    }

    /** Double unlink a dart with alpha i.
     * alphai(\em adart) is i-unlinked and \em adart is i-unlinked.
     * The attributes are not updated, thus we can obtain a non-valid map
     * with darts belonging to different orbits and having the same
     * attributes.
     * @param adart a dart.
     * @param i the dimension of the alpha.
     */
    template<unsigned int i>
    void unlink_alpha(Dart_handle adart)
    {
      CGAL_assertion(!this->template is_free<i>(adart));
      CGAL_assertion(i<=dimension);
      this->template dart_unlink_alpha<i>(alpha<i>(adart));
      this->template dart_unlink_alpha<i>(adart);
    }
    void unlink_alpha(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(!is_free(adart,i));
      CGAL_assertion(i<=dimension);
      dart_unlink_alpha(alpha(adart, i), i);
      dart_unlink_alpha(adart, i);
    }

    /** Test if it is possible to sew by alphai the two given darts
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff \em adart1 can be i-sewn with \em adart2.
     */
    template<unsigned int i>
    bool is_sewable(Dart_const_handle adart1, Dart_const_handle adart2) const
    {
      return CGAL::internal::
          GMap_is_sewable_functor<Self, i>::run(*this, adart1, adart2);
    }

    /** Topological sew by alphai two given darts plus all the required darts
     * to satisfy the generalized map validity: but do not update attributes
     * thus the map can be non valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre i<=dimension.
     * @pre is_sewable<i>(adart1, adart2).
     */
    template<unsigned int i>
    void topo_sew(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( i<=Self::dimension );
      CGAL_assertion( (is_sewable<i>(adart1,adart2)) );

      CGAL::GMap_dart_iterator_of_involution<Self,i> I1(*this, adart1);
      CGAL::GMap_dart_iterator_of_involution<Self,i> I2(*this, adart2);
      for ( ; I1.cont();  ++I1, ++I2 )
      {
        basic_link_alpha<i>(I1, I2);
      }
    }

    /** Sew by alphai the two given darts plus all the required darts
     * to satisfy the generalized map validity, and updates enabled
     * attributes when necessary so that the final map is valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable<i>(adart1, adart2).
     * @pre i<=dimension.
     * @post is_valid()
     */
    template<unsigned int i>
    void sew(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( i<=dimension );
      CGAL_assertion( (is_sewable<i>(adart1,adart2)) );

      size_type amark=get_new_mark();

      CGAL::GMap_dart_iterator_basic_of_involution<Self, i>
          I1(*this, adart1, amark);
      CGAL::GMap_dart_iterator_basic_of_involution<Self, i>
          I2(*this, adart2, amark);

      // This first loop do not modify the map, but only the attributes
      // (by calling when required the onmerge functors).
      for ( ; I1.cont(); ++I1, ++I2 )
      {
        Helper::template Foreach_enabled_attributes_except
            <CGAL::internal::GMap_group_attribute_functor<Self, i>, i>::
            run(*this, I1, I2);
      }

      // Now we update the alpha links.
      negate_mark( amark );
      for ( I1.rewind(), I2.rewind(); I1.cont(); ++I1, ++I2 )
      {
        basic_link_alpha<i>(I1, I2);
      }

      negate_mark( amark );
      CGAL_assertion( is_whole_map_unmarked(amark) );
      free_mark(amark);
    }

    /** Sew by alphai the two given darts plus all the required darts
     * to satisfy the generalized map validity. Enabled attributes
     * are updated only if update_attributes==true.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @param update_attributes a boolean to update the enabled attributes
     *       (deprecated, now we use are_attributes_automatically_managed())
     * @pre is_sewable<i>(adart1, adart2).
     */
    template<unsigned int i>
    void sew(Dart_handle adart1, Dart_handle adart2, bool update_attributes)
    {
      if ( update_attributes ) sew<i>(adart1, adart2);
      else topo_sew<i>(adart1, adart2);
    }

    /** Topological unsew by alphai the given dart plus all the required darts
     * to satisfy the generalized map validity: but do not update attributes
     * thus the map can be non valid
     * @param adart first dart.
     * @pre !adart->is_free(i).
     * @pre i<=dimension.
     */
    template<unsigned int i>
    void topo_unsew(Dart_handle adart)
    {
      CGAL_assertion( !(is_free<i>(adart)) );
      CGAL_assertion( i<=Self::dimension );

      for ( CGAL::GMap_dart_iterator_of_involution<Self,i> it(*this, adart);
            it.cont(); ++it )
      { unlink_alpha<i>(*it); }
    }

    /** Unsew by alphai the given dart plus all the required darts
     * to satisfy the generalized map validity, and update enabled
     * attributes when necessary so that the final map is valid.
     * @param adart first dart.
     * @pre !adart->is_free(i).
     * @post is_valid()
     * @pre i<=dimension
     */
    template<unsigned int i>
    void unsew(Dart_handle adart)
    {
      CGAL_assertion(i<=Self::dimension);
      CGAL_assertion( !this->template is_free<i>(adart) );

      std::deque<Dart_handle> modified_darts;

      for ( CGAL::GMap_dart_iterator_of_involution<Self, i> it(*this, adart);
            it.cont(); ++it )
      {
        modified_darts.push_back(it);
        modified_darts.push_back(alpha<i>(it));
        unlink_alpha<i>(it);
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::GMap_test_split_attribute_functor<Self, i>, i>::
          run(*this, modified_darts);
    }

    /** Unsew by alphai the given dart plus all the required darts
     * to satisfy the generalized map validity. Enabled attributes
     * are updated only if update_attributes==true.
     * @param adart first dart.
     * @param update_attributes a boolean to update the enabled attributes
     *       (deprecated, now we use are_attributes_automatically_managed())
     * @pre !adart->is_free(i).
     */
    template<unsigned int i>
    void unsew(Dart_handle adart, bool update_attributes)
    {
      if ( update_attributes ) unsew<i>(adart);
      else topo_unsew<i>(adart);
    }

    /** Count the marked cells (at least one marked dart).
     * @param amark the mark to consider.
     * @param avector containing the dimensions of the cells to count.
     * @return a vector containing the number of cells.
     */
    std::vector<unsigned int>
    count_marked_cells(size_type amark, const std::vector<unsigned int>& acells) const
    {
      std::vector<unsigned int> res(dimension+2);
      std::vector<size_type> marks(dimension+2);

      // Initialization of the result
      for ( unsigned int i=0; i<dimension+2; ++i)
      {
        res[i]=0;
        marks[i]=INVALID_MARK;
      }

      // Mark reservation
      for ( unsigned int i=0; i<acells.size(); ++i)
      {
        CGAL_assertion(acells[i]<=dimension+1);
        if ( marks[acells[i]]==INVALID_MARK )
        {
          marks[acells[i]] = get_new_mark();
          assert(is_whole_map_unmarked(marks[acells[i]]));
        }
      }

      // Counting and marking cells
      for ( typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        if ( is_marked(it, amark) )
        {
          CGAL::internal::Foreach_static
            <CGAL::internal::Count_cell_functor<Self>,dimension+1>::
            run(*this, it, marks, res);
        }
      }

      // Unmarking darts
      std::vector<size_type> tounmark;
      for ( unsigned int i=0; i<acells.size(); ++i)
      {
        if ( is_whole_map_marked(marks[acells[i]]) ||
             is_whole_map_unmarked(marks[acells[i]]))
        {
          free_mark(marks[acells[i]]);
        }
        else
        {
          tounmark.push_back(marks[acells[i]]);
        }
      }

      if ( tounmark.size() > 0 )
      {
        for ( typename Dart_range::const_iterator it(darts().begin()),
               itend(darts().end()); it!=itend; ++it)
        {
          for ( unsigned int i=0; i<tounmark.size(); ++i)
            unmark(it, tounmark[i]);
        }
        for ( unsigned int i=0; i<tounmark.size(); ++i)
        {
          CGAL_assertion(is_whole_map_unmarked(tounmark[i]));
          free_mark(tounmark[i]);
        }
      }

      return res;
    }

    /** Count the number of given cells
     * @param avector containing the dimensions of the cells to count.
     * @return a vector containing the number of cells.
     */
    std::vector<unsigned int>
    count_cells(const std::vector<unsigned int>& acells) const
    {
      std::vector<unsigned int> res;
      size_type m = get_new_mark();
      negate_mark(m); // We mark all the cells.

      res = count_marked_cells(m, acells);

      negate_mark(m); // We unmark the cells
      free_mark(m);

      return res;
    }

    /** Count the number of cells in each dimension.
     * @return a vector containing the number of cells.
     */
    std::vector<unsigned int> count_all_cells() const
    {
      std::vector<unsigned int> dim(dimension+2);

      for ( unsigned int i=0; i<dimension+2; ++i)
        dim[i]=i;

      return count_cells(dim);
    }

  protected:
    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(Dart_const_handle adart,
                   const std::bitset<NB_MARKS> & amarks) const
    { set_dart_marks(adart, amarks ^ mmask_marks); }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(Dart_const_handle adart) const
    { return get_dart_marks(adart) ^ mmask_marks; }

    /** Get the mask associated to a given mark.
     * @param amark the mark.
     * @return the mask associated to mark amark.
     */
    bool get_mask_mark(size_type amark) const
    {
      CGAL_assertion(amark>=0 && amark<NB_MARKS);
      return mmask_marks[amark];
    }

  public:
    /** Erase marked darts from the map.
     * Marked darts are unlinked before to be removed, thus surviving darts
     * are correctly linked, but the map is not necessarily valid depending
     * on the configuration of removed darts. User must check carefully marked
     * darts before calling this method.
     * @param amark the mark of darts to erase.
     * @return the number of removed darts.
     */
    unsigned int erase_marked_darts(size_type amark)
    {
      unsigned int res = 0, i = 0;
      Dart_handle d;
      for ( typename Dart_range::iterator it(darts().begin()),
             itend(darts().end()); it!=itend; )
      {
        d = it++;
        if (is_marked(d, amark))
        {
          for ( i = 0; i <= dimension; ++i)
          { if (!is_free(d, i)) unlink_beta(d, i); }
          erase_dart(d); ++res;
        }
      }
      return res;
    }

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    //**************************************************************************
    // Dart_of_orbit_basic_range
    template<unsigned int ... Alpha>
    struct Dart_of_orbit_basic_range : public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_basic_of_orbit<Self,Alpha...>,
     CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,Alpha...> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_basic_of_orbit<Self,Alpha...>,
       CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,Alpha...> > Base;

      Dart_of_orbit_basic_range(Self &amap, Dart_handle adart, size_type amark=INVALID_MARK):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_basic_const_range
    template<unsigned int ... Alpha>
    struct Dart_of_orbit_basic_const_range : public CGAL::CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,Alpha...> >
    {
      typedef CGAL::CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,Alpha...> >
      Base;

      Dart_of_orbit_basic_const_range(const Self &amap, Dart_const_handle
                                      adart, size_type amark=INVALID_MARK):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_range
    template<unsigned int ... Alpha>
    struct Dart_of_orbit_range : public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_of_orbit<Self,Alpha...>,
     CGAL::GMap_dart_const_iterator_of_orbit<Self,Alpha...> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_of_orbit<Self,Alpha...>,
       CGAL::GMap_dart_const_iterator_of_orbit<Self,Alpha...> > Base;

      Dart_of_orbit_range(Self &amap, Dart_handle adart) : Base(amap,adart)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_const_range
    template<unsigned int ... Alpha>
    struct Dart_of_orbit_const_range : public CGAL::CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_of_orbit<Self,Alpha...> >
    {
      typedef CGAL::CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_of_orbit<Self,Alpha...> > Base;

      Dart_of_orbit_const_range(const Self &amap, Dart_const_handle adart):
        Base(amap,adart)
      {}
    };
    //**************************************************************************
    /// @return a range on all the darts of the given orbit
    template<unsigned int ... Alpha>
    Dart_of_orbit_range<Alpha...> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<Alpha...>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Alpha>
    Dart_of_orbit_const_range<Alpha...>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<Alpha...>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Alpha>
    Dart_of_orbit_basic_range<Alpha...> darts_of_orbit_basic(Dart_handle adart,
                                                            size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<Alpha...>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Alpha>
    Dart_of_orbit_basic_const_range<Alpha...>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<Alpha...>(*this,adart,amark); }
    //**************************************************************************
#else
    //**************************************************************************
    // Dart_of_orbit_basic_range
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>
    struct Dart_of_orbit_basic_range: public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
        B8,B9>,
     CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
        B8,B9> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
      B8,B9>,
       CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,
                                               B6,B7,B8,B9> > Base;

      Dart_of_orbit_basic_range(Self &amap, Dart_handle adart,
                                size_type /*amark*/=INVALID_MARK):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_basic_const_range
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>
    struct Dart_of_orbit_basic_const_range: public CMap_const_range
    <Self,
     CGAL::GMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
        B8,B9> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_basic_of_orbit
       <Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> > Base;

      Dart_of_orbit_basic_const_range(const Self &amap,
                                      Dart_const_handle adart, size_type amark=INVALID_MARK):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_range
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>
    struct Dart_of_orbit_range: public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9>,
     CGAL::GMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
      B8,B9>,
       CGAL::GMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
      B8,B9> >
      Base;

      Dart_of_orbit_range(Self &amap, Dart_handle adart):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_const_range
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>
    struct Dart_of_orbit_const_range: public CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,
        B8,B9> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_of_orbit
       <Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> > Base;

      Dart_of_orbit_const_range(const Self &amap, Dart_const_handle adart):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    /// @return a range on all the darts of the given orbit
    Dart_of_orbit_range<> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_range<B1> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_range<B1,B2> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_range<B1,B2,B3> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_range<B1,B2,B3,B4> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_range<B1,B2,B3,B4,B5> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4,B5>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_range<B1,B2,B3,B4,B5,B6> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8> darts_of_orbit
    (Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
    darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>(*this,adart); }
    //--------------------------------------------------------------------------
    // Const versions.
    Dart_of_orbit_const_range<> darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_const_range<B1> darts_of_orbit(Dart_const_handle
                                                 adart) const
    { return Dart_of_orbit_const_range<B1>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_const_range<B1,B2> darts_of_orbit(Dart_const_handle
                                                    adart) const
    { return Dart_of_orbit_const_range<B1,B2>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_const_range<B1,B2,B3> darts_of_orbit
    (Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_const_range<B1,B2,B3,B4>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_const_range<B1,B2,B3,B4,B5>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4,B5>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8>(*this,adart); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
        (*this,adart); }
    //--------------------------------------------------------------------------
    // Basic versions
    Dart_of_orbit_basic_range<> darts_of_orbit_basic(Dart_handle adart,
                                                     size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    Dart_of_orbit_basic_const_range<> darts_of_orbit_basic
    (Dart_const_handle adart,size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_basic_range<B1> darts_of_orbit_basic(Dart_handle adart,
                                                       size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_basic_const_range<B1> darts_of_orbit_basic
    (Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_basic_range<B1,B2> darts_of_orbit_basic(Dart_handle adart,
                                                          size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_basic_const_range<B1,B2> darts_of_orbit_basic
    (Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_basic_range<B1,B2,B3> darts_of_orbit_basic(Dart_handle adart,
                                                             size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_basic_const_range<B1,B2,B3> darts_of_orbit_basic
    (Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_basic_range<B1,B2,B3,B4> darts_of_orbit_basic
    (Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5> darts_of_orbit_basic
    (Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6> darts_of_orbit_basic
    (Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7> darts_of_orbit_basic
    (Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8> darts_of_orbit
    (Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
    darts_of_orbit_basic(Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
    darts_of_orbit_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
        (*this,adart,amark); }
    //**************************************************************************
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    //**************************************************************************
    // Dart_of_cell_basic_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_basic_range: public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_basic_of_cell<Self,i,dim>,
     CGAL::GMap_dart_const_iterator_basic_of_cell<Self,i,dim> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_basic_of_cell<Self,i,dim>,
       CGAL::GMap_dart_const_iterator_basic_of_cell<Self,i,dim> > Base;

      Dart_of_cell_basic_range(Self &amap, Dart_handle adart, size_type amark=INVALID_MARK) :
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_basic_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_basic_const_range: public CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_basic_of_cell<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_basic_of_cell<Self,i,dim> > Base;

      Dart_of_cell_basic_const_range(const Self &amap, Dart_const_handle adart,
                                     size_type amark=INVALID_MARK) :
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_range: public CGAL::CMap_range
    <Self,GMap_dart_iterator_of_cell<Self,i,dim>,
     CGAL::GMap_dart_const_iterator_of_cell<Self,i,dim> >
    {
      typedef CGAL::CMap_range
      <Self,GMap_dart_iterator_of_cell<Self,i,dim>,
       CGAL::GMap_dart_const_iterator_of_cell<Self,i,dim> > Base;

      Dart_of_cell_range(Self &amap, Dart_handle adart) :
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_const_range: public CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_of_cell<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_of_cell<Self,i,dim> > Base;

      Dart_of_cell_const_range(const Self &amap, Dart_const_handle adart) :
        Base(amap, adart)
      {}
    };
    //--------------------------------------------------------------------------
    /// @return a range on all the darts of the given i-cell
    template<unsigned int i, int dim>
    Dart_of_cell_basic_range<i,dim> darts_of_cell_basic(Dart_handle adart,
                                                        size_type amark=INVALID_MARK)
    { return Dart_of_cell_basic_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i, int dim>
    Dart_of_cell_basic_const_range<i,dim> darts_of_cell_basic
    (Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_cell_basic_const_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_basic_range<i,Self::dimension>
    darts_of_cell_basic(Dart_handle adart, size_type amark=INVALID_MARK)
    { return darts_of_cell_basic<i,Self::dimension>(adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_basic_const_range<i,Self::dimension>
    darts_of_cell_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return darts_of_cell_basic<i,Self::dimension>(adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i, int dim>
    Dart_of_cell_range<i,dim> darts_of_cell(Dart_handle adart)
    { return Dart_of_cell_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i, int dim>
    Dart_of_cell_const_range<i,dim> darts_of_cell(Dart_const_handle adart) const
    { return Dart_of_cell_const_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_range<i,Self::dimension> darts_of_cell(Dart_handle adart)
    { return darts_of_cell<i,Self::dimension>(adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_const_range<i,Self::dimension>
    darts_of_cell(Dart_const_handle adart) const
    { return darts_of_cell<i,Self::dimension>(adart); }
    //**************************************************************************
    // Dart_of_involution_basic_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_basic_range: public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_basic_of_involution<Self,i,dim>,
     CGAL::GMap_dart_const_iterator_basic_of_involution<Self,i,dim> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_basic_of_involution<Self,i,dim>,
       CGAL::GMap_dart_const_iterator_basic_of_involution<Self,i,dim> > Base;

      Dart_of_involution_basic_range(Self &amap, Dart_handle adart,
                                     size_type amark=INVALID_MARK):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_basic_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_basic_const_range: public CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_basic_of_involution<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_basic_of_involution<Self,i,dim> >
      Base;

      Dart_of_involution_basic_const_range(const Self &amap,
                                           Dart_const_handle adart,
                                           size_type amark=INVALID_MARK) :
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    template<unsigned int i,int dim>
    Dart_of_involution_basic_range<i,dim>
    darts_of_involution_basic(Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_involution_basic_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i,int dim>
    Dart_of_involution_basic_const_range<i,dim>
    darts_of_involution_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_involution_basic_const_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_basic_range<i,Self::dimension>
    darts_of_involution_basic(Dart_handle adart, size_type amark=INVALID_MARK)
    { return Dart_of_involution_basic_range<i,Self::dimension>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_basic_const_range<i,Self::dimension>
    darts_of_involution_basic(Dart_const_handle adart, size_type amark=INVALID_MARK) const
    { return Dart_of_involution_basic_const_range<i,Self::dimension>
        (*this,adart,amark); }
    //**************************************************************************
    // Dart_of_involution_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_range: public CGAL::CMap_range
    <Self, CGAL::GMap_dart_iterator_of_involution<Self,i,dim>,
     CGAL::GMap_dart_const_iterator_of_involution<Self,i,dim> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_dart_iterator_of_involution<Self,i,dim>,
       CGAL::GMap_dart_const_iterator_of_involution<Self,i,dim> > Base;

      Dart_of_involution_range(Self &amap, Dart_handle adart) :
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_const_range: public CMap_const_range
    <Self, CGAL::GMap_dart_const_iterator_of_involution<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_dart_const_iterator_of_involution<Self,i,dim> > Base;

      Dart_of_involution_const_range(const Self &amap,
                                     Dart_const_handle adart):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    template<unsigned int i,int dim>
    Dart_of_involution_range<i,dim>
    darts_of_involution(Dart_handle adart)
    { return Dart_of_involution_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i,int dim>
    Dart_of_involution_const_range<i,dim>
    darts_of_involution(Dart_const_handle adart) const
    { return Dart_of_involution_const_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_range<i,Self::dimension>
    darts_of_involution(Dart_handle adart)
    { return Dart_of_involution_range<i,Self::dimension>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_const_range<i,Self::dimension>
    darts_of_involution(Dart_const_handle adart) const
    { return Dart_of_involution_const_range<i,Self::dimension>(*this,adart); }
    //**************************************************************************
    // Dart_basic_range
    struct Dart_basic_range {
      typedef CGAL::GMap_dart_iterator_basic_of_all<Self> iterator;
      typedef CGAL::GMap_dart_const_iterator_basic_of_all<Self> const_iterator;
      Dart_basic_range(Self &amap) : mmap(amap)
      {}
      iterator begin() { return iterator(mmap); }
      iterator end()   { return iterator(mmap,mmap.null_handle); }
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
      size_type size() const
      { return mmap.number_of_darts(); }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      Self & mmap;
    };
    //**************************************************************************
    // Dart_basic_const_range
    struct Dart_basic_const_range {
      typedef CGAL::GMap_dart_const_iterator_basic_of_all<Self> const_iterator;
      Dart_basic_const_range(Self &amap) : mmap(amap)
      {}
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
      size_type size() const
      { return mmap.number_of_darts(); }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      const Self & mmap;
    };
    //**************************************************************************
    Dart_basic_range darts_basic()
    { return Dart_basic_range(*this); }
    //--------------------------------------------------------------------------
    Dart_basic_const_range darts_basic() const
    { return Dart_basic_const_range(*this); }
    //**************************************************************************
    // One_dart_per_incident_cell_range
    template<unsigned int i,unsigned int j,int dim=Self::dimension>
    struct One_dart_per_incident_cell_range: public CGAL::CMap_range
    <Self, CGAL::GMap_one_dart_per_incident_cell_iterator<Self,i,j,dim>,
     CGAL::GMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
    {
      typedef CGAL::CMap_range
      <Self, CGAL::GMap_one_dart_per_incident_cell_iterator<Self,i,j,dim>,
       CGAL::GMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
      Base;

      One_dart_per_incident_cell_range(Self &amap, Dart_handle adart):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // One_dart_per_incident_cell_const_range
    template<unsigned int i,unsigned int j,int dim=Self::dimension>
    struct One_dart_per_incident_cell_const_range: public CMap_const_range
    <Self, CGAL::GMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
    {
      typedef CMap_const_range
      <Self, CGAL::GMap_one_dart_per_incident_cell_const_iterator
      <Self,i,j,dim> > Base;

      One_dart_per_incident_cell_const_range(const Self &amap,
                                             Dart_const_handle adart) :
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // One_dart_per_cell_range
    template<unsigned int i,int dim=Self::dimension>
    struct One_dart_per_cell_range {
      typedef CGAL::GMap_one_dart_per_cell_iterator<Self,i,dim> iterator;
      typedef CGAL::GMap_one_dart_per_cell_const_iterator<Self,i,dim>
      const_iterator;
      One_dart_per_cell_range(Self &amap) : mmap(amap), msize(0)
      {}
      iterator begin() { return iterator(mmap); }
      iterator end()   { return iterator(mmap,mmap.null_handle); }
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
      size_type size() const
      {
        if (msize==0)
          for ( const_iterator it=begin(), itend=end(); it!=itend; ++it)
            ++msize;
        return msize;
      }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      Self & mmap;
      mutable size_type msize;
    };
    //**************************************************************************
    // One_dart_per_cell_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct One_dart_per_cell_const_range {
      typedef CGAL::GMap_one_dart_per_cell_const_iterator<Self,i,dim>
      const_iterator;
      One_dart_per_cell_const_range(const Self &amap) : mmap(amap), msize(0)
      {}
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
      size_type size() const
      {
        if (msize==0)
          for ( const_iterator it=begin(), itend=end(); it!=itend; ++it)
            ++msize;
        return msize;
      }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      const Self & mmap;
      mutable size_type msize;
    };
    //**************************************************************************
    /// @return a range on the i-cells incindent to the given j-cell.
    template<unsigned int i, unsigned int j, int dim>
    One_dart_per_incident_cell_range<i,j,dim>
    one_dart_per_incident_cell(Dart_handle adart)
    { return One_dart_per_incident_cell_range<i,j,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i, unsigned int j, int dim>
    One_dart_per_incident_cell_const_range<i,j,dim>
    one_dart_per_incident_cell(Dart_const_handle adart) const
    { return One_dart_per_incident_cell_const_range<i,j,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i, unsigned int j>
    One_dart_per_incident_cell_range<i,j,Self::dimension>
    one_dart_per_incident_cell(Dart_handle adart)
    { return one_dart_per_incident_cell<i,j,Self::dimension>(adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i, unsigned int j>
    One_dart_per_incident_cell_const_range<i,j,Self::dimension>
    one_dart_per_incident_cell(Dart_const_handle adart) const
    { return one_dart_per_incident_cell<i,j,Self::dimension>(adart); }
    //--------------------------------------------------------------------------
    /// @return a range on all the i-cells
    template<unsigned int i, int dim>
    One_dart_per_cell_range<i,dim> one_dart_per_cell()
    { return One_dart_per_cell_range<i,dim>(*this); }
    //--------------------------------------------------------------------------
    template<unsigned int i, int dim>
    One_dart_per_cell_const_range<i,dim> one_dart_per_cell() const
    { return One_dart_per_cell_const_range<i,dim>(*this); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    One_dart_per_cell_range<i,Self::dimension> one_dart_per_cell()
    { return one_dart_per_cell<i,Self::dimension>(); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    One_dart_per_cell_const_range<i,Self::dimension> one_dart_per_cell() const
    { return one_dart_per_cell<i,Self::dimension>(); }
    //--------------------------------------------------------------------------

  public:

    /** Compute the dual of a Generalized_map.
     * @param amap the gmap in which we build the dual of this map.
     * @param adart a dart of the initial map, NULL by default.
     * @return adart of the dual map, the dual of adart if adart!=NULL,
     *         any dart otherwise.
     * As soon as we don't modify this map and amap map, we can iterate
     * simultaneously through all the darts of the two maps and we have
     * each time of the iteration two "dual" darts.
     */
    Dart_handle dual(Self& amap, Dart_handle adart=null_handle)
    {
      CGAL_assertion( is_without_boundary(dimension) );

      CGAL::Unique_hash_map< Dart_handle, Dart_handle,
        typename Self::Hash_function > dual;
      Dart_handle d, d2, res = amap.null_handle;

      // We clear amap. TODO return a new amap ?
      amap.clear();

      // We create a copy of all the dart of the map.
      for ( typename Dart_range::iterator it=darts().begin();
            it!=darts().end(); ++it)
      {
        dual[it] = amap.create_dart();
        internal::Copy_dart_info_functor<Refs, Refs>::
          run(static_cast<Refs&>(amap), static_cast<Refs&>(*this),
              it, dual[it]);
        if ( it==adart && res==amap.null_handle ) res = dual[it];
      }

      // Then we link the darts by using the dual formula :
      // G(B,a0,a1,...,a(n-1),an) =>
      //    dual(G)=(B, an, a(n-1),...,a1, a0)
      // We suppose darts are run in the same order for both maps.
      typename Dart_range::iterator it2=amap.darts().begin();
      for ( typename Dart_range::iterator it=darts().begin();
            it!=darts().end(); ++it, ++it2)
      {
        d = it2; // The supposition on the order allows to avoid d=dual[it];
        CGAL_assertion( it2==dual[it] );

        for ( unsigned int i=0; i<=dimension; ++i)
        {
          if ( !is_free(it, i) && amap.is_free(d,dimension-i) )
            amap.basic_link_alpha(d, dual[alpha(it, i)], dimension-i);
        }
      }

      //  CGAL_postcondition(amap2.is_valid());

      if ( res==amap.null_handle ) res = amap.darts().begin();
      return res;
    }


    /** Test if the connected component of gmap containing dart dh1 is
     *  isomorphic to the connected component of map2 containing dart dh2,
     *  starting from dh1 and dh2.
     * @param dh1  initial dart for this map
     * @param map2 the second generalized map
     * @param dh2  initial dart for map2
     * @param testDartInfo Boolean to test the equality of dart info (true)
     *                     or not (false)
     * @param testAttributes Boolean to test the equality of attributes (true)
     *                       or not (false)
     * @param testPoint Boolean to test the equality of points (true)
     *                     or not (false) (used for LCC)
     * @return true iff the cc of map is isomorphic to the cc of map2 starting
     *     from dh1 and dh2; by testing the equality of dartinfo and/or
     *     attributes and/or points.
     */
    template <unsigned int d2, typename Refs2, typename Items2, class Alloc2,
              class Storage2>
    bool are_cc_isomorphic(Dart_const_handle dh1,
                           const Generalized_map_base
                           <d2,Refs2,Items2,Alloc2, Storage2>& map2,
                           typename Generalized_map_base
                           <d2,Refs2,Items2,Alloc2, Storage2>::Dart_const_handle dh2,
                           bool testDartInfo=true,
                           bool testAttributes=true,
                           bool testPoint=true) const
    {
      typedef Generalized_map_base<d2,Refs2,Items2,Alloc2, Storage2> Map2;

      bool match = true;

      // Two stacks used to run through the two maps.
      std::deque< Dart_const_handle > toTreat1;
      std::deque< typename Map2::Dart_const_handle > toTreat2;

       // A dart of this map is marked with m1 if its bijection was set
      // (and similarly for mark m2 and darts of map2)
      size_type m1 = get_new_mark();
      size_type m2 = map2.get_new_mark();

      // A dart of this map is marked with markpush if it was already pushed
      // in the queue toTreat1.
      size_type markpush = get_new_mark();

      toTreat1.push_back(dh1);
      toTreat2.push_back(dh2);

      Dart_const_handle current;
      typename Map2::Dart_const_handle other;

      unsigned int i = 0;
      CGAL::Unique_hash_map<Dart_const_handle,
                            typename Map2::Dart_const_handle,
                            typename Self::Hash_function> bijection;

      while (match && !toTreat1.empty())
      {
        // Next dart
        current = toTreat1.front();
        toTreat1.pop_front();
        other = toTreat2.front();
        toTreat2.pop_front();

        if (!is_marked(current, m1))
        {
          if (map2.is_marked(other, m2))
            match=false;
          else
          {
            bijection[current] = other;

            mark(current, m1);
            map2.mark(other, m2);

            // We first test info of darts
            if (match && testDartInfo)
              match=internal::Test_is_same_dart_info_functor<Self, Map2>::
                  run(*this, map2, current, other);

            // We need to test in both direction because
            // Foreach_enabled_attributes only test non void attributes
            // of Self. Functor Test_is_same_attribute_functor will modify
            // the value of match to false if attributes do not match
            if (testAttributes)
            {
              if (match)
                Helper::template Foreach_enabled_attributes
                    < internal::Test_is_same_attribute_functor<Self, Map2> >::
                    run(*this, map2, current, other, match);
              if (match)
                Map2::Helper::template Foreach_enabled_attributes
                    < internal::Test_is_same_attribute_functor<Map2, Self> >::
                    run(map2, *this, other, current, match);
            }

            if (match && testPoint)
            {
              // Only point of 0-attributes are tested. TODO test point of all
              // attributes ?
              match=internal::Test_is_same_attribute_point_functor
                  <Self, Map2, 0>::run(*this, map2, current, other);
            }

            // We test if the injection is valid with its neighboors.
            // We go out as soon as it is not satisfied.
            for (i = 0; match && i <= dimension; ++i)
            {
              if ( i>map2.dimension )
              {
                if (!is_free(current,i))
                { match=false; }
              }
              else
              {
                if (is_free(current,i))
                {
                  if (!map2.is_free(other,i))
                  { match = false; }
                }
                else
                {
                  if (map2.is_free(other,i))
                  { match = false; }
                  else
                  {
                    if (is_marked(alpha(current,i), m1) !=
                        map2.is_marked(map2.alpha(other,i), m2))
                    { match = false; }
                    else
                    {
                      if (!is_marked(alpha(current,i), m1))
                      {
                        if (!is_marked(alpha(current,i), markpush))
                        {
                          toTreat1.push_back(alpha(current,i));
                          toTreat2.push_back(map2.alpha(other,i));
                          mark(alpha(current,i), markpush);
                        }
                      }
                      else
                      {
                        if (bijection[alpha(current,i)]!=map2.alpha(other,i))
                        { match = false; }
                      }
                    }
                  }
                }
              }
            }
            // Now we test if the second map has more alpha links than the first
            for ( i=dimension+1; match && i<=map2.dimension; ++i )
            {
              if (!map2.is_free(other,i))
              { match=false; }
            }
          }
        }
        else
        {
          if (!map2.is_marked(other, m2))
          { match=false; }
        }
      }

      // Here we test if both queue are empty
      if ( !toTreat1.empty() || !toTreat2.empty() )
      { match=false; }

      // Here we unmark all the marked darts.
      toTreat1.clear();
      toTreat2.clear();

      toTreat1.push_back(dh1);
      toTreat2.push_back(dh2);

      unmark(dh1, m1);
      unmark(dh1, markpush);
      map2.unmark(dh2, m2);

      while (!toTreat1.empty())
      {
        current = toTreat1.front();
        toTreat1.pop_front();
        other = toTreat2.front();
        toTreat2.pop_front();

        for (i = 0; i <= dimension; ++i)
        {
          if (!is_free(current,i) && is_marked(alpha(current,i), markpush))
          {
            toTreat1.push_back(alpha(current,i));
            toTreat2.push_back(map2.alpha(other,i));
            unmark(alpha(current,i), m1);
            unmark(alpha(current,i), markpush);
            map2.unmark(map2.alpha(other,i), m2);
          }
        }
      }

      assert(is_whole_map_unmarked(m1));
      assert(is_whole_map_unmarked(markpush));
      assert(map2.is_whole_map_unmarked(m2));
      free_mark(m1);
      free_mark(markpush);
      map2.free_mark(m2);

      return match;
    }

    /** Test if this gmap is isomorphic to map2.
     * @pre gmap is connected.
     * @param map2 the second generalized map
     * @param testDartInfo Boolean to test the equality of dart info (true)
     *                     or not (false)
     * @param testAttributes Boolean to test the equality of attributes (true)
     *                       or not (false)
     * @param testPoint Boolean to test the equality of points (true)
     *                     or not (false) (used for LCC)
     * @return true iff this map is isomorphic to map2, testing the equality
     *         of attributes if testAttributes is true
     */
    template <unsigned int d2, typename Refs2, typename Items2, class Alloc2,
              class Storage2>
    bool is_isomorphic_to(const Generalized_map_base
                          <d2,Refs2,Items2,Alloc2, Storage2>& map2,
                          bool testDartInfo=true,
                          bool testAttributes=true,
                          bool testPoint=true) const
    {
      if (is_empty() && map2.is_empty()) return true;
      if (is_empty() || map2.is_empty()) return false;

      Dart_const_handle d1=darts().begin();

      for (typename Generalized_map_base<d2,Refs2,Items2,Alloc2, Storage2>::
             Dart_range::const_iterator it(map2.darts().begin()),
             itend(map2.darts().end()); it!=itend; ++it)
      {
        if (are_cc_isomorphic(d1, map2, it, testDartInfo, testAttributes,
                              testPoint))
        {
          return true;
        }
      }

      return false;
    }

    /** @return true iff the connected component containing dh is orientable
     *  If amark!=INVALID_MARK, mark all darts in the cc with amark
     *  If aorientationmark!=INVALID_MARK, and if the cc is orientable, mark
     *       all darts of one orientation of the cc with aorientationmark
     *       (i.e. one out of two darts)
     *  Note it is not allow to have accmark!=INVALID_MARK and
     *  aorientationmark==INVALID_MARK.
     */
    bool is_cc_orientable(Dart_const_handle dh,
                          size_type accmark=INVALID_MARK,
                          size_type aorientationmark=INVALID_MARK) const
    {
      if (accmark!=INVALID_MARK && aorientationmark==INVALID_MARK)
      {
        std::cerr<<"Error for is_cc_orientable: you cannot use accmark"
                 <<" different from INVALID_MARK and aorientationmark "
                 <<" equal to INVALID_MARK"<<std::endl;
        accmark=INVALID_MARK;
      }

      size_type ccmark=(accmark==INVALID_MARK?get_new_mark():accmark);
      size_type orientationmark=
        (aorientationmark==INVALID_MARK?get_new_mark():aorientationmark);
      bool stop_if_nonorientable =
        (accmark==INVALID_MARK && aorientationmark==INVALID_MARK);

      std::queue<Dart_const_handle> toTreat;
      bool orientable=true;
      Dart_const_handle cur;

      // Iterator through the connected component
      toTreat.push(dh);
      mark(dh, ccmark);
      mark(dh, orientationmark);

      while((!stop_if_nonorientable || orientable) && !toTreat.empty())
      {
        cur=toTreat.front();
        toTreat.pop();

        for (unsigned int i=0; i<=dimension; ++i)
        {
          if (is_marked(cur, orientationmark))
          {
            if (!is_free(cur, i) &&
                is_marked(alpha(cur, i), orientationmark))
              orientable=false;
          }
          else
          {
            if (!is_free(cur, i))
              mark(alpha(cur, i), orientationmark);
          }

          if (!is_marked(alpha(cur, i), ccmark))
          {
            mark(alpha(cur, i), ccmark);
            toTreat.push(alpha(cur, i));
          }
        }
      }

      // Now we need to unmark the marked darts for accmark and/or for
      // aorientationmark (when they are different from INVALID_MARK).
      if (stop_if_nonorientable)
      {
        toTreat.push(dh);
        unmark(dh, ccmark);
        if (aorientationmark==INVALID_MARK) unmark(dh, orientationmark);

        while(!toTreat.empty())
        {
          cur=toTreat.front();
          toTreat.pop();

          for (unsigned int i=0; i<=dimension; ++i)
          {
            if (aorientationmark==INVALID_MARK)
              unmark(alpha(cur, i), orientationmark);

            if (!is_marked(alpha(cur, i), ccmark))
            {
              unmark(alpha(cur, i), ccmark);
              toTreat.push(alpha(cur, i));
            }
          }
        }

        if (aorientationmark==INVALID_MARK)
        {
          assert(is_whole_map_marked(orientationmark));
          free_mark(orientationmark);
        }

        assert(is_whole_map_marked(ccmark));
        free_mark(ccmark);
      }

      return orientable;
    }

    /** @return true iff the GMap is orientable
     *  If aorientationmark!=INVALID_MARK, and if the cc is orientable, mark
     *       all darts of one orientation of the cc with aorientationmark
     *       (i.e. one out of two darts)
     */
    bool is_orientable(size_type aorientationmark=INVALID_MARK) const
    {
      size_type ccmark=get_new_mark();
      size_type orientationmark=
        (aorientationmark==INVALID_MARK?get_new_mark():aorientationmark);

      bool stop_if_nonorientable = (aorientationmark==INVALID_MARK);
      bool orientable=true;

      for (typename Dart_const_range::const_iterator it=darts().begin(),
             itend=darts().end();
           (!stop_if_nonorientable || orientable) && it!=itend; ++it)
      {
        if (!is_marked(it, ccmark))
        {
          if (!is_cc_orientable(it, ccmark, orientationmark))
            orientable=false;
        }
      }

      for (typename Dart_const_range::const_iterator it=darts().begin(),
             itend=darts().end();
           number_of_marked_darts(ccmark)>0 && it!=itend; ++it)
      {
        unmark(it, ccmark);
        if (aorientationmark==INVALID_MARK) unmark(it, orientationmark);
      }

      if (aorientationmark==INVALID_MARK)
      {
        assert(is_whole_map_unmarked(orientationmark));
        free_mark(orientationmark);
      }

      assert(is_whole_map_unmarked(ccmark));
      free_mark(ccmark);

      return orientable;
    }

    /** Test if the attributes of this map are automatically updated.
     * @return true iff the boolean automatic_attributes_management is set to true.
     */
    bool are_attributes_automatically_managed() const
    {
      return automatic_attributes_management;
    }

    /** Sets the automatic_attributes_management boolean.
     */
    void set_automatic_attributes_management(bool newval)
    {
      if (this->automatic_attributes_management == false && newval == true)
      {
        correct_invalid_attributes();
      }

      this->automatic_attributes_management = newval;
    }

    /** Create an half-edge.
     * @return a dart of the new half-edge.
     */
    Dart_handle make_half_edge()
    {
      Dart_handle d1=create_dart();
      basic_link_alpha<0>(d1, create_dart());
      return d1;
    }

    /** Create an edge.
     * if closed==true, the edge has no 2-free dart.
     * @return a dart of the new edge.
     */
    Dart_handle make_edge(bool closed=false)
    {
      Dart_handle d1=create_dart();
      Dart_handle d2=create_dart();
      basic_link_alpha<0>(d1, d2);
      if (closed)
      {
        Dart_handle d3=create_dart();
        Dart_handle d4=create_dart();
        basic_link_alpha<0>(d3, d4);
        basic_link_alpha<2>(d1, d3);
        basic_link_alpha<2>(d2, d4);
      }

      return d1;
    }

    /** Create an edge given 2 Attribute_handle<0>.
     * Note that this function can be used only if 0-attributes are non void
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * if closed==true, the edge has no 2-free dart.
     * @return the dart of the new edge incident to h0.
     */
    Dart_handle make_segment(typename Attribute_handle<0>::type h0,
                             typename Attribute_handle<0>::type h1,
                             bool closed=false)
    {
      Dart_handle d1 = this->make_edge(closed);

      set_dart_attribute<0>(d1,h0);
      set_dart_attribute<0>(this->alpha<0>(d1),h1);
      if (closed)
      {
        set_dart_attribute<0>(this->alpha<2>(d1),h0);
        set_dart_attribute<0>(this->alpha<0,2>(d1),h1);
      }

      return d1;
    }

    /** Create a combinatorial polygon of length alg
     * (a cycle of alg edges alpha1 links together).
     * @return a new dart.
     */
    Dart_handle make_combinatorial_polygon(unsigned int alg)
    {
      CGAL_assertion(alg>0);

      Dart_handle start = make_edge();
      Dart_handle prev = alpha<0>(start);
      for ( unsigned int nb=1; nb<alg; ++nb )
      {
        Dart_handle cur = make_edge();
        basic_link_alpha<1>(prev, cur);
        prev=alpha<0>(cur);
      }

      basic_link_alpha<1>(prev, start);
      return start;
    }

    /** Test if a face is a combinatorial polygon of length alg
     *  (a cycle of alg edges alpha1 links together).
     * @param adart an intial dart
     * @return true iff the face containing adart is a polygon of length alg.
     */
    bool is_face_combinatorial_polygon(Dart_const_handle adart,
                                       unsigned int alg) const
    {
      CGAL_assertion(alg>0);

      unsigned int nb = 0;
      Dart_const_handle cur = adart;
      do
      {
        ++nb;
        if ( is_free(cur, 0) || is_free(alpha(cur, 0), 1) )
          return false; // Open face
        cur = alpha(cur,0,1);
      }
      while( cur!=adart );
      return (nb==alg);
    }

    /** Create a triangle given 3 Attribute_handle<0>.
     * @param h0 the first handle.
     * @param h1 the second handle.
     * @param h2 the third handle.
     * Note that this function can be used only if 0-attributes are non void
     * @return the dart of the new triangle incident to h0 and to edge h0h1.
     */
    Dart_handle make_triangle(typename Attribute_handle<0>::type h0,
                              typename Attribute_handle<0>::type h1,
                              typename Attribute_handle<0>::type h2)
    {
      Dart_handle d1 = this->make_combinatorial_polygon(3);

      set_dart_attribute<0>(d1,h0);
      set_dart_attribute<0>(this->alpha<1>(d1),h0);
      set_dart_attribute<0>(this->alpha<0>(d1),h1);
      set_dart_attribute<0>(this->alpha<0,1>(d1),h1);
      set_dart_attribute<0>(this->alpha<1,0>(d1),h2);
      set_dart_attribute<0>(this->alpha<1,0,1>(d1),h2);

      return d1;
    }

    /** Create a quadrangle given 4 Vertex_attribute_handle.
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * @param h2 the third vertex handle.
     * @param h3 the fourth vertex handle.
     * Note that this function can be used only if 0-attributes are non void
     * @return the dart of the new quadrilateral incident to h0 and to edge h0h1.
     */
    Dart_handle make_quadrangle(typename Attribute_handle<0>::type h0,
                                typename Attribute_handle<0>::type h1,
                                typename Attribute_handle<0>::type h2,
                                typename Attribute_handle<0>::type h3)
    {
      Dart_handle d1 = this->make_combinatorial_polygon(4);

      set_dart_attribute<0>(d1,h0);
      set_dart_attribute<0>(this->alpha<1>(d1),h0);
      set_dart_attribute<0>(this->alpha<0>(d1),h1);
      set_dart_attribute<0>(this->alpha<0,1>(d1),h1);
      set_dart_attribute<0>(this->alpha<0,1,0>(d1),h2);
      set_dart_attribute<0>(this->alpha<0,1,0,1>(d1),h2);
      set_dart_attribute<0>(this->alpha<1,0>(d1),h3);
      set_dart_attribute<0>(this->alpha<1,0,1>(d1),h3);

      return d1;
    }

    /** Create a combinatorial tetrahedron from 4 triangles.
     * @param d1 a dart onto a first triangle.
     * @param d2 a dart onto a second triangle.
     * @param d3 a dart onto a third triangle.
     * @param d4 a dart onto a fourth triangle.
     * @return d1.
     */
    Dart_handle make_combinatorial_tetrahedron(Dart_handle d1,
                                               Dart_handle d2,
                                               Dart_handle d3,
                                               Dart_handle d4)
    {
      topo_sew<2>(d1, alpha(d2, 0));
      topo_sew<2>(d3, alpha(d2, 1));
      topo_sew<2>(alpha(d1, 0, 1), alpha(d3, 1));
      topo_sew<2>(alpha(d4, 0), alpha(d2, 0, 1));
      topo_sew<2>(alpha(d4, 1), alpha(d3, 0, 1));
      topo_sew<2>(alpha(d4, 0, 1), alpha(d1, 1));

      return d1;
    }

    /** Test if a volume is a combinatorial tetrahedron.
     * @param adart an intial dart
     * @return true iff the volume containing adart is a combinatorial tetrahedron.
     */
    bool is_volume_combinatorial_tetrahedron(Dart_const_handle d1)
    {
      Dart_const_handle d2 = alpha(d1, 0, 2);
      Dart_const_handle d3 = alpha(d2, 1, 2);
      Dart_const_handle d4 = alpha(d2, 0, 1, 0, 2);

      if ( !is_face_combinatorial_polygon(d1, 3) ||
           !is_face_combinatorial_polygon(d2, 3) ||
           !is_face_combinatorial_polygon(d3, 3) ||
           !is_face_combinatorial_polygon(d4, 3) ) return false;

      // TODO do better with marks (?).
      if ( belong_to_same_cell<Self,2,1>(*this, d1, d2) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d3) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d4) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d3) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d4) ||
           belong_to_same_cell<Self,2,1>(*this, d3, d4) ) return false;

      if ( alpha(d1, 0,1,2)!=alpha(d3,1) ||
           alpha(d4, 1,2)!=alpha(d3,0,1) ||
           alpha(d4, 0,1,2)!=alpha(d1,1) ) return false;

      return true;
    }

    /** Create a new combinatorial tetrahedron.
     * @return a new dart.
     */
    Dart_handle make_combinatorial_tetrahedron()
    {
      Dart_handle d1 = make_combinatorial_polygon(3);
      Dart_handle d2 = make_combinatorial_polygon(3);
      Dart_handle d3 = make_combinatorial_polygon(3);
      Dart_handle d4 = make_combinatorial_polygon(3);

      return make_combinatorial_tetrahedron(d1, d2, d3, d4);
    }

    /** Create a combinatorial hexahedron from 6 quadrilaterals.
     * @param d1 a dart onto a first quadrilateral.
     * @param d2 a dart onto a second quadrilateral.
     * @param d3 a dart onto a third quadrilateral.
     * @param d4 a dart onto a fourth quadrilateral.
     * @param d5 a dart onto a fifth quadrilateral.
     * @param d6 a dart onto a sixth quadrilateral.
     * @return d1.
     */
    Dart_handle make_combinatorial_hexahedron(Dart_handle d1,
                                              Dart_handle d2,
                                              Dart_handle d3,
                                              Dart_handle d4,
                                              Dart_handle d5,
                                              Dart_handle d6)
    {
      topo_sew<2>(d1, alpha(d4,1,0,1));
      topo_sew<2>(alpha(d1,0,1), alpha(d6,1));
      topo_sew<2>(alpha(d1,1,0,1), d2);
      topo_sew<2>(alpha(d1,1), d5);

      topo_sew<2>(d3, alpha(d2,1,0,1));
      topo_sew<2>(alpha(d3,0,1,0), alpha(d6,0,1));
      topo_sew<2>(alpha(d3,1,0,1), d4);
      topo_sew<2>(alpha(d3,1,0), alpha(d5,1,0,1));

      topo_sew<2>(d6, alpha(d4,0,1,0));
      topo_sew<2>(alpha(d6,1,0,1), alpha(d2,0,1));

      topo_sew<2>(alpha(d5,1,0), alpha(d4,1));
      topo_sew<2>(alpha(d5,0,1), alpha(d2,1));
      return d1;
    }

    /** Test if a volume is a combinatorial hexahedron.
     * @param adart an intial dart
     * @return true iff the volume containing adart is a combinatorial hexahedron.
     */
    bool is_volume_combinatorial_hexahedron(Dart_const_handle d1)
    {
      Dart_const_handle d2 = alpha(d1,1,0,1,2);
      Dart_const_handle d3 = alpha(d2,1,0,1,2);
      Dart_const_handle d4 = alpha(d3,1,0,1,2);
      Dart_const_handle d5 = alpha(d1,1,2);
      Dart_const_handle d6 = alpha(d4,0,1,0,2);

      if (!is_face_combinatorial_polygon(d1, 4) ||
          !is_face_combinatorial_polygon(d2, 4) ||
          !is_face_combinatorial_polygon(d3, 4) ||
          !is_face_combinatorial_polygon(d4, 4) ||
          !is_face_combinatorial_polygon(d5, 4) ||
          !is_face_combinatorial_polygon(d6, 4) ) return false;

      // TODO do better with marks.
      if ( belong_to_same_cell<Self,2,1>(*this, d1, d2) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d3) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d4) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d5) ||
           belong_to_same_cell<Self,2,1>(*this, d1, d6) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d3) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d4) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d5) ||
           belong_to_same_cell<Self,2,1>(*this, d2, d6) ||
           belong_to_same_cell<Self,2,1>(*this, d3, d4) ||
           belong_to_same_cell<Self,2,1>(*this, d3, d5) ||
           belong_to_same_cell<Self,2,1>(*this, d3, d6) ||
           belong_to_same_cell<Self,2,1>(*this, d4, d5) ||
           belong_to_same_cell<Self,2,1>(*this, d4, d6) ||
           belong_to_same_cell<Self,2,1>(*this, d5, d6) )
        return false;

      if ( alpha(d1,2)      !=alpha(d4,1,0,1) ||
           alpha(d1,0,1,2)  !=alpha(d6,1)     ||
           alpha(d3,0,1,0,2)!=alpha(d6,0,1)   ||
           alpha(d3,1,0,2)  !=alpha(d5,1,0,1) ||
           alpha(d6,1,0,1,2)!=alpha(d2,0,1)   ||
           alpha(d5,1,0,2)  !=alpha(d4,1)     ||
           alpha(d5,0,1,2)  !=alpha(d2,1) ) return false;

      return true;
    }

    /** Create a new combinatorial hexahedron.
     * @return a new dart.
     */
    Dart_handle make_combinatorial_hexahedron()
    {
      Dart_handle d1 = make_combinatorial_polygon(4);
      Dart_handle d2 = make_combinatorial_polygon(4);
      Dart_handle d3 = make_combinatorial_polygon(4);
      Dart_handle d4 = make_combinatorial_polygon(4);
      Dart_handle d5 = make_combinatorial_polygon(4);
      Dart_handle d6 = make_combinatorial_polygon(4);

      return make_combinatorial_hexahedron(d1, d2, d3, d4, d5, d6);
    }

    /** Test if an i-cell can be removed.
     *  An i-cell can be removed if i==dimension or i==dimension-1,
     *     or if there are at most two (i+1)-cell incident to it.
     * @param adart a dart of the i-cell.
     * @return true iff the i-cell can be removed.
     */
    template < unsigned int i >
    bool is_removable(Dart_const_handle adart) const
    { return CGAL::Is_removable_functor_gmap<Self, i>::run(*this, adart); }

    /** Remove an i-cell, 0<=i<=dimension.
     * @param adart a dart of the i-cell to remove.
     * @param update_attributes a boolean to update the enabled attributes
     * @return the number of deleted darts.
     */
    template < unsigned int i >
    size_t remove_cell(Dart_handle adart, bool update_attributes = true)
    {
      return CGAL::Remove_cell_functor_gmap<Self,i,Self::dimension-i>::
        run(*this,adart,update_attributes);
    }

    /** Test if an i-cell can be contracted.
     *  An i-cell can be contracted if i==1
     *     or if there are at most two (i-1)-cell incident to it.
     * @param adart a dart of the i-cell.
     * @return true iff the i-cell can be contracted.
     */
    template < unsigned int i >
    bool is_contractible(Dart_const_handle adart) const
    { return CGAL::Is_contractible_functor_gmap<Self, i>::run(*this,adart); }

    /** Contract an i-cell, 1<=i<=dimension.
     * @param adart a dart of the i-cell to remove.
     * @return the number of deleted darts.
     */
    template < unsigned int i >
    size_t contract_cell(Dart_handle adart, bool update_attributes = true)
    {
      return CGAL::Contract_cell_functor_gmap<Self,i>::
        run(*this, adart, update_attributes);
    }

    /** Insert a vertex in a given edge.
     * @param adart a dart of the edge (!=NULL).
     * @return a dart of the new vertex.
     */
    Dart_handle insert_cell_0_in_cell_1( Dart_handle adart,
                                         typename Attribute_handle<0>::type
                                         ah=null_handle,
                                         bool update_attributes=true )
    {
      Dart_handle d1;
      size_type amark=get_new_mark();

      // 1) We store all the darts of the edge.
      std::deque<Dart_handle> vect;
      size_type m=get_new_mark();
      {
        for ( typename Dart_of_cell_basic_range<1>::iterator
                it=darts_of_cell_basic<1>(adart, m).begin();
              it != darts_of_cell_basic<1>(adart, m).end(); ++it )
          vect.push_back(it);
      }

      // 2) For each dart of the cell, we modify link of neighbors.
      typename std::deque<Dart_handle>::iterator it = vect.begin();
      for (; it != vect.end(); ++it)
      {
        d1 = create_dart();

        if (!(this->template is_free<0>(*it)) &&
            is_marked(alpha<0>(*it), amark))
          basic_link_alpha<1>(d1, alpha<0,0>(*it));

        basic_link_alpha<0>(*it, d1);
        mark(*it, amark);

        for ( unsigned int dim=2; dim<=dimension; ++dim )
        {
          if (!is_free(*it, dim) && is_marked(alpha(*it, dim), amark))
          {
            basic_link_alpha(alpha(*it, dim, 0), d1, dim);
          }
        }

        if (are_attributes_automatically_managed() && update_attributes)
        {
          // We copy all the attributes except for dim=0
          Helper::template Foreach_enabled_attributes_except
            <CGAL::internal::GMap_group_attribute_functor_of_dart<Self>, 0>::
            run(*this,*it,d1);
        }
        if (ah != null_handle)
        {
          // We initialise the 0-atttrib to ah
          CGAL::internal::Set_i_attribute_of_dart_functor<Self, 0>::
            run(*this, d1, ah);
          mark(*it, amark);
        }
      }

      for (it = vect.begin(); it != vect.end(); ++it)
      {
        unmark(*it, m);
        unmark(*it, amark);
      }

      CGAL_assertion(is_whole_map_unmarked(m));
      CGAL_assertion(is_whole_map_unmarked(amark));

      free_mark(m);
      free_mark(amark);

      if (are_attributes_automatically_managed() && update_attributes)
      {
        if ( !(this->template is_free<1>(alpha<0>(adart))) )
        {
          CGAL::internal::GMap_degroup_attribute_functor_run<Self, 1>::
              run(*this, adart, alpha<0,1>(adart));
        }
      }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
      CGAL_assertion( is_valid() );
#endif

      return alpha<0, 1>(adart);
    }

    /** Insert a vertex in the given 2-cell which is splitted in triangles,
     * once for each inital edge of the facet.
     * @param adart a dart of the facet to triangulate.
     * @return A dart incident to the new vertex.
     */
    Dart_handle insert_cell_0_in_cell_2( Dart_handle adart,
                                         typename Attribute_handle<0>::type
                                         ah=null_handle,
                                         bool update_attributes=true )
    {
      CGAL_assertion(adart!=null_handle);

      Dart_handle d1=null_handle, d2=null_handle;

      // Mark used to mark darts already treated.
      size_type treated = get_new_mark();
      size_type m = get_new_mark();
      size_type edge_pushed = get_new_mark();

      // Stack of darts of the face
      std::deque<Dart_handle> vect;
      {
        for ( typename Dart_of_cell_basic_range<2>::iterator
                it=darts_of_cell_basic<2>(adart,m).begin();
              it!=darts_of_cell_basic<2>(adart,m).end(); ++it )
          vect.push_back(it);
      }

      // Stack of darts to degroup (one dart per edge of the face)
      std::deque<Dart_handle> todegroup;
      if (are_attributes_automatically_managed() && update_attributes)
      {
        for ( typename Dart_of_cell_basic_range<2,2>::iterator
              it=darts_of_cell_basic<2,2>(adart).begin();
              it!=darts_of_cell_basic<2,2>(adart).end(); ++it )
          if ( it!=adart && it!=alpha<0>(adart) && !is_marked(it, edge_pushed))
          {
            todegroup.push_back(it);
            mark(it, edge_pushed);
            mark(this->template alpha<0>(it), edge_pushed);
          }
      }

    // 2) For each dart of the cell, we modify link of neighbors.
    typename std::deque<Dart_handle>::iterator it = vect.begin();
    for (; it != vect.end(); ++it)
    {
      d1 = create_dart();
      d2 = create_dart();
      basic_link_alpha<0>(d1, d2);
      mark(*it, treated);

      if (!(this->template is_free<0>(*it)) &&
          is_marked(this->template alpha<0>(*it), treated))
        basic_link_alpha<1>(d2, this->template alpha<0,1,0>(*it));

      if (!(this->template is_free<1>(*it)) &&
          is_marked(this->template alpha<1>(*it), treated))
      {
        basic_link_alpha<2>(d1, this->template alpha<1,1>(*it));
        basic_link_alpha<2>(d2, this->template alpha<1,1,0>(*it));
      }

      basic_link_alpha<1>(*it, d1);

      for ( unsigned int dim=3; dim<=dimension; ++dim )
      {
        if (!is_free(*it, dim) && is_marked(alpha(*it, dim), treated))
        {
          basic_link_alpha(alpha(*it, dim, 1), d1, dim);
          basic_link_alpha(alpha(*it, dim, 1, 0), d2, dim);
        }
      }

      if (are_attributes_automatically_managed() && update_attributes)
      {
        // We copy all the attributes except for dim=1
        Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::GMap_group_attribute_functor_of_dart<Self>, 1>::
          run(*this,*it,d1);
        Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::GMap_group_attribute_functor_of_dart<Self>, 0>::
          run(*this,d1,d2);

        // We initialise the 0-atttrib to ah
        CGAL::internal::Set_i_attribute_of_dart_functor<Self, 0>::
          run(*this, d2, ah);
      }
    }

    for (it = vect.begin(); it != vect.end(); ++it)
    {
      unmark(*it, m);
      unmark(*it, treated);
      unmark(*it, edge_pushed);
    }

    CGAL_assertion(is_whole_map_unmarked(m));
    CGAL_assertion(is_whole_map_unmarked(treated));
    CGAL_assertion(is_whole_map_unmarked(edge_pushed));
    free_mark(m);
    free_mark(treated);
    free_mark(edge_pushed);

    if (are_attributes_automatically_managed() && update_attributes)
    {
      for (it = todegroup.begin(); it != todegroup.end(); ++it)
      {
        CGAL::internal::GMap_degroup_attribute_functor_run<Self, 2>::
          run(*this, adart, *it);
      }
    }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
    CGAL_assertion( is_valid() );
#endif

    return alpha<1,0>(adart);
  }

  /** Test if an edge can be inserted onto a 2-cell between two given darts.
   * @param adart1 a first dart.
   * @param adart2 a second dart.
   * @return true iff an edge can be inserted between adart1 and adart2.
   */
  bool is_insertable_cell_1_in_cell_2(Dart_const_handle adart1,
                                      Dart_const_handle adart2)
  {
    if ( adart1==adart2 || adart1==this->template alpha<0>(adart2) )
      return false;
    for ( CGAL::GMap_dart_const_iterator_of_orbit<Self,0,1> it(*this,adart1);
          it.cont(); ++it )
    {
      if ( it==adart2 ) return true;
    }
    return false;
  }

  /** Insert an edge in a 2-cell between two given darts.
   * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
   * @param adart2 a second dart of the facet. If NULL insert a dangling edge.
   * @return a dart of the new edge, and not incident to the
   *         same vertex than adart1.
   */
  Dart_handle insert_cell_1_in_cell_2(Dart_handle adart1,
                                      Dart_handle adart2,
                                      bool update_attributes=true,
                                      typename Attribute_handle<0>::type
                                      ah=null_handle)
  {
    if ( adart2!=null_handle)
    {
      CGAL_assertion(is_insertable_cell_1_in_cell_2(adart1, adart2));
    }

    /* CGAL::GMap_dart_iterator_basic_of_involution<Self,1> will contain all
     * alpha_i except alpha_0, alpha_1 and alpha_2, i.e. this is
     * <alpha_3,...,alpha_d>
     */
    size_type m1=get_new_mark();
    CGAL::GMap_dart_iterator_basic_of_involution<Self,1> it1(*this, adart1, m1);

    size_type m2=get_new_mark();
    CGAL::GMap_dart_iterator_basic_of_involution<Self,1> it2(*this, adart2, m2);

    Dart_handle d1=null_handle;
    Dart_handle d2=null_handle;
    Dart_handle d3=null_handle;
    Dart_handle d4=null_handle;

    size_type treated=get_new_mark();
    bool isfree1 = (this->template is_free<1>(adart1));

    for ( ; it1.cont(); ++it1)
    {
      d1 = create_dart();
      d2 = create_dart();
      mark(it1,treated);

      if (!isfree1)
      {
        d3 = create_dart();
        d4 = create_dart();
        this->template basic_link_alpha<2>(d1, d3);
        this->template basic_link_alpha<2>(d2, d4);
      }

      for ( unsigned int dim=3; dim<=dimension; ++dim)
      {
        if ( !is_free(it1, dim) &&
             is_marked(alpha(it1, dim), treated) )
        {
          basic_link_alpha(alpha(it1, dim, 1), d1, dim);
          basic_link_alpha(alpha(d1, dim, 0), d2, dim);

          if (!isfree1)
          {
            basic_link_alpha(alpha(it1, 1, dim, 1), d3, dim);
            basic_link_alpha(alpha(d3, dim, 0), d4, dim);
          }
        }
      }

      if (!isfree1)
      {
        this->template link_alpha<1>(this->template alpha<1>(it1), d3);
        if ( adart2!=null_handle )
        {
          CGAL_assertion (it2.cont());
          this->template link_alpha<1>(this->template alpha<1>(it2), d4);
        }
      }

      this->template link_alpha<1>(it1, d1);
      if ( adart2!=null_handle )
      {
        CGAL_assertion (it2.cont());
        this->template link_alpha<1>(it2, d2);
        ++it2;
      }
      else
      {
        if (are_attributes_automatically_managed() &&
            update_attributes && ah!=NULL)
        {
          internal::Set_i_attribute_of_dart_functor<Self, 0>::run(*this, d2, ah);
          if (!isfree1)
          {
            internal::Set_i_attribute_of_dart_functor<Self, 0>::run(*this, d4, ah);
          }
        }
      }

      // We do the link_alpha<0> after the link_alpha<1> to update the
      // possible attributes of d2.
      this->template link_alpha<0>(d1, d2);
      if (!isfree1)
      { this->template link_alpha<0>(d3, d4); }
    }

    if (are_attributes_automatically_managed() && update_attributes)
    {
      if ( !this->template is_free<2>(d1) && d2!=null_handle )
        CGAL::internal::GMap_degroup_attribute_functor_run<Self, 2>::
          run(*this, d1, this->template alpha<2>(d1));
    }

    negate_mark(m1);
    it1.rewind();

    if ( adart2!=null_handle )
    { it2.rewind(); negate_mark(m2); }

    for (; it1.cont(); ++it1)
    {
      mark(it1,m1);
      unmark(it1,treated);
      if ( adart2!=null_handle )
      { mark(it2,m2); ++it2; }
    }
    negate_mark(m1);
    CGAL_assertion( is_whole_map_unmarked(m1) );
    CGAL_assertion( is_whole_map_unmarked(treated) );
    free_mark(m1);
    free_mark(treated);

    if ( adart2!=null_handle )
    {
      negate_mark(m2);
      CGAL_assertion( is_whole_map_unmarked(m2) );
    }
    free_mark(m2);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
    CGAL_assertion( is_valid() );
#endif

    return this->template alpha<1,0>(adart1);
  }

  /** Insert a dangling edge in a 2-cell between given by a dart.
   * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
   * @param update_attributes a boolean to update the enabled attributes
   * @return a dart of the new edge, not incident to the vertex of adart1.
   */
  Dart_handle insert_dangling_cell_1_in_cell_2( Dart_handle adart1,
                                                typename Attribute_handle<0>::
                                                type ah=null_handle,
                                                bool update_attributes=true )
  { return insert_cell_1_in_cell_2(adart1, NULL, update_attributes, ah); }

  /** Test if a 2-cell can be inserted onto a given 3-cell along
   *  a path of edges.
   *  @param afirst iterator on the begining of the path.
   *  @param alast  iterator on the end of the path.
   *  @return true iff a 2-cell can be inserted along the path.
   *          the path is a sequence of dartd, one per edge
   *          where the face will be inserted.
   */
  template <class InputIterator>
  bool is_insertable_cell_2_in_cell_3(InputIterator afirst,
                                      InputIterator alast)
  {
    CGAL_assertion( dimension>= 3 );

    // The path must have at least one dart.
    if (afirst==alast) return false;
    Dart_const_handle prec = null_handle;

    for (InputIterator it(afirst); it!=alast; ++it)
    {
      // The path must contain only non empty darts.
      if (*it == null_handle) return false;

      if (this->template is_free<0>(*it)) return false;

      // Two consecutive darts of the path must belong to two edges
      // incident to the same vertex of the same volume.
      if (prec!=null_handle)
      {
        // prec and *it must belong to the same vertex of the same volume
        if ( !CGAL::belong_to_same_cell<Self, 0, 2>(*this,  prec, *it) )
          return false;
      }
      prec = this->template alpha<0>(*it);
    }

    // The path must be closed.
    if (!CGAL::belong_to_same_cell<Self, 0, 2>(*this, prec, *afirst))
      return false;

    return true;
  }

  /** Insert a 2-cell in a given 3-cell along a path of darts.
   * @param amap the used generalized map.
   * @param afirst iterator on the begining of the path.
   * @param alast  iterator on the end of the path.
   *         the path is a sequence of darts, one per edge
   *         where the face will be inserted.
   * @return a dart of the new 2-cell.
   */
  template<class InputIterator>
  Dart_handle insert_cell_2_in_cell_3(InputIterator afirst,
                                      InputIterator alast,
                                      bool update_attributes=true)
  {
    CGAL_assertion(is_insertable_cell_2_in_cell_3(afirst,alast));

    Dart_handle prec = null_handle, d = null_handle,
        dd = null_handle, first = null_handle, ddd = null_handle,
        dddd = null_handle, it0, d0, oldb2;

    bool withAlpha3 = false;

    {
      for (InputIterator it(afirst); !withAlpha3 && it!=alast; ++it)
      {
        if (!this->template is_free<2>(*it)) withAlpha3 = true;
      }
    }

    {
      for (InputIterator it(afirst); it!=alast; ++it)
      {
        d = create_dart();
        d0 = create_dart();
        basic_link_alpha<0>(d, d0);

        Helper::template Foreach_enabled_attributes_except
          <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
          run(*this, d,*it);

        it0=alpha<0>(*it);
        Helper::template Foreach_enabled_attributes_except
          <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
          run(*this, d0, it0);

        if ( withAlpha3 )
        {
          dd = create_dart();
          d0 = create_dart();
          basic_link_alpha<0>(dd, d0);

          basic_link_alpha<3>(d, dd);
          basic_link_alpha<3>(alpha<0>(d), alpha<0>(dd));

          Helper::template Foreach_enabled_attributes_except
            <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
            run(*this,dd,d);

          Helper::template Foreach_enabled_attributes_except
            <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
            run(*this,d0, it0);
        }

        if ( prec!=null_handle )
        {
          basic_link_alpha<1>(prec, d);
          if (withAlpha3)
            basic_link_alpha<1>(alpha<3>(prec), dd);
        }
        else first=d;

        if ( !this->template is_free<2>(*it) )
        {
          oldb2=alpha<2>(*it);
          basic_link_alpha<2>(oldb2, dd);
          basic_link_alpha<2>(alpha<0>(oldb2), alpha<0>(dd));
        }
        else
          oldb2=NULL;

        basic_link_alpha<2>(*it, d);
        basic_link_alpha<2>(alpha<0>(*it), alpha<0>(d));

        // Make copies of the new facet for dimension >=4
        for ( unsigned int dim=4; dim<=dimension; ++dim )
        {
          if ( !is_free(*it, dim) )
          {
            ddd = create_dart();
            d0 = create_dart();
            basic_link_alpha<0>(ddd, d0);

            basic_link_alpha(d, ddd, dim);
            basic_link_alpha(alpha<0>(d), d0, dim);

            basic_link_alpha<2>(alpha(*it, dim), ddd);
            basic_link_alpha<2>(alpha(*it, 0, dim), d0);

            Helper::template Foreach_enabled_attributes_except
              <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
              run(*this, ddd, *it);

            it0=alpha<0>(*it);
            Helper::template Foreach_enabled_attributes_except
              <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
              run(*this, d0, it0);

            if ( withAlpha3 )
            {
              dddd = create_dart();
              d0 = create_dart();
              basic_link_alpha<0>(dddd, d0);

              basic_link_alpha(dd, dddd, dim);
              basic_link_alpha(alpha<0>(dd), d0, dim);

              if (oldb2!=NULL)
              {
                basic_link_alpha<2>(alpha(oldb2, dim), dddd);
                basic_link_alpha<2>(alpha(oldb2, 0, dim), d0);
              }

              Helper::template Foreach_enabled_attributes_except
                <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
                run(*this, dddd, *it);
              Helper::template Foreach_enabled_attributes_except
                <internal::GMap_group_attribute_functor_of_dart<Self, 2>, 2>::
                run(*this, d0, it0);
            }

            if ( prec!=null_handle )
            {
              basic_link_alpha<1>(alpha(prec, dim), ddd);
              if (withAlpha3)
                basic_link_alpha<1>(alpha(prec, 3, dim), dddd);
            }

          }
        }

        prec = alpha<0>(d);
      }
    }

    basic_link_alpha<1>(prec, first);
    if ( withAlpha3 )
    {
      basic_link_alpha<1>(alpha<3>(prec), alpha<3>(first));
    }

    for ( unsigned int dim=4; dim<=dimension; ++dim )
    {
      if ( !is_free(first, dim) )
      {
        basic_link_alpha<1>(alpha(prec, dim), alpha(first, dim));
        if ( withAlpha3 )
        {
          basic_link_alpha<1>(alpha(prec, 3, dim), alpha(first, 3, dim));
        }
      }
    }

    // Degroup the attributes
    if ( withAlpha3 )
    { // Here we cannot use Degroup_attribute_functor_run as new darts do not
      // have their 3-attribute
      if (are_attributes_automatically_managed() && update_attributes)
      {
        CGAL::internal::GMap_degroup_attribute_functor_run<Self, 3>::
            run(*this, first, alpha<3>(first));
      }
    }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
    CGAL_assertion( is_valid() );
#endif

    return first;
  }

  protected:
    /// Number of times each mark is reserved. 0 if the mark is free.
    mutable size_type mnb_times_reserved_marks[NB_MARKS];

    /// Mask marks to know the value of unmark dart, for each index i.
    mutable std::bitset<NB_MARKS> mmask_marks;

    /// Number of used marks.
    mutable size_type mnb_used_marks;

    /// Index of each mark, in mfree_marks_stack or in mfree_marks_stack.
    mutable size_type mindex_marks[NB_MARKS];

    /// "Stack" of free marks.
    mutable size_type mfree_marks_stack[NB_MARKS];

    /// "Stack" of used marks.
    mutable size_type mused_marks_stack[NB_MARKS];

    /// Number of marked darts for each used marks.
    mutable size_type mnb_marked_darts[NB_MARKS];

    /// Automatic management of the attributes:
    /// true means attributes are always maintained updated during operations.
    bool automatic_attributes_management;

    /// Tuple of unary and binary functors (for all non void attributes).
    typename Helper::Split_functors m_onsplit_functors;
    typename Helper::Merge_functors m_onmerge_functors;
  };

  template < unsigned int d_,
             class Items_=Generic_map_min_items,
             class Alloc_=CGAL_ALLOCATOR(int),
             class Storage_= Generalized_map_storage_1<d_, Items_, Alloc_> >
  class Generalized_map :
    public Generalized_map_base<d_,
                                  Generalized_map<d_,Items_,Alloc_, Storage_>,
                                  Items_, Alloc_, Storage_ >
  {
  public:
    typedef Generalized_map<d_, Items_,Alloc_, Storage_>  Self;
    typedef Generalized_map_base<d_, Self, Items_, Alloc_, Storage_> Base;

    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;
    typedef typename Base::Alloc Alloc;
    typedef typename Base::Exception_no_more_available_mark
    Exception_no_more_available_mark;

    Generalized_map() : Base()
    {}

    Generalized_map(const Self & amap) : Base(amap)
    {}

    template < class Gmap >
    Generalized_map(const Gmap & amap) : Base(amap)
    {}

    template < class Gmap, typename Converters >
    Generalized_map(const Gmap & amap, const Converters& converters) :
      Base(amap, converters)
    {}

    template < class Gmap, typename Converters, typename DartInfoConverter >
    Generalized_map(const Gmap & amap, const Converters& converters,
                    const DartInfoConverter& dartinfoconverter) :
      Base(amap, converters, dartinfoconverter)
    {}

    template < class Gmap, typename Converters, typename DartInfoConverter,
               typename PointConverter >
    Generalized_map(const Gmap & amap, const Converters& converters,
                    const DartInfoConverter& dartinfoconverter,
                    const PointConverter& pointconverter) :
      Base(amap, converters, dartinfoconverter, pointconverter)
    {}
  };

} // namespace CGAL

#if  (BOOST_GCC >= 40900)
 _Pragma("GCC diagnostic pop")
#endif

#endif // CGAL_GENERALIZED_MAP_H //
// EOF //
