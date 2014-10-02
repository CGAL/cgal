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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_STORAGES_H
#define CGAL_COMBINATORIAL_MAP_STORAGES_H 1

#include <CGAL/Handle_hash_function.h>

#include <CGAL/Compact_container_with_index_2.h>

namespace CGAL {

  /** @file Combinatorial_map_storages.h
   * Definition of storages for dD Combinatorial map.
   */

  struct Index_hash_function {
    typedef std::size_t result_type;
    template <class H>
    std::size_t operator() (const H& h) const {
      return h;
    }
  };

  // Storage of darts with compact container, beta with handles
  template<unsigned int d_, class Items_, class Alloc_ >
  class Combinatorial_map_storage_1
  {
  public:
    typedef Combinatorial_map_storage_1<d_, Items_, Alloc_> Self;
    typedef CGAL::Tag_false Use_index;

    typedef internal::Combinatorial_map_helper<Self> Helper;

    typedef typename Items_::template Dart_wrapper<Self>  Dart_wrapper;
    typedef typename Dart_wrapper::Dart                   Dart;
    typedef typename Alloc_::template rebind<Dart>::other Dart_allocator;

    typedef Compact_container<Dart,Dart_allocator>  Dart_container;

    typedef typename Dart_container::iterator       Dart_handle;
    typedef typename Dart_container::const_iterator Dart_const_handle;
    typedef typename Dart_container::size_type      size_type;

    typedef CGAL::Void* Null_handle_type;
    static Null_handle_type null_handle;

    typedef Items_ Items;
    typedef Alloc_ Alloc;

    template <typename T>
    struct Container_for_attributes :
        public Compact_container<T, typename Alloc_::template rebind<T>::other>
    {};

    /// Typedef for attributes
    typedef typename Dart_wrapper::Attributes Attributes;

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

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d_;

    typedef Handle_hash_function Hash_function;

    // Init
    void init_storage()
    {
#ifdef CGAL_CMAP_DEPRECATED
      // We must do this ony once, but problem because null_dart_handle
      // is static !
      if ( mnull_dart_container.empty() )
#endif // CGAL_CMAP_DEPRECATED
      { // emplace null_dart; initialized in Combinatorial_map class
        null_dart_handle = mnull_dart_container.emplace();
      }
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
      return dh->mbeta[i]==null_dart_handle;
    }
    bool is_free(Dart_const_handle dh, unsigned int i) const
    {
      CGAL_assertion( dh!=NULL );
      CGAL_assertion(i <= dimension);
      return dh->mbeta[i]==null_dart_handle;
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
    bool get_dart_mark(Dart_const_handle ADart, int amark) const
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_handle ADart, int amark, bool avalue) const
    {
      CGAL_assertion( ADart!=NULL );
      ADart->set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_handle ADart, int amark) const
    {
      CGAL_assertion( ADart!=NULL );
      ADart->flip_mark(amark);
    }

    // Access to beta maps
    Dart_handle get_beta(Dart_handle ADart, int B1)
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return ADart->mbeta[B1];
    }
    Dart_const_handle get_beta(Dart_const_handle ADart, int B1) const
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mbeta[B1];
    }
    template<int B1>
    Dart_handle get_beta(Dart_handle ADart)
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mbeta[B1];
    }
    template<int B1>
    Dart_const_handle get_beta(Dart_const_handle ADart) const
    {
      CGAL_assertion(ADart!=NULL && B1>=0 && B1<=(int)dimension);
      return  ADart->mbeta[B1];
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

    Dart & get_dart(Dart_handle ah)
    {
      CGAL_assertion( ah!=NULL );
      return *ah;
    }
    const Dart & get_dart(Dart_const_handle ah) const
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
    Dart_handle & dart(Dart_handle adart)
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

    void display_dart(Dart_const_handle ADart) const
    { std::cout<<&*ADart; }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_handle<i>::type ah) const
    { std::cout<<&*ah; }

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
      adart->mbeta[i] = adart2;
    }
    void dart_link_beta(Dart_handle adart, Dart_handle adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=NULL && adart2!=NULL);
      CGAL_assertion(adart!=null_dart_handle);
      adart->mbeta[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_beta(Dart_handle adart)
    {
      CGAL_assertion(adart!=NULL && i <= dimension);
      adart->mbeta[i] = null_dart_handle;
    }
    void dart_unlink_beta(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(adart!=NULL && i <= dimension);
      adart->mbeta[i] = null_dart_handle;
    }

  public:
    /// Void dart. A dart d is i-free if beta_i(d)=null_dart_handle.
#ifdef CGAL_CMAP_DEPRECATED
    static
#endif // CGAL_CMAP_DEPRECATED
    Dart_handle null_dart_handle; // Todo Dart_const_handle ??

  protected:
    /// Dart container.
    Dart_container mdarts;

    /// Container for the null_dart_handle, static data member.
#ifdef CGAL_CMAP_DEPRECATED
    static
#endif // CGAL_CMAP_DEPRECATED
    Dart_container mnull_dart_container;

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

  /// null_handle
  template < unsigned int d_, class Items_, class Alloc_ >
  typename Combinatorial_map_storage_1<d_, Items_, Alloc_>::Null_handle_type
  Combinatorial_map_storage_1<d_, Items_, Alloc_>::null_handle = NULL;

#ifdef CGAL_CMAP_DEPRECATED
  /// Allocation of static data members
  /// mnull_dart_container
  template<unsigned int d_, class Items_, class Alloc_ >
  typename Combinatorial_map_storage_1<d_, Items_, Alloc_>::Dart_container
  Combinatorial_map_storage_1<d_, Items_, Alloc_>::mnull_dart_container;

  /// null_dart_handle
  template < unsigned int d_, class Items_, class Alloc_ >
  typename Combinatorial_map_storage_1<d_, Items_, Alloc_>::Dart_handle
  Combinatorial_map_storage_1<d_, Items_, Alloc_>::null_dart_handle;
  // =  mnull_dart_container.emplace( std::bitset<NB_MARKS>() );
  // Does not work on windows => segfault
  // Thus we initialize null_dart_handle in the Combinatorial_map constructor
#endif // CGAL_CMAP_DEPRECATED

template<typename T>
class MyIndex
{
public:
 MyIndex(size_t s=-1) : idx(s)
 {}
 operator size_t() const
 { return idx; }
 MyIndex<T>& operator++()
 { ++idx; return *this; }
 MyIndex<T> operator++(int)
 { MyIndex<T> res(*this); ++idx; return res; }
 MyIndex<T>& operator--()
 { --idx; return *this; }
 MyIndex<T> operator--(int)
 {MyIndex<T> res(*this); --idx; return res;}
private:
 T idx;
};

  // Storage with combinatorial maps using index
  template<unsigned int d_, class Items_, class Alloc_>
  class Combinatorial_map_storage_2
  {
  public:
    typedef Combinatorial_map_storage_2<d_, Items_, Alloc_> Self;

    typedef internal::Combinatorial_map_helper<Self> Helper;

    typedef typename Items_::template Dart_wrapper<Self>  Dart_wrapper;
    typedef typename Dart_wrapper::Dart                   Dart;
    typedef typename Alloc_::template rebind<Dart>::other Dart_allocator;

    typedef Compact_container_with_index_2<Dart,Dart_allocator,
    Constant_size_policy_for_cc_with_size<1024> >
    Dart_container;

    typedef CGAL::Tag_true Use_index;
    typedef typename Dart_container::size_type      size_type;

    typedef Items_ Items;
    typedef Alloc_ Alloc;

    template <typename T>
    struct Container_for_attributes : public
        Compact_container_with_index_2<T,
        typename Alloc_::template rebind<T>::other,
        Constant_size_policy_for_cc_with_size<1024> >
    {};

    /// Typedef for attributes
    typedef typename Dart_wrapper::Attributes Attributes;

    template<int i>
    struct Attribute_type: public Helper::template Attribute_type<i>
    {};
    template<int i>
    struct Attribute_handle: public CGAL::MyIndex<unsigned int> // public Helper::template Attribute_handle<i>
    {
      Attribute_handle(size_t s=-1) : MyIndex<unsigned int>(s)
      {}
    };
    template<int i>
    struct Attribute_const_handle:public CGAL::MyIndex<unsigned int>//Helper::template Attribute_const_handle<i>
    {
      Attribute_const_handle(size_t s=-1) : MyIndex<unsigned int>(s)
      {}
    };
    template<int i>
    struct Attribute_range: public Helper::template Attribute_range<i>
    {};
    template<int i>
    struct Attribute_const_range:
      public Helper::template Attribute_const_range<i>
    {};

    /// Number of marks
    static const size_type NB_MARKS = 32;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d_;

    // typedef unsigned int Dart_index;
    // typedef MyIndex<unsigned int> Dart_index;
    struct Dart_index : public CGAL::MyIndex<unsigned int>
    {
      Dart_index(size_t s=-1) : MyIndex<unsigned int>(s)
      {}
    };

    // Definition of old types, for backward compatibility.
    typedef Dart_index Dart_handle;
    typedef Dart_index Dart_const_handle;

    /// Value of null handle (!= null_dart_handle !!)
    typedef Dart_index Null_handle_type;
    static const Dart_index null_handle; //=0;

    typedef Index_hash_function Hash_function;

    void init_storage()
    {}

    /** Return if this dart is free for adimension.
     * @param dh a dart handle
     * @param i the dimension.
     * @return true iff dh is linked with NULL for \em adimension.
     */
    template<unsigned int i>
    bool is_free(Dart_const_handle dh) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mbeta[i]==null_dart_handle;
    }
    bool is_free(Dart_const_handle dh, unsigned int i) const
    {
      CGAL_assertion(i <= dimension);
      return mdarts[dh].mbeta[i]==null_dart_handle;
    }

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
    bool get_dart_mark(Dart_const_handle ADart, int amark) const
    {
      return mdarts[ADart].get_mark(amark);
    }

    /// Set the mark of a given mark number to a given value.
    void set_dart_mark(Dart_const_handle ADart, int amark, bool avalue) const
    {
      mdarts[ADart].set_mark(amark, avalue);
    }

    /// Flip the mark of a given mark number to a given value.
    void flip_dart_mark(Dart_const_handle ADart, int amark) const
    {
      mdarts[ADart].flip_mark(amark);
    }

    // Access to beta maps
    Dart_handle get_beta(Dart_handle ADart, int B1)
    {
      CGAL_assertion(B1>=0 && B1<=dimension);
      return mdarts[ADart].mbeta[B1];
    }
    Dart_const_handle get_beta(Dart_const_handle ADart, int B1) const
    {
      CGAL_assertion(B1>=0 && B1<=dimension);
      return  mdarts[ADart].mbeta[B1];
    }
    template<int B1>
    Dart_handle get_beta(Dart_handle ADart)
    {
      CGAL_assertion(B1>=0 && B1<=dimension);
      return  mdarts[ADart].mbeta[B1];
    }
    template<int B1>
    Dart_const_handle get_beta(Dart_const_handle ADart) const
    {
      CGAL_assertion(B1>=0 && B1<=dimension);
      return  mdarts[ADart].mbeta[B1];
    }

    // return a handle on the i-attribute
    template<unsigned int i>
    typename Attribute_handle<i>::type attribute(Dart_handle ADart)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_handles);
    }
    template<unsigned int i>
    typename Attribute_const_handle<i>::type
    attribute(Dart_const_handle ADart) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                     "attribute<i> called but i-attributes are disabled.");
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_handles);
    }

    template<unsigned int i>
    typename Attribute_type<i>::type&
    get_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }

    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah];
    }

    Dart & get_dart(Dart_handle ah)
    {
      CGAL_assertion( ah!=null_handle );
      return mdarts[ah];
    }
    const Dart & get_dart(Dart_const_handle ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return mdarts[ah];
    }

    // Get the attribute of a dart
    template<unsigned int i>
    typename Attribute_type<i>::type& get_attribute_of_dart(Dart_handle adart)
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( attribute<i>(adart)!=null_handle );
      return get_attribute<i>(attribute<i>(adart));
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type&
    get_attribute_of_dart(Dart_const_handle adart) const
    {
      CGAL_assertion( adart!=null_handle );
      CGAL_assertion( attribute<i>(adart)!=null_handle );
      return get_attribute<i>(attribute<i>(adart));
    }

    // Get the dart of the given attribute
    template<unsigned int i>
    Dart_handle dart_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah].dart();
    }
    template<unsigned int i>
    Dart_const_handle
    dart_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah].dart();
    }

    // Set the dart of the given attribute
    template<unsigned int i>
    void set_dart_of_attribute(typename Attribute_handle<i>::type ah,
                               Dart_handle adart)
    {
      CGAL_assertion( ah!=null_handle );
       CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah].set_dart(adart);
    }

    // Get the info of the given attribute
    template<unsigned int i>
    typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_handle<i>::type ah)
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah].info();
    }
    template<unsigned int i>
    const typename Attribute_type<i>::type::Info &
    info_of_attribute(typename Attribute_const_handle<i>::type ah) const
    {
      CGAL_assertion( ah!=null_handle );
      return CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers)[ah].info();
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
    Dart_handle & dart(Dart_handle adart)
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

    void display_dart(Dart_const_handle ADart) const
    { std::cout<<ADart; }

    template<unsigned int i>
    void display_attribute(typename Attribute_const_handle<i>::type ah) const
    { std::cout<<ah; }

  protected:
    // Set the handle on the i th attribute
    template<unsigned int i>
    void basic_set_dart_attribute(Dart_handle ADart,
                                  typename Attribute_handle<i>::type ah)
    {
      CGAL::cpp11::get<Helper::template Dimension_index<i>::value>
        (mdarts[ADart].mattribute_handles) = ah;
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
      CGAL_assertion(adart!=null_dart_handle);
      mdarts[adart].mbeta[i] = adart2;
    }
    void dart_link_beta(Dart_handle adart, Dart_handle adart2, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      CGAL_assertion(adart!=null_dart_handle);
      mdarts[adart].mbeta[i] = adart2;
    }

    /** Unlink a dart for a given dimension.
     * @param adart a dart.
     * @param i the dimension.
     */
    template<unsigned int i>
    void dart_unlink_beta(Dart_handle adart)
    {
      CGAL_assertion(i <= dimension);
      mdarts[adart].mbeta[i] = null_dart_handle;
    }
    void dart_unlink_beta(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(i <= dimension);
      mdarts[adart].mbeta[i] = null_dart_handle;
    }

  public:
    /// Void dart. A dart d is i-free if beta_i(d)=null_dart_handle.
    static const Dart_index null_dart_handle; //=0;

  protected:
    /// Dart container.
    Dart_container mdarts;

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

  /// null_dart_handle
  template<unsigned int d_, class Items_, class Alloc_>
  const typename Combinatorial_map_storage_2<d_, Items_, Alloc_>::Dart_index
  Combinatorial_map_storage_2<d_, Items_, Alloc_>::null_dart_handle(0);


} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_H //
// EOF //
