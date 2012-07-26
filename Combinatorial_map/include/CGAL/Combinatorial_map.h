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
#ifndef CGAL_COMBINATORIAL_MAP_H
#define CGAL_COMBINATORIAL_MAP_H 1

#include <CGAL/Compact_container.h>
#include <CGAL/internal/Combinatorial_map_utility.h>
#include <CGAL/internal/Combinatorial_map_functors.h>
#include <CGAL/Combinatorial_map_min_items.h>
#include <CGAL/Dart_const_iterators.h>
#include <CGAL/Cell_const_iterators.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <bitset>
#include <vector>

// suppress bogus warning when compiling with gcc 4.3 or 4.4
#if (__GNUC__ == 4 && (__GNUC_MINOR__ == 3 || __GNUC_MINOR__ == 4))
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif 

namespace CGAL {
  
  /** @file Combinatorial_map.h
   * Definition of generic dD Combinatorial map.
   */

  /** Generic definition of combinatorial map in dD.
   * The Combinatorial_map class describes an dD combinatorial map. It allows
   * mainly to create darts, to use marks onto these darts, to get and set
   * the beta links, and to manage enabled attributes.
   */
  template < unsigned int d_, class Refs,
             class Items_=Combinatorial_map_min_items<d_>,
             class Alloc_=CGAL_ALLOCATOR(int) >
  class Combinatorial_map_base
  {    
    template<class Map, unsigned int i, unsigned int nmi>
    friend struct Remove_cell_functor;    

    template<class Map>
    friend typename Map::Dart_handle 
    insert_cell_0_in_cell_1(Map& amap, typename Map::Dart_handle adart);
    
    template<class CMap>
    friend typename CMap::Dart_handle 
    insert_cell_1_in_cell_2(CMap& amap, typename CMap::Dart_handle adart1);

    template<class CMap>
    friend
    typename CMap::Dart_handle 
    insert_cell_1_in_cell_2(CMap& amap,
                            typename CMap::Dart_handle adart1,
                            typename CMap::Dart_handle adart2);
    
    template<class Map, class InputIterator>
    friend
    typename Map::Dart_handle
    insert_cell_2_in_cell_3(Map& amap, InputIterator afirst, 
                            InputIterator alast);

    template < class Map >
    friend
    typename Map::Dart_handle 
    insert_cell_0_in_cell_2(Map& amap, typename Map::Dart_handle adart);

    template<typename CMap, unsigned int i, typename Type_attr, typename Range>
    friend struct internal::Degroup_one_attribute_of_dart_functor;

    template <typename CMap, unsigned int i, typename Type_attr>
    friend struct internal::Degroup_one_attribute_functor;

    template<typename Map>
    friend struct internal::Test_is_valid_attribute_functor;

    template<typename CMap, unsigned int i>
    friend struct internal::Group_attribute_functor_of_dart_run;

    template<typename Map,unsigned int i>
    friend struct internal::Group_attribute_functor_run;

    template <typename CMap, unsigned int i, typename Type_attr>
    friend struct internal::Group_one_attribute_functor;

    template<typename Map,unsigned int i>
    friend struct internal::Degroup_attribute_functor_run;

    template<typename Map>
    friend struct internal::Update_dart_of_attribute_functor;    

    template<typename Map>
    friend struct internal::Decrease_attribute_functor;
    
  public:
    /// Types definition
    typedef Combinatorial_map_base<d_, Refs, Items_,Alloc_>  Self;

    typedef Items_ Items;
    typedef Alloc_ Alloc;

    typedef typename Items::template Dart_wrapper<Refs> Dart_wrapper;
    typedef typename Dart_wrapper::Dart                 Dart;

    typedef typename Alloc::template rebind<Dart>::other Dart_allocator;
    typedef Compact_container<Dart,Dart_allocator>       Dart_container;

    typedef typename Dart_container::iterator       Dart_handle;
    typedef typename Dart_container::const_iterator Dart_const_handle;
    typedef typename Dart_container::size_type      size_type;

    /// The dimension of the combinatorial map.
    static const unsigned int dimension = d_;

    /// Number of marks
    static const size_type NB_MARKS = 32;

    typedef internal::Combinatorial_map_helper<Self> Helper;

    /// Typedef for Dart_range, a range through all the darts of the map.
    typedef Dart_container       Dart_range;
    typedef const Dart_container Dart_const_range;

    typedef typename Dart_wrapper::Attributes Attributes;

    /// Typedef for attributes
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
    
  public:
    /** Default Combinatorial_map constructor.
     * The map is empty.
     */
    Combinatorial_map_base()
    {
      CGAL_static_assertion_msg(Dart::dimension==dimension,
                  "Dimension of dart different from dimension of map");

      CGAL_static_assertion_msg(Helper::nb_attribs<=dimension+1,
                  "Too many attributes in the tuple Attributes_enabled");

      this->mnb_used_marks = 0;
      this->mmask_marks.reset();

      for (size_type i = 0; i < NB_MARKS; ++i)
      {
        this->mfree_marks_stack[i]        = (int)i;
        this->mindex_marks[i]             = i;
        this->mnb_marked_darts[i]         = 0;
        this->mnb_times_reserved_marks[i] = 0;
      }

      // We must do this ony once, but problem because null_dart_handle 
      // is static !
      if ( mnull_dart_container.empty() )
      {
        null_dart_handle =
          mnull_dart_container.emplace( std::bitset<NB_MARKS>() );
          
        for (unsigned int i=0; i<=dimension; ++i)
        {
          null_dart_handle->unlink_beta(i);
        }
      }

      CGAL_assertion(number_of_darts()==0);
    }

    /** Clear the combinatorial map. Remove all darts and all attributes.
     *  Note that reserved marks are not free.
     */
    void clear()
    {
      mdarts.clear();
      for (unsigned int i = 0; i < NB_MARKS; ++i)
        this->mnb_marked_darts[i]  = 0;

      internal::Clear_all::run(mattribute_containers);
    }

    /** Test if the map is empty.
     *  @return true iff the map is empty.
     */
    bool is_empty() const
    { return mdarts.empty(); }  

    /** Create a new dart and add it to the map.
     * The marks of the darts are initialised with mmask_marks, i.e. the dart
     * is unmarked for all the marks.
     * @return a Dart_handle on the new dart.
     */
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template < typename... Args >
    Dart_handle create_dart(const Args&... args)
    { return mdarts.emplace(mmask_marks, args...); }
#else
    Dart_handle create_dart()
    { return mdarts.emplace(mmask_marks); }
    template < typename T1 >
    Dart_handle create_dart(const T1 &t1)
    { return mdarts.emplace(mmask_marks, t1); }
    template < typename T1, typename T2 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2)
    { return mdarts.emplace(mmask_marks, t1, t2); }
    template < typename T1, typename T2, typename T3 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3)
    { return mdarts.emplace(mmask_marks, t1, t2, t3); }
    template < typename T1, typename T2, typename T3, typename T4 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
    { return mdarts.emplace(mmask_marks, t1, t2, t3, t4); }
    template < typename T1, typename T2, typename T3, typename T4, typename T5 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                            const T5 &t5)
    { return mdarts.emplace(mmask_marks, t1, t2, t3, t4, t5); }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                            const T5 &t5, const T6 &t6)
    { return mdarts.emplace(mmask_marks, t1, t2, t3, t4, t5, t6); }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6, typename T7 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                            const T5 &t5, const T6 &t6, const T7 &t7)
    { return mdarts.emplace(mmask_marks, t1, t2, t3, t4, t5, t6, t7); }
    template < typename T1, typename T2, typename T3, typename T4, typename T5,
               typename T6, typename T7, typename T8 >
    Dart_handle create_dart(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4,
                            const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
    { return mdarts.emplace(mmask_marks, t1, t2, t3, t4, t5, t6, t7, t8); }
#endif

    /** Erase a dart from the list of darts.
     * @param adart the dart to erase.
     */
    void erase_dart(Dart_handle adart)
    {
      // 1) We update the number of marked darts.
      for (unsigned int i = 0; i < mnb_used_marks; ++i)
      {
        if (is_marked(adart, mused_marks_stack[i]))
          --mnb_marked_darts[mused_marks_stack[i]];
      }

      // 2) We update the attribute_ref_counting.
      Helper::template Foreach_enabled_attributes
        <internal::Decrease_attribute_functor<Self> >::run(this,adart);

      // 3) We erase the dart.
      mdarts.erase(adart);
    }

    /// @return a Dart_range (range through all the darts of the map).
    Dart_range& darts()             { return mdarts;}
    Dart_const_range& darts() const { return mdarts; }

    /** Get the first dart of this map.
     * @return the first dart.
     */
    Dart_handle first_dart()
    {
      if (darts().begin() == darts().end()) return null_dart_handle;
      return mdarts.begin();
    }
    Dart_const_handle first_dart() const
    {
      if (darts().begin() == darts().end()) return null_dart_handle;
      return mdarts.begin();
    }

    /// @return the Dart_handle corresponding to the given dart.
    Dart_handle dart_handle(Dart& adart)
    { return mdarts.iterator_to(adart); }
    Dart_const_handle dart_handle(const Dart& adart) const
    { return mdarts.iterator_to(adart); }

    /// @return the betas of ADart (beta are used in the same order than
    ///         they are given as parameters)
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template<typename ... Betas>
    Dart_handle beta(Dart_handle ADart, Betas... betas) const
    { return internal::Beta_functor<Dart_handle, Betas ...>::
        run(ADart, betas...); }
    template<typename ... Betas>
    Dart_const_handle beta(Dart_const_handle ADart, Betas... betas) const
    { return internal::Beta_functor<Dart_const_handle, Betas ...>::
        run(ADart, betas...); }
#else
    Dart_handle beta(Dart_handle ADart, int B1)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4,
                                                      B5); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5,
                     int B6)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4, B5,
                                                      B6); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5,
                     int B6, int B7)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4, B5,
                                                      B6, B7); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5,
                     int B6, int B7, int B8)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4, B5,
                                                      B6, B7, B8); }
    Dart_handle beta(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5,
                     int B6, int B7, int B8, int B9)
    { return internal::Beta_functor<Dart_handle>::run(ADart, B1, B2, B3, B4, B5,
                                                      B6, B7, B8, B9); }
#endif
    
    /** Count the number of used marks.
     * @return the number of used marks.
     */
    size_type number_of_used_marks() const
    { return mnb_used_marks; }

    /** Test if a given mark is reserved.
     *  @return true iff the mark is reserved (ie in used).
     */ 
    bool is_reserved(int amark) const
    { 
      CGAL_assertion(amark>=0 && (size_type)amark<NB_MARKS);
      return (mnb_times_reserved_marks[(size_type)amark]!=0);
    }

    /**  Count the number of marked darts for a given mark.
     * @param amark the mark index.
     * @return the number of marked darts for amark.
     */
    size_type number_of_marked_darts(int amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      return mnb_marked_darts[(size_type)amark];
    }

    /**  Count the number of unmarked darts for a given mark.
     * @param amark the mark index.
     * @return the number of unmarked darts for amark.
     */
    size_type number_of_unmarked_darts(int amark) const
    { 
      CGAL_assertion( is_reserved(amark) );
      return number_of_darts() - number_of_marked_darts(amark); 
    }

    /** Test if all the darts are unmarked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are unmarked for amark.
     */
    bool is_whole_map_unmarked(int amark) const
    { return number_of_marked_darts(amark) == 0; }

    /** Test if all the darts are marked for a given mark.
     * @param amark the mark index.
     * @return true iff all the darts are marked for amark.
     */
    bool is_whole_map_marked(int amark) const
    {  return number_of_marked_darts(amark) == number_of_darts(); }

    /** Reserve a new mark.
     * Get a new free mark and return its index.
     * All the darts are unmarked for this mark.
     * @return the index of the new mark.
     * @pre mnb_used_marks < NB_MARKS
     */
    int get_new_mark() const
    {
      if (mnb_used_marks == NB_MARKS)
      {
        std::cerr << "Not enough Boolean marks: "
          "increase NB_MARKS in item class." << std::endl;
        return -1;
      }

      int m = mfree_marks_stack[mnb_used_marks];
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
    void share_a_mark(int amark) const
    {
      CGAL_assertion( is_reserved(amark) );
      ++mnb_times_reserved_marks[amark];
    }

    /** @return the number of times a mark is reserved.
     *  @param amark the mark to share.
     */
    size_type get_number_of_times_mark_reserved(int amark) const
    { return mnb_times_reserved_marks[amark]; }
    
    /** Negate the mark of all the darts for a given mark.
     * After this call, all the marked darts become unmarked and all the
     * unmarked darts become marked (in constant time operation).
     * @param amark the mark index
     */
    void negate_mark(int amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      mnb_marked_darts[amark] =
        number_of_darts() - mnb_marked_darts[amark];

      mmask_marks.flip((size_type)amark);
    }

    /** Test if a given dart is marked for a given mark.
     * @param adart the dart to test.
     * @param amark the given mark.
     * @return true iff adart is marked for the mark amark.
     */
    bool is_marked(Dart_const_handle adart, int amark) const
    {
      CGAL_assertion( adart != null_dart_handle );
      CGAL_assertion( is_reserved(amark) );

      return adart->get_mark(amark)!=mmask_marks[(size_type)amark];
    }

    /** Set the mark of a given dart to a state (on or off).
     * @param adart the dart.
     * @param amark the given mark.
     * @param astate the state of the mark (on or off).
     */
    void set_mark_to(Dart_const_handle adart, int amark, 
                     bool astate) const
    {
      CGAL_assertion( adart != null_dart_handle );
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark) != astate)
      {
        if (astate) ++mnb_marked_darts[(size_type)amark];
        else --mnb_marked_darts[(size_type)amark];

        adart->set_mark(amark, astate ^ mmask_marks[(size_type)amark]);
      }
    }

    /** Mark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void mark(Dart_const_handle adart, int amark) const
    {
      CGAL_assertion( adart != null_dart_handle );
      CGAL_assertion( is_reserved(amark) );

      if (is_marked(adart, amark)) return;

      ++mnb_marked_darts[(size_type)amark];
      adart->set_mark(amark, !mmask_marks[(size_type)amark]);
    }

    /** Unmark the given dart.
     * @param adart the dart.
     * @param amark the given mark.
     */
    void unmark(Dart_const_handle adart, int amark) const
    {
      CGAL_assertion( adart != null_dart_handle );
      CGAL_assertion( is_reserved(amark) );

      if (!is_marked(adart, amark)) return;

      --mnb_marked_darts[(size_type)amark];
      adart->set_mark(amark, mmask_marks[(size_type)amark]);      
    }

    /** Unmark all the darts of the map for a given mark.
     * If all the darts are marked or unmarked, this operation takes O(1)
     * operations, otherwise it traverses all the darts of the map.
     * @param amark the given mark.
     */
    void unmark_all(int amark) const
    {
      CGAL_assertion( is_reserved(amark) );

      if (is_whole_map_unmarked(amark)) return;

      if (is_whole_map_marked(amark))
      {
        negate_mark(amark);
      }
      else
      {
        for (typename Dart_range::const_iterator it(darts().begin()), 
               itend(darts().end()); it!=itend; ++it)
          unmark(it, amark);
      }
      CGAL_assertion(is_whole_map_unmarked(amark));
    }

    /** Free a given mark, previously calling unmark_all_darts.
     * @param amark the given mark.
     */
    void free_mark(int amark) const
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
      mused_marks_stack[mindex_marks[(size_type)amark]] =
        mused_marks_stack[--mnb_used_marks];
      mindex_marks[mused_marks_stack[mnb_used_marks]] =
        mindex_marks[(size_type)amark];

      // 2) We add amark in the array mfree_marks_stack and update its index.
      mfree_marks_stack[ mnb_used_marks ] = amark;
      mindex_marks[(size_type)amark] = mnb_used_marks;

      mnb_times_reserved_marks[amark]=0;
    }

    /** Test if this map is without boundary for a given dimension.
     * @param i the dimension.
     * @return true iff all the darts are not i-free.
     * @pre 1<=i<=n
     */
    bool is_without_boundary(unsigned int i) const
    {
      CGAL_assertion(1<=i && i<=dimension);
      for (typename Dart_const_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
        if (it->is_free(i)) return false;
      return true;
    }

    /** Test if this map is without boundary for all the dimensions.
     * @return true iff all the darts are non free.
     */
    bool is_without_boundary() const
    {
      for (typename Dart_const_range::const_iterator it(darts().begin()), 
             itend(darts().end()); it!=itend; ++it)
        for (unsigned int i = 1; i<=dimension; ++i)
          if (it->is_free(i)) return false;
      return true;
    }

    /** Close the combinatorial map for a given dimension.
     *  @param i the dimension to close
     *  @return the number of new darts.
     *  @pre 2<=i<=n (TODO case i==1)
     */
    unsigned int close(unsigned int i)
    {
      CGAL_assertion(2<=i && i<=dimension);
      unsigned int res = 0;
      Dart_handle d, d2;

      for (typename Dart_range::iterator it(darts().begin());
           it!=darts().end(); ++it)
      {
        if ( it->is_free(i) )
        {
          d = create_dart();
          ++res;
          link_beta(it, d, i);

          // Special cases for 0 and 1
          if ( !it->is_free(1) && !it->beta(1)->is_free(i) )
            link_beta<1>(it->beta(1)->beta(i),d);
          if ( !it->is_free(0) && !it->beta(0)->is_free(i) )
            link_beta<0>(it->beta(0)->beta(i),d);
          // General case for 2...dimension
          for (unsigned int j=2; j<=dimension; ++j)
          {
            if ( j+1!=i && j!=i && j!=i+1 && 
                 !it->is_free(j) && !it->beta(j)->is_free(i) )
            {
              link_beta(it->beta(j)->beta(i), d, j);
            }
          }

          d2 = it;
          while (d2 != null_dart_handle && !d2->is_free(i-1))
          { d2 = d2->beta(i-1)->beta(i); }
          if (d2 != null_dart_handle) 
          {
            if (i==2) link_beta<1>(d2, d);
            else link_beta(d2, d, i-1);
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
      std::vector<int> marks(dimension+1);
      for (i=0; i<=dimension; ++i)
        marks[i] = -1;

      Helper::template
        Foreach_enabled_attributes<internal::Reserve_mark_functor<Self> >::
        run(this,&marks);

      for (typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        if ( !valid )
        { // We continue the traversal to mark all the darts.              
          for (i=0; i<=dimension; ++i)
            if (marks[i]!=-1) mark(it,marks[i]);
        }
        else
        {
          // beta0 must be the inverse of beta1
          if ((!it->is_free(0) && it->beta(0)->beta(1)!=it) ||
              (!it->is_free(1) && it->beta(1)->beta(0)!=it ))
          {
            std::cerr << "Map not valid: beta(0) "
              "is not the inverse of beta(1) for "
                      <<&(*it) << std::endl;
            valid = false;
          }
              
          // Each beta(i>=2) must be an involution
          for (i = 2; i <= dimension; ++i)
            if (!it->is_free(i) && it->beta(i)->beta(i)!=it)
            {
              std::cerr << "Map not valid: beta(" << i
                        << ") is not an involution for "
                        <<&(*it) << std::endl;
              valid = false;
            }
              
          // beta1 o betai and beta0 o betai (i>=3) must be involutions
          if (!it->is_free(0))
          {
            for (i = 3; i <= dimension; ++i)
              if ((it->is_free(i) != it->beta(0)->is_free(i)) ||
                  (!it->is_free(i) &&
                   it->beta(0)->beta(i)!=it->beta(i)->beta(1)))
              {
                std::cerr << "Map not valid: beta(0) o beta(" << i
                          << ") is not an involution for "
                          <<&(*it) << std::endl;
                valid = false;
              }
          }
          if (!it->is_free(1))
          {
            for (i = 3; i <= dimension; ++i)
              if ((it->is_free(i) != it->beta(1)->is_free(i)) ||
                  (!it->is_free(i) &&
                   it->beta(1)->beta(i)!=it->beta(i)->beta(0)))
              {
                std::cerr << "Map not valid: beta(1) o beta(" << i
                          << ") is not an involution for " 
                          <<&(*it)<< std::endl;
                valid = false;
              }
          }
              
          // beta(i>=2) o beta(j>=i+2) must be an involution
          for (i = 2; i <= dimension; ++i)
          {
            if (!it->is_free(i))
            {
              for (j = i + 2; j <= dimension; ++j)
                if ((it->is_free(j)!=it->beta(i)->is_free(j)) ||
                    (!it->is_free(j) &&
                     it->beta(i)->beta(j)!=it->beta(j)->beta(i)))
                {
                  std::cerr << "Map not valid: beta(" << i
                            << ") o beta(" << j 
                            << ") is not an involution for "
                            << &(*it)<< std::endl;
                  valid = false;
                }
            }
          }
          Helper::template Foreach_enabled_attributes
            <internal::Test_is_valid_attribute_functor<Self> >::
            run(this,it,&marks,&valid);
        }
      }
      for (i=0; i<=dimension; ++i)
        if ( marks[i]!=-1 ) 
        {
          CGAL_assertion( is_whole_map_marked(marks[i]) );
          free_mark(marks[i]);
        }

      return valid;
    }

    /// @return the number of darts.
    size_type number_of_darts() const
    { return mdarts.size(); }

    /// @return an estimation of the bytes used by the combinatorial map.
    size_type bytes() const
    {
      return mdarts.capacity() * sizeof(Dart) +
        internal::Count_bytes_all_attributes_functor<Self>::run(*this);
    }

    /** Write the content of the map: each dart and each beta links.
     * @param os the ostream.
     * @return the ostream.
     */
    std::ostream& display_darts(std::ostream & os) const
    {
      unsigned int nb = 0;
      for (typename Dart_range::const_iterator it=darts().begin();
           it!=darts().end(); ++it)
      {
        os << " dart " << &(*it) << "; beta[i]=";
        for (unsigned int i=0; i<=dimension; ++i)
        {
          os << &(*it->beta(i)) << ",\t"; 
          if (it->is_free(i))os << "\t";
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
      int amark = get_new_mark();
      for (typename Dart_range::const_iterator it1(darts().begin()),
             itend(darts().end()); it1!=itend; ++it1)
      {
        if ( !is_marked(it1, amark) )
        {
          ++nb;
          for ( Ite it2(*this, it1, amark); it2.cont(); ++it2 )
          {
            aos << &(**it2) << " - " << std::flush;
            mark(*it2, amark);
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
      return display_orbits<CMap_dart_const_iterator_basic_of_cell<Self,i> >
        (aos); 
    }

    /** Write the number of darts and cells of the map into a given ostream.
     * @param os the ostream.
     * @return the ostream.
     */
    std::ostream& display_characteristics(std::ostream & os) const
    {
      std::vector<unsigned int> cells(dimension+2);
      for (unsigned int i=0; i<=dimension+1; ++i)
      { cells[i]=i; }
  
      std::vector<unsigned int> res = count_cells(cells);
  
      os << "#Darts=" << number_of_darts();
      for (unsigned int i=0; i<=dimension; ++i)
        os<<", #"<<i<<"-cells="<<res[i];
      os<<", #ccs="<<res[dimension+1];
      
      return os;
    }

    /// Create a new attribute.
    /// @return a handle on the new attribute.
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template<unsigned int i, typename... Args>
    typename Attribute_handle<i>::type create_attribute(const Args&... args)
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(args...); 
    }
#else
    template<unsigned int i>
    typename Attribute_handle<i>::type
    create_attribute()
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(); 
    }
    template<unsigned int i, typename T1>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1)
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1);
    }
    template<unsigned int i, typename T1, typename T2>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2)
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2); 
    }
    template<unsigned int i, typename T1, typename T2, typename T3>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3)
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3); 
    }
    template<unsigned int i, typename T1, typename T2, typename T3, typename T4>
    typename Attribute_handle<i>::type
    create_attribute(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
    {      
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "create_attribute<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
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
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
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
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
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
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
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
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).emplace(t1, t2, t3, t4, t5, t6, t7, t8); 
    }
#endif

    /// Erase an attribute.
    /// @param h a handle to the attribute to erase.
    template<unsigned int i>
    void erase_attribute(typename Attribute_handle<i>::type h)
    { 
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "erase_attribute<i> but i-attributes are disabled");
      CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).erase(h); 
    }

    /// @return the number of attributes.
    template <unsigned int i>
    size_type number_of_attributes() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "number_of_attributes<i> but i-attributes are disabled");
      return  CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers).size(); 
    }

    /// @return a Attributes_range<i> (range through all the 
    /// attributes<i> of the map).
    template<unsigned int i>  
    typename Attribute_range<i>::type & attributes() 
    { 
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "attributes<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers);
    }
    
    template<unsigned int i>  
    typename Attribute_const_range<i>::type & attributes() const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "attributes<i> but i-attributes are disabled");
      return CGAL::cpp0x::get<Helper::template Dimension_index<i>::value>
        (mattribute_containers); 
    }

    /** Double link a dart with beta 0 to a second dart.
     * \em adart1 is 0-linked to \em adart2 and \em adart2 is 1-linked
     * with \em adart1. Attributes are not updated, thus we can obtain
     * a non-valid map with darts belonging to a same orbit and having
     * different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    void basic_link_beta_0(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      adart1->basic_link_beta(adart2, 0);
      adart2->basic_link_beta(adart1, 1);
    }

    /** Double link a dart with beta 0 to a second dart.
     * \em adart1 is 0-linked to \em adart2 and \em adart2 is 1-linked
     * with \em adart1. Attributes are not updated, thus we can obtain
     * a non-valid map with darts belonging to a same orbit and having
     * different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    void basic_link_beta_1(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      adart1->basic_link_beta(adart2, 1);
      adart2->basic_link_beta(adart1, 0);
    }
    
    /** Double link a dart with beta i to a second dart, when i>=2.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i-linked
     * with \em adart1. Attributes are not updated, thus we can obtain
     * a non-valid map with darts belonging to a same orbit and having
     * different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param i the dimension of the beta.
     */
    void basic_link_beta_for_involution(Dart_handle adart1, Dart_handle adart2,
                                        unsigned int i)
    {
      CGAL_assertion( i>=2 && i<=dimension );
      CGAL_assertion(adart1 != NULL && adart2 != NULL && adart1!=adart2);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      adart1->basic_link_beta(adart2, i);
      adart2->basic_link_beta(adart1, i);
    }

    /** Double link a dart with betai to a second dart.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i^-1-linked
     * with \em adart1. Attributes are not updated, thus we can obtain
     * a non-valid map with darts belonging to a same orbit and having
     * different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    template<unsigned int i>
    void basic_link_beta(Dart_handle adart1, Dart_handle adart2)
    {
      if ( i==0 ) basic_link_beta_0(adart1, adart2);
      else if ( i==1 ) basic_link_beta_1(adart1, adart2);
      else basic_link_beta_for_involution(adart1, adart2, i);
    }
    void basic_link_beta(Dart_handle adart1, Dart_handle adart2,
                         unsigned int i)
    {
      if ( i==0 ) basic_link_beta_0(adart1, adart2);
      else if ( i==1 ) basic_link_beta_1(adart1, adart2);
      else basic_link_beta_for_involution(adart1, adart2, i);
    }

    /** Double link two darts, and update the NULL attributes.
     * \em adart1 is 0-linked to \em adart2 and \em adart2 is 1-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa.
     * We can obtain a non-valid map with darts belonging to a same cell
     * and having different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    void link_beta_0(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      Helper::template Foreach_enabled_attributes
        <internal::Group_attribute_functor_of_dart<Self> >::
        run(this,adart1,adart2,0);
      adart1->basic_link_beta(adart2, 0);
      adart2->basic_link_beta(adart1, 1);
    }

    /** Double link two darts, and update the NULL attributes.
     * \em adart1 is 1-linked to \em adart2 and \em adart2 is 0-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa.
     * We can obtain a non-valid map with darts belonging to a same cell
     * and having different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    void link_beta_1(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      Helper::template Foreach_enabled_attributes
        <internal::Group_attribute_functor_of_dart<Self> >::
        run(this,adart1,adart2,1);
      adart1->basic_link_beta(adart2, 1);
      adart2->basic_link_beta(adart1, 0);
    }

    /** Double link two darts, and update the NULL attributes.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i^-1-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa.
     * We can obtain a non-valid map with darts belonging to a same cell
     * and having different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param i the dimension of the beta.
     * @pre 2<=i<=dimension.
     */
    void link_beta_for_involution(Dart_handle adart1, Dart_handle adart2,
                                  unsigned int i)
    {
      CGAL_assertion(adart1 != NULL && adart2 != NULL && adart1!=adart2 );
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);
      CGAL_assertion( 2<=i && i<=dimension );
      Helper::template Foreach_enabled_attributes
        <internal::Group_attribute_functor_of_dart<Self> >::
        run(this,adart1,adart2,i);
      adart1->basic_link_beta(adart2, i);
      adart2->basic_link_beta(adart1, i);
    }

    /** Double link two darts, and update the NULL attributes.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i^-1-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa.
     * We can obtain a non-valid map with darts belonging to a same cell
     * and having different attributes.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     */
    template<unsigned int i>
    void link_beta(Dart_handle adart1, Dart_handle adart2)
    {
      if ( i==0 ) link_beta_0(adart1, adart2);
      else if ( i==1 ) link_beta_1(adart1, adart2);
      else link_beta_for_involution(adart1, adart2, i);
    }
    void link_beta(Dart_handle adart1, Dart_handle adart2, unsigned int i)
    {
      if ( i==0 ) link_beta_0(adart1, adart2);
      else if ( i==1 ) link_beta_1(adart1, adart2);
      else link_beta_for_involution(adart1, adart2, i);
    }

    /** Double link a dart with betai to a second dart.
     * \em adart1 is i-linked to \em adart2 and \em adart2 is i^-1-linked
     * with \em adart1. The NULL attributes of \em adart1 are updated to
     * non NULL attributes associated to \em adart2, and vice-versa, only .
     * if update_attributes==true.
     * @param adart1 a first dart.
     * @param adart2 a second dart.
     * @param update_attributes a boolean to update the enabled attributes.
     */
    template<unsigned int i>
    void link_beta(Dart_handle adart1, Dart_handle adart2, 
                   bool update_attributes)
    {
      if ( update_attributes ) link_beta<i>(adart1, adart2);
      else basic_link_beta<i>(adart1, adart2);
    }

    /** Double unlink a dart with beta 0.
     * beta0(\em adart) is 1-unlinked and \em adart is 0-unlinked.
     * The attributes are not updated, thus we can obtain a non-valid map
     * with darts belonging to different orbits and having the same
     * attributes.
     * @param adart a dart.
     */
    void unlink_beta_0(Dart_handle adart)
    {
      CGAL_assertion(adart != NULL && !adart->is_free(0));
      adart->beta(0)->unlink_beta(1);
      adart->unlink_beta(0);
    }

    /** Double unlink a dart with beta 1.
     * beta1(\em adart) is 0-unlinked and \em adart is 1-unlinked.
     * The attributes are not updated, thus we can obtain a non-valid map
     * with darts belonging to different orbits and having the same
     * attributes.
     * @param adart a dart.
     */
    void unlink_beta_1(Dart_handle adart)
    {
      CGAL_assertion(adart != NULL && !adart->is_free(1));
      adart->beta(1)->unlink_beta(0);
      adart->unlink_beta(1);
    }

    /** Double unlink a dart with beta i, for i>=2.
     * betai(\em adart) is i-unlinked and \em adart is i-unlinked.
     * The attributes are not updated, thus we can obtain a non-valid map
     * with darts belonging to different orbits and having the same
     * attributes.
     * @param adart a dart.
     * @param i the dimension of the beta.
     */
    void unlink_beta_for_involution(Dart_handle adart, unsigned int i)
    {
      CGAL_assertion(adart!=NULL && adart!=null_dart_handle && 
                     !adart->is_free(i));
      CGAL_assertion(2<=i && i<=dimension);
      adart->beta(i)->unlink_beta(i);
      adart->unlink_beta(i);
    }

    /** Double unlink a dart with beta i.
     * betai(\em adart) is i-1-unlinked and \em adart is i-unlinked.
     * The attributes are not updated, thus we can obtain a non-valid map
     * with darts belonging to different orbits and having the same
     * attributes.
     * @param adart a dart.
     * @param i the dimension of the beta.
     */
    template<unsigned int i>
    void unlink_beta(Dart_handle adart)
    {
      if ( i==0 ) unlink_beta_0(adart);
      else if ( i==1 ) unlink_beta_1(adart);
      else unlink_beta_for_involution(adart, i);
    }
    void unlink_beta(Dart_handle adart, unsigned int i)
    {
      if ( i==0 ) unlink_beta_0(adart);
      else if ( i==1 ) unlink_beta_1(adart);
      else unlink_beta_for_involution(adart, i);
    }

    /** Test if it is possible to sew by beta1 the two given darts 
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff \em adart1 can be 1-sewn with \em adart2.
     */
    bool is_sewable_1(Dart_const_handle adart1, Dart_const_handle adart2) const
    {
      CGAL_assertion(adart1!=NULL && adart2!=NULL);
      
      if ( !adart1->is_free(1) || !adart2->is_free(0) )
        return false;
      
      if ( adart1 == adart2 ) return true;
      
      CMap_dart_const_iterator_of_involution    <Self,1> I1(*this, adart1);
      CMap_dart_const_iterator_of_involution_inv<Self,1> I2(*this, adart2);
      bool res = true;
      while (res && I1.cont() && I2.cont())
      {
        // We can remove this constraint which is not required for 
        // combinatorial map definition, but which imposes quite "normal"
        // configurations
        if ( I1==adart2 || I2==adart1 ) res=false;
        
        for (unsigned int j=3;res && j<=Self::dimension; ++j)
        {
          if ( I1->is_free(j)!=I2->is_free(j) )
          {
            res = false;
          }
        }
        ++I1; ++I2;
      }
      if (I1.cont() != I2.cont()) 
        res = false;
      
      return res;
    }
    
    /** Test if it is possible to sew by beta0 the two given darts 
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff \em adart1 can be 0-sewn with \em adart2.
     */
    bool is_sewable_0(Dart_const_handle adart1, Dart_const_handle adart2) const
    { return is_sewable_1(adart2, adart1); }

    /** Test if it is possible to sew by betai the two given darts
     * for 2<=i<=dimension.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff \em adart1 can be 1-sewn with \em adart2.
     */
    template<unsigned int i>
    bool is_sewable_for_involution(Dart_const_handle adart1,
                                   Dart_const_handle adart2) const
    {
      CGAL_assertion(2<=i && i<=Self::dimension);
      CGAL_assertion(adart1!=NULL && adart2!=NULL);
      
      if ( !adart1->is_free(i) || !adart2->is_free(i) || adart1==adart2 )
        return false;
      
      CMap_dart_const_iterator_of_involution<Self,i>     I1(*this, adart1);
      CMap_dart_const_iterator_of_involution_inv<Self,i> I2(*this, adart2);
      bool res = true;
      while (res && I1.cont() && I2.cont())
      {
        // We can remove this constraint which is not required for 
        // combinatorial map definition, but which is quite "normal"
        if ( I1==adart2 || I2==adart1 ) res=false;
        
        // Special case to consider beta0 and beta1
        if ( i>2 )
        {
          if ( I1->is_free(0)!=I2->is_free(1) )      res = false;
          else if ( I1->is_free(1)!=I2->is_free(0) ) res = false;
        }
        
        // General case
        for (unsigned int j=2;res && j<=Self::dimension; ++j)
        {
          if ( j+1!=i && j!=i && j!=i+1 && 
               I1->is_free(j)!=I2->is_free(j) )
          { res = false; }
        }
        ++I1; ++I2;
      }
      if (I1.cont() != I2.cont()) 
        res = false;
      
      return res;
    }
    
    /** Test if it is possible to sew by betai the two given darts 
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff \em adart1 can be i-sewn with \em adart2.
     */
    template<unsigned int i>
    bool is_sewable(Dart_const_handle adart1, Dart_const_handle adart2) const
    {
      if ( i==0 ) return is_sewable_0(adart1, adart2);
      else if ( i==1 ) return is_sewable_1(adart1, adart2);
      else return is_sewable_for_involution<i>(adart1, adart2);
    }    

    /** Topological sew by beta1 the two given darts plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes 
     * thus the map can be non valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable_1(adart1, adart2).
     */
    void topo_sew_1(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( (is_sewable_1(adart1,adart2)) );
      
      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart1,m);
           it.cont(); ++it)
      {
        mark(it,m);
        dartv.push_back(it);
      }

      CMap_dart_iterator_of_involution<Self,1>     I1(*this, adart1);
      CMap_dart_iterator_of_involution_inv<Self,1> I2(*this, adart2);
      while ( I1.cont() )        
      {
        if ( is_marked(I1,m) )
          basic_link_beta_1(I1, I2);
        else
          basic_link_beta_0(I1, I2);
        ++I1; ++I2;
      }

      for (typename std::vector<Dart_handle>::iterator 
             it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }
      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }

    /** Topological sew by beta0 two given darts plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes 
     * thus the map can be non valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable_0(adart1, adart2).
     */
    void topo_sew_0(Dart_handle adart1, Dart_handle adart2)
    { topo_sew_1(adart2, adart1); }

    /** Topological sew by betai two given darts plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes 
     * thus the map can be non valid.   
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre 2<=i<=dimension.
     * @pre is_sewable_for_involution<i>(adart1, adart2).
     */
    template<unsigned int i>
    void topo_sew_for_involution(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(2<=i && i<=Self::dimension);
      CGAL_assertion( (is_sewable_for_involution<i>(adart1,adart2)) );
      
      CMap_dart_iterator_of_involution<Self,i>     I1(*this, adart1);
      CMap_dart_iterator_of_involution_inv<Self,i> I2(*this, adart2);
      while ( I1.cont() )        
      {
        basic_link_beta_for_involution(I1, I2, i);
        ++I1; ++I2;
      }
    }

    /** Topological sew by betai two given darts plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes 
     * thus the map can be non valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable<i>(adart1, adart2).
     */
    template<unsigned int i>
    void topo_sew(Dart_handle adart1, Dart_handle adart2)
    {
      if ( i==0 ) topo_sew_0(adart1, adart2);
      else if ( i==1 ) topo_sew_1(adart1, adart2);
      else topo_sew_for_involution<i>(adart1, adart2);
    }

    /** Sew by beta0 the two given darts plus all the required darts
     * to satisfy the combinatorial map validity, and updates enabled 
     * attributes when necessary so that the final map is valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable_0(adart1, adart2).
     * @post is_valid()
     */
    void sew_0(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( (is_sewable_0(adart1,adart2)) );

      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart1,m);
           it.cont(); ++it)
      {
        mark(it,m);
        dartv.push_back(it);
      }

      CMap_dart_iterator_of_involution<Self,1>     I1(*this, adart1);
      CMap_dart_iterator_of_involution_inv<Self,1> I2(*this, adart2);
      while ( I1.cont() )
      {
        Dart_handle od1=I1->other_extremity();
        Dart_handle od2=I2->other_extremity();
        if (od1!=NULL && od2!=NULL)
          group_all_attributes_except(od1, od2, 1);
        ++I1; ++I2;	  
      }

      I1.rewind(); I2.rewind();      
      while ( I1.cont() )
      {
        if ( is_marked(I1,m) )
          basic_link_beta_0(I1, I2);
        else
          basic_link_beta_1(I1, I2);
        ++I1; ++I2;
      }

      for (typename std::vector<Dart_handle>::iterator 
           it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }
      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }
    
    /** Sew by beta1 the two given darts plus all the required darts
     * to satisfy the combinatorial map validity, and updates enabled 
     * attributes when necessary so that the final map is valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable_1(adart1, adart2).
     * @post is_valid()
     */
    void sew_1(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion( (is_sewable_1(adart1,adart2)) );
      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart1,m); 
           it.cont(); ++it)
      {
        mark(it,m);
        dartv.push_back(it);
      }

      CMap_dart_iterator_of_involution<Self,1>     I1(*this, adart1);
      CMap_dart_iterator_of_involution_inv<Self,1> I2(*this, adart2);
      while ( I1.cont() )
      {
        group_all_attributes_except(I1,I2,1);
        ++I1; ++I2;
      }

      I1.rewind(); I2.rewind();      
      while ( I1.cont() )
      {
        if ( is_marked(I1,m) )
          basic_link_beta_1(I1, I2);
        else
          basic_link_beta_0(I1, I2);
        ++I1; ++I2;
      }

      for (typename std::vector<Dart_handle>::iterator 
             it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }
      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }
    
    /** Sew by betai the two given darts plus all the required darts
     * to satisfy the combinatorial map validity, and updates enabled 
     * attributes when necessary so that the final map is valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable<i>(adart1, adart2).
     * @pre 2<=i<=dimension.
     * @post is_valid()
     */
    template<unsigned int i>
    void sew_for_involution(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_assertion(2<=i && i<=dimension);
      CGAL_assertion( (is_sewable_for_involution<i>(adart1,adart2)) );
      
      CMap_dart_iterator_of_involution<Self,i>     I1(*this, adart1);
      CMap_dart_iterator_of_involution_inv<Self,i> I2(*this, adart2);
      while ( I1.cont() )        
      {
        group_all_attributes_except(I1,I2,i);
        ++I1; ++I2;
      }
      
      I1.rewind(); I2.rewind();      
      while ( I1.cont() )
      {
        basic_link_beta_for_involution(I1, I2, i);
        ++I1; ++I2;
      }
    }

    /** Sew by betai the two given darts plus all the required darts
     * to satisfy the combinatorial map validity, and updates enabled 
     * attributes when necessary so that the final map is valid.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @pre is_sewable<i>(adart1, adart2).
     * @post is_valid()
     */
    template<unsigned int i>
    void sew(Dart_handle adart1, Dart_handle adart2)
    {
      if ( i==0 ) sew_0(adart1, adart2);
      else if ( i==1 ) sew_1(adart1, adart2);
      else sew_for_involution<i>(adart1, adart2);
    }    

    /** Sew by betai the two given darts plus all the required darts
     * to satisfy the combinatorial map validity. Enabled attributes 
     * are updated only if update_attributes==true.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @param update_attributes a boolean to update the enabled attributes
     * @pre is_sewable<i>(adart1, adart2).
     */
    template<unsigned int i>
    void sew(Dart_handle adart1, Dart_handle adart2, bool update_attributes)
    {
      if ( update_attributes ) sew<i>(adart1, adart2);
      else topo_sew<i>(adart1, adart2);
    }

    /** Topological unsew by beta1 the given dart plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes
     * thus the map can be non valid
     * @param adart first dart.
     * @pre !adart->is_free(1).
     */
    void topo_unsew_1(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL && !adart->is_free(1) );

      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart,m);
           it.cont(); ++it)
      {
        mark(*it,m);
        dartv.push_back(*it);
      }

      {
        CMap_dart_iterator_of_involution<Self,1> it(*this, adart);
        while ( it.cont() )
        {
          if ( is_marked(*it,m) ) basic_unlink_beta_1(*it);
          else basic_unlink_beta_0(*it);
          ++it;
        }
      }

      for (typename std::vector<Dart_handle>::iterator
             it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }
      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }

    /** Topological unsew by beta0 the given dart plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes
     * thus the map can be non valid
     * @param adart first dart.
     * @pre !adart->is_free(0).
     */
    void topo_unsew_0(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL && !adart->is_free(0) );
      topo_unsew_1(adart->beta(0));
    }

    /** Topological unsew by betai the given dart plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes
     * thus the map can be non valid
     * @param adart first dart.
     * @pre !adart->is_free(i).
     * @pre 2<=i<=dimension.
     */
    template<unsigned int i>
    void topo_unsew_for_involution(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL && !adart->is_free(i) );
      CGAL_assertion(2<=i && i<=Self::dimension);

      CMap_dart_iterator_of_involution<Self,i> it(*this, adart);
      while ( it.cont() )
      {
        unlink_beta(*it, i);
        ++it;
      }
    }

    /** Topological unsew by betai the given dart plus all the required darts
     * to satisfy the combinatorial map validity: but do not update attributes
     * thus the map can be non valid
     * @param adart first dart.
     * @pre !adart->is_free(i).
     */
    template<unsigned int i>
    void topo_unsew(Dart_handle adart)
    {
      if ( i==0 ) topo_unsew_0(adart);
      else if ( i==1 ) topo_unsew_1(adart);
      else topo_unsew_for_involution<i>(adart);
    }

    /** Unsew by beta0 the given dart plus all the required darts
     * to satisfy the combinatorial map validity, and update enabled
     * attributes when necessary so that the final map is valid.
     * @param adart first dart.
     * @pre !adart->is_free(0).
     * @post is_valid()
     */
    void unsew_0(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL && !adart->is_free(0) );
      Dart_handle d2 = NULL;

      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart,m);
           it.cont(); ++it)
      {
        mark(it,m);
        dartv.push_back(it);        
      }

      {
        CMap_dart_iterator_of_involution<Self,1> it(*this, adart);
        while ( it.cont() )
        {
          if ( is_marked(it,m) )
          {
            d2 = it->beta(0);
            unlink_beta_0(it);
          }
          else
          {
            d2 = it->beta(1);
            unlink_beta_1(it);
          }

          // TODO do the degroup after the loop (cf unsew_for_involution)
          Dart_handle od1=it->other_extremity();
          Dart_handle od2=d2->other_extremity();
          if ( od1!=NULL && od2!=NULL )
            degroup_all_attributes_except(od1,od2,1);

          ++it;
        }
      }

      for (typename std::vector<Dart_handle>::iterator
             it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }     

      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }

    /** Unsew by beta1 the given dart plus all the required darts
     * to satisfy the combinatorial map validity, and update enabled
     * attributes when necessary so that the final map is valid.
     * @param adart first dart.
     * @pre !adart->is_free(1).
     * @post is_valid()
     */
    void unsew_1(Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL && !adart->is_free(1) );
      Dart_handle d2 = NULL;

      int m = get_new_mark();
      std::vector<Dart_handle> dartv;
      for (CMap_dart_iterator_basic_of_cell<Self,0> it(*this,adart,m);
           it.cont(); ++it)
      {
        mark(it,m);
        dartv.push_back(it);
      }

      {
        CMap_dart_iterator_of_involution<Self,1> it(*this, adart);
        while ( it.cont() )
        {
          if ( is_marked(it,m) )
          { d2 = it->beta(1); unlink_beta_1(it); }
          else
          { d2 = it->beta(0); unlink_beta_0(it); }
          // TODO do the degroup after the loop (cf unsew_for_involution)
          degroup_all_attributes_except(it,d2,1);
          ++it;
        }
      }

      for (typename std::vector<Dart_handle>::iterator
             it=dartv.begin(); it!=dartv.end(); ++it)
      { unmark(*it,m); }
      CGAL_assertion( is_whole_map_unmarked(m) );
      free_mark(m);
    }

    /** Unsew by betai the given dart plus all the required darts
     * to satisfy the combinatorial map validity, and update enabled
     * attributes when necessary so that the final map is valid.
     * @param adart first dart.
     * @pre !adart->is_free(i).
     * @post is_valid()
     * @pre 2<=i<=dimension
     */
    template<unsigned int i>
    void unsew_for_involution(Dart_handle adart)
    {
      CGAL_assertion(2<=i && i<=Self::dimension);
      CGAL_assertion( adart!=NULL && !adart->is_free(i) );

      std::stack<internal::Couple_dart_and_dim<Dart_handle> > todegroup;

      CMap_dart_iterator_of_involution<Self,i> it(*this, adart);
      while ( it.cont() )
      {
        todegroup.push(internal::Couple_dart_and_dim<Dart_handle>
                       (it,it->beta(i),i));
        unlink_beta_for_involution(it,i);
        ++it;
      }

      while (!todegroup.empty() )
      {
        internal::Couple_dart_and_dim<Dart_handle> c=todegroup.top();
        todegroup.pop();
        degroup_all_attributes_except(c.d1,c.d2,c.dim);
      }
    }

    /** Unsew by betai the given dart plus all the required darts
     * to satisfy the combinatorial map validity, and update enabled 
     * attributes when necessary so that the final map is valid.
     * @param adart first dart.
     * @pre !adart->is_free(i).
     * @post is_valid()
     */
    template<unsigned int i>
    void unsew(Dart_handle adart)
    {
      if ( i==0 ) unsew_0(adart);
      else if ( i==1 ) unsew_1(adart);
      else unsew_for_involution<i>(adart);
    }

    /** Unsew by betai the given dart plus all the required darts
     * to satisfy the combinatorial map validity. Enabled attributes 
     * are updated only if update_attributes==true.
     * @param adart first dart.
     * @param update_attributes a boolean to update the enabled attributes
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
    count_marked_cells(int amark, const std::vector<unsigned int>& acells) const
    {
      std::vector<unsigned int> res(dimension+2);
      std::vector<int> marks(dimension+2);
      
      // Initialization of the result
      for (unsigned int i=0; i<dimension+2; ++i)
      {
        res[i]=0;
        marks[i]=-1;
      }

      // Mark reservation
      for (unsigned int i=0; i<acells.size(); ++i) 
      {
        CGAL_assertion(acells[i]<=dimension+1);
        if ( marks[acells[i]]==-1 )
        {
          marks[acells[i]] = get_new_mark();
        }
      }

      // Counting and marking cells      
      for (typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        if ( is_marked(it, amark) )
        {
          internal::Foreach_static
            <internal::Count_cell_functor<Self>,dimension+1>::
            run(this, it, &marks, &res);
        }
      }

      // Unmarking darts
      std::vector<unsigned int> tounmark;
      for (unsigned int i=0; i<acells.size(); ++i) 
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
        for (typename Dart_range::const_iterator it(darts().begin()),
               itend(darts().end()); it!=itend; ++it)
        {
          for (unsigned int i=0; i<tounmark.size(); ++i) 
            unmark(it, tounmark[i]);
        }
        for (unsigned int i=0; i<tounmark.size(); ++i) 
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
      int m = get_new_mark();
      negate_mark(m); // We mark all the cells.

      res = count_marked_cells(m,acells);

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
      
      for (unsigned int i=0; i<dimension+2; ++i)
        dim[i]=i;

      return count_cells(dim);
    }
    
  protected:
    /** Set simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @param amarks the marks to set.
     */
    void set_marks(Dart_handle adart,
                   const std::bitset<NB_MARKS> & amarks) const
    {
      CGAL_assertion(adart != NULL && adart!=null_dart_handle);
      adart->set_marks(amarks ^ mmask_marks);
    }

    /** Get simultaneously all the marks of a given dart.
     * @param adart the dart.
     * @return allt the marks of adart.
     */
    std::bitset<NB_MARKS> get_marks(Dart_handle adart) const
    {
      CGAL_assertion(adart != NULL && adart!=null_dart_handle);
      return adart->get_marks() ^ mmask_marks;
    }

    /** Get the mask associated to a given mark.
     * @param amark the mark.
     * @return the mask associated to mark amark.
     */
    bool get_mask_mark(int amark) const
    {
      CGAL_assertion(amark>=0 && (size_type)amark<NB_MARKS);
      return mmask_marks[(size_type)amark];
    }

    /* Decrease the cell attribute reference counting of the given dart.
     * The cell is removed if there is no more darts linked with it.
     */
    template<unsigned int i>
    void decrease_attribute_ref_counting(Dart_handle adart)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "decrease_attribute_ref_counting<i> but "
                                "i-attributes are disabled");
      if ( adart->template attribute<i>()!=NULL )
      { 
        adart->template attribute<i>()->dec_nb_refs();
        if ( adart->template attribute<i>()->get_nb_refs()==0 ) 
          erase_attribute<i>(adart->template attribute<i>());
      }
    }

    /** Update the dart of the given i-cell attribute onto a non marked dart.
     * @param ah the attribute to update.
     * @param amark the mark.
     */
    template<unsigned int i>
    void update_dart_of_attribute(Dart_handle ah, int amark)
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "update_dart_of_attribute<i> but "
                                "i-attributes are disabled");
      CGAL_assertion(ah!=NULL && ah!=null_dart_handle);

      if ( ah->template attribute<i>()==NULL || 
           ah->template attribute<i>()->dart()==NULL ||
           !is_marked(ah->template attribute<i>()->dart(),amark) )
        return;

      for (CMap_dart_iterator_of_cell<Self,i> it(*this, ah); it.cont(); ++it)
      {
        if (!is_marked(it,amark))
        {
          ah->template attribute<i>()->set_dart(it);
          return;
        }
      }
      ah->template attribute<i>()->set_dart(NULL);      
    }
  
    /** Update the dart of all the cell-attributes incident to ah onto a non
     *  marked dart. This method is used before to remove a cell (which is
     *  marked), to put all the darts of the enabled cells onto surviving dart.
     * @param ah the dart to update.
     * @param amark the mark.
     */
    void update_dart_of_all_attributes(Dart_handle ah, int amark)
    { 
      Helper::template Foreach_enabled_attributes
        <internal::Update_dart_of_attribute_functor<Self> >::run(this,ah,amark);
    }  

    /** Group the i cell-attributes of two darts.
     * If the two i cell-attribute of \em adart1 and \em adart2 are different,
     * we set the i cell-attribute of each dart belonging to the i-cell orbit
     * of \em adart2 onto the i-cell of \em adart1.
     * We use the functor On_merge and possibly remove the second 
     * cell-attribute if there is no more darts linked to it.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     */
    template<unsigned int i, class Type_attr>
    void group_enabled_attribute( Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "group_enabled_attribute<i> but "
                                "i-attributes are disabled");
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);

      typename Attribute_handle<i>::type a1=adart1->template attribute<i>();
      typename Attribute_handle<i>::type a2=adart2->template attribute<i>();

      // If the two attributes are equal, nothing to do.
      if ( a1 == a2 ) return;

      Dart_handle toSet = NULL;

      // If the attribute associated to adart1 is NULL, set it with
      // the attribute associated to adart2 (necessarily != NULL)
      if (a1 == NULL)
      { toSet  = adart1; a1 = a2; }
      else
      {
        toSet = adart2;
        if (a2 != NULL)
        {
          internal::Apply_cell_functor<Type_attr,
            typename Type_attr::On_merge>::run(*a1,*a2);
        }
      }
      set_attribute<i>(toSet, a1);
    }

    /** Group all the attributes of adart1 and adart2. If both dart have a 
     * i-attribute, the attribute associated to adart1 is kept.
     * @param adart1 the first dart.
     * @param adart1 the second dart.
     */
    template<unsigned int i>
    void group_attribute(Dart_handle adart1, Dart_handle adart2)
    {
      internal::Group_one_attribute_functor<Self,i,
        typename Attribute_type<i>::type>::run(this,adart1,adart2);
    }

    /** Group the i cell-attributes of two darts.
     * If the two i cell-attribute of \em adart1 and \em adart2 are different,
     * we set the i cell-attribute of adart2 onto the attribute of adart1.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     */
    template<unsigned int i, class Type_attr>
    void group_enabled_attribute_of_dart( Dart_handle dh1, Dart_handle dh2)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "group_enabled_attribute_of_dart<i> but "
                                "i-attributes are disabled");
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=null_dart_handle && dh2!=null_dart_handle );

      typename Attribute_handle<i>::type a1=dh1->template attribute<i>();
      typename Attribute_handle<i>::type a2=dh2->template attribute<i>();

      // If the two attributes are equal, nothing to do.
      if ( a1 == a2 ) return;

      if ( a1==NULL ) set_attribute_of_dart<i>(dh1, a2);
      else            set_attribute_of_dart<i>(dh2, a1);
    }

    /** Group all the dart attributes of adart1 and adart2, except the
     *  adim-cell attribute.
     * @param adart1 the first dart.
     * @param adart1 the second dart.
     * @param adim   the dimension to not group (-1 to group all dimensions).
     * note that 0-attr are always grouped if adart1-> other_extremity()!=NULL.
     */
    void group_all_dart_attributes_except(Dart_handle adart1,
                                          Dart_handle adart2, int adim)
    {
      CGAL_assertion( adim==-1 || (1<=adim && (unsigned int)adim<=dimension) );
      Helper::template Foreach_enabled_attributes
        <internal::Group_attribute_functor_of_dart<Self> >::
        run(this,adart1,adart2,adim);
    }

    /** Group all the cells attributes of adart1 and adart2, except the
     *  adim-cell attribute.
     * @param adart1 the first dart.
     * @param adart1 the second dart.
     * @param adim   the dimension to not group (-1 to group all dimensions).
     * note that 0-attr are always grouped if adart1-> other_extremity()!=NULL.
     */
    void group_all_attributes_except(Dart_handle adart1, Dart_handle adart2,
                                     int adim)
    { 
      CGAL_assertion( adim==-1 || (1<=adim && (unsigned int)adim<=dimension) );
      Helper::template Foreach_enabled_attributes
        <internal::Group_attribute_functor<Self> >::
        run(this,adart1,adart2,adim);
    }
  
    /** Degroup the i-cell attribute of the two given darts, if required.
     * If the two darts are incident to the same attribute and do not belong to
     * the same i-cell, we create a new attibute and link each dart
     * belonging to the i-cell of adart2 onto this new attribute.
     * The new attribute is initialized by copying the old one then by using
     * the functor On_split from the original attribute.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff the attribute is split.
     */
    template<unsigned int i, class Type_attr>
    bool degroup_enabled_attribute(Dart_handle adart1, Dart_handle adart2)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "degroup_enabled_attribute<i> but "
                                "i-attributes are disabled");
      CGAL_assertion(adart1 != NULL && adart2 != NULL);
      CGAL_assertion(adart1 != null_dart_handle && adart2 != null_dart_handle);

      typename Attribute_handle<i>::type a1=adart1->template attribute<i>();
      typename Attribute_handle<i>::type a2=NULL;

      // If the two attributes are not equal, nothing to do.
      if ( a1 != adart2->template attribute<i>() || a1 == NULL) return false;

      // If the two darts belong to the same cell, nothing to do.
      if ( belong_to_same_cell<Self,i,dimension>(*this, adart1, adart2) )
        return false;

      // Here we create a new attribute
      a2 = create_attribute<i>(*a1);

      // We call the on_split functor
      internal::Apply_cell_functor<Type_attr,
        typename Type_attr::On_merge>::run(*a1,*a2);

      // We set the dart of the cell a1 onto adart1.
      a1->set_dart(adart1);

      // And we set all the dart of the cell of adart2 to v1.
      set_attribute<i>(adart2, a2);

      return true;
    }

    template<unsigned int i>
    bool degroup_attribute(Dart_handle adart1, Dart_handle adart2)
    {
      return internal::Degroup_one_attribute_functor<Self,i,
        typename Attribute_type<i>::type>::run(this,adart1,adart2);
    }

    /** Degroup all the cells attributes of adart1 and adart2, except the
     *  adim-cell attribute.
     * @param adart1 the first dart.
     * @param adart1 the second dart.
     * @param adim the dimension to not degroup (-1 to degroup all).
     */
    void degroup_all_attributes_except(Dart_handle adart1, Dart_handle adart2,
                                       int adim)
    { 
      CGAL_assertion( adim==-1 || (1<=adim && (unsigned int)adim<=dimension) );
      Helper::template Foreach_enabled_attributes
        <internal::Degroup_attribute_functor<Self> >::
        run(this,adart1,adart2,adim);
    }
  
    /** Degroup all the cells attributes of adart1 and adart2.
     * @param adart1 the first dart.
     * @param adart1 the second dart.
     */
    void degroup_all_attributes(Dart_handle adart1, Dart_handle adart2)
    { degroup_all_attributes_except(adart1, adart2, -1); }
  
    /** Degroup the i-cell attribute of the two given darts, if required.
     * We create a new attibute and link dart2 onto this new attribute.
     * The new attribute is initialized by copying the old one then by using
     * the functor On_split from the original attribute.
     * @param adart1 the first dart.
     * @param adart2 the second dart.
     * @return true iff the attribute is split.
     */
    template<unsigned int i, class Type_attr, typename Range>
    bool degroup_enabled_attribute_of_dart( Dart_handle dh1, Dart_handle dh2)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "group_enabled_attribute_of_dart<i> but "
                                "i-attributes are disabled");
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=null_dart_handle && dh2!=null_dart_handle );

      typename Attribute_handle<i>::type a1=dh1->template attribute<i>();
      typename Attribute_handle<i>::type a2=dh2->template attribute<i>();

      // If the two attributes are equal, nothing to do.
      if ( a1==NULL || a1 != a2 ) return false;
      
      a2 = create_attribute<i>(*a1);

      // We call the on_split functor
      //      internal::Apply_cell_functor<Type_attr,
      //        typename Type_attr::On_split>::run(*a1,*a2);

      // We set the attribute of dh2 to a2.
      for (typename Range::iterator it=Range(*this,dh2).begin(), 
             itend=Range(*this,dh2).end(); it!=itend; ++it)
      {
        set_attribute_of_dart<i>(it, a2);
      }
      return true;
    }

    template<unsigned int i,typename Range>
    bool degroup_attribute_of_dart(Dart_handle adart1, Dart_handle adart2)
    {
      return internal::Degroup_one_attribute_of_dart_functor<Self,i,
        typename Attribute_type<i>::type, Range>::run(this,adart1,adart2);
    }

    /** Test the validity of a i-cell-attribute.
     * ie all the darts belonging to a i-cell are linked to the same attribute.
     * @param adart a dart.
     * @param amark a mark used to mark darts of the i-cell.
     * @return true iff all the darts of the i-cell link to the same attribute.
     */
    template<unsigned int i>
    bool is_valid_attribute(Dart_const_handle adart, 
                            unsigned int amark) const
    {
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "is_valid_attribute<i> but i-attributes"
                                " are disabled");
      if ( is_marked(adart, amark) ) return true;
      bool valid = true;
      bool found_dart = false;

      typename Attribute_const_handle<i>::type
        a=adart->template attribute<i>();
      
      unsigned int nb = 0;
      for (CMap_dart_const_iterator_basic_of_cell<Self,i> 
             it(*this,adart,amark); it.cont(); ++it)
      {
        if ( it->template attribute<i>() != a )
          valid = false; 

        if ( a!=NULL && it==a->dart() ) found_dart = true;

        mark(it, amark);
        ++nb;
      }

      if ( a!=NULL && a->get_nb_refs()!=nb )
        valid = false;

      if ( a!=NULL && a->dart()!=NULL && !found_dart )
        valid = false;

      return valid;
    }

    /** Erase marked darts from the map.
     * Marked darts are unlinked before to be removed, thus surviving darts
     * are correctly linked, but the map is not necessarily valid depending
     * on the configuration of removed darts. User must check carefully marked
     * darts before calling this method.
     * @param amark the mark of darts to erase.
     * @return the number of removed darts.
     */
    unsigned int erase_marked_darts(int amark)
    {
      unsigned int res = 0, i = 0;
      Dart_handle d;
      for (typename Dart_range::iterator it(darts().begin()),
             itend(darts().end()); it!=itend; )
      {
        d = it++;
        if (is_marked(d, amark))
        {
          for (i = 0; i <= dimension; ++i)
          { if (!d->is_free(i)) unlink_beta(d, i); }
          erase_dart(d); ++res;
        }
      }
      return res;
    }

  public:
    /** Set the i th attribute of the given dart.
     * @param adart a dart.
     * @param ah the attribute to set.
     */
    template<unsigned int i>
    void set_attribute_of_dart(Dart_handle adart, 
                               typename Attribute_handle<i>::type ah)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                                "set_attribute_of_dart<i> but "
                                "i-attributes are disabled");
      CGAL_assertion( adart!=NULL && adart!=null_dart_handle && ah!=NULL );
      if ( adart->template attribute<i>()==ah ) return;

      decrease_attribute_ref_counting<i>(adart);

      adart->template set_attribute<i>(ah);
      ah->set_dart(adart);
    }

    /** Set the i th attribute of all the darts of a given i-cell.
     * @param adart a dart of the i-cell.
     * @param ah the vertex to set.
     */
    template<unsigned int i>
    void set_attribute(Dart_handle adart, 
                       typename Attribute_handle<i>::type ah)
    {
      CGAL_static_assertion(i<=dimension);
      CGAL_static_assertion_msg(Helper::template Dimension_index<i>::value>=0,
                  "set_attribute<i> but i-attributes are disabled");
      CGAL_assertion( adart!=NULL && adart!=null_dart_handle && ah!=NULL );
      for (CMap_dart_iterator_of_cell<Self,i> it(*this, adart); 
           it.cont(); ++it)
      {
        if ( it->template attribute<i>()!=ah )
        {
          decrease_attribute_ref_counting<i>(it);
          it->template set_attribute<i>(ah);
        }
      }
      ah->set_dart(adart);
    }

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    //**************************************************************************
    // Dart_of_orbit_basic_range
    template<unsigned int ... Beta>
    struct Dart_of_orbit_basic_range : public CMap_range
    <Self, CMap_dart_iterator_basic_of_orbit<Self,Beta...>,
     CMap_dart_const_iterator_basic_of_orbit<Self,Beta...> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_basic_of_orbit<Self,Beta...>,
       CMap_dart_const_iterator_basic_of_orbit<Self,Beta...> > Base;
    
      Dart_of_orbit_basic_range(Self &amap, Dart_handle adart, int amark=-1):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_basic_const_range
    template<unsigned int ... Beta>
    struct Dart_of_orbit_basic_const_range : public CMap_const_range
    <Self, CMap_dart_const_iterator_basic_of_orbit<Self,Beta...> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_basic_of_orbit<Self,Beta...> > Base;
      
      Dart_of_orbit_basic_const_range(const Self &amap, Dart_const_handle adart,
                                      int amark=-1):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_range
    template<unsigned int ... Beta>
    struct Dart_of_orbit_range : public CMap_range
    <Self, CMap_dart_iterator_of_orbit<Self,Beta...>,
     CMap_dart_const_iterator_of_orbit<Self,Beta...> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_of_orbit<Self,Beta...>,
       CMap_dart_const_iterator_of_orbit<Self,Beta...> > Base;
    
      Dart_of_orbit_range(Self &amap, Dart_handle adart) : Base(amap,adart)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_const_range
    template<unsigned int ... Beta>
    struct Dart_of_orbit_const_range : public CMap_const_range
    <Self, CMap_dart_const_iterator_of_orbit<Self,Beta...> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_of_orbit<Self,Beta...> > Base;
      
      Dart_of_orbit_const_range(const Self &amap, Dart_const_handle adart):
        Base(amap,adart)
      {}
    };
    //**************************************************************************
    /// @return a range on all the darts of the given orbit
    template<unsigned int ... Beta>
    Dart_of_orbit_range<Beta...> darts_of_orbit(Dart_handle adart)
    { return Dart_of_orbit_range<Beta...>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Beta>
    Dart_of_orbit_const_range<Beta...> 
    darts_of_orbit(Dart_const_handle adart) const
    { return Dart_of_orbit_const_range<Beta...>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Beta>
    Dart_of_orbit_basic_range<Beta...> darts_of_orbit_basic(Dart_handle adart,
                                                            int amark=-1)
    { return Dart_of_orbit_basic_range<Beta...>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int ... Beta>
    Dart_of_orbit_basic_const_range<Beta...> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<Beta...>(*this,adart,amark); }
    //**************************************************************************
#else
    //**************************************************************************
    // Dart_of_orbit_basic_range    
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>    
    struct Dart_of_orbit_basic_range: public CMap_range
    <Self, CMap_dart_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9>,
     CMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9>,
       CMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,
                                               B6,B7,B8,B9> > Base;
      
      Dart_of_orbit_basic_range(Self &amap, Dart_handle adart, int amark=-1): 
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_basic_const_range    
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>    
    struct Dart_of_orbit_basic_const_range: public CMap_const_range
    <Self,
     CMap_dart_const_iterator_basic_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_basic_of_orbit
       <Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> > Base;
      
      Dart_of_orbit_basic_const_range(const Self &amap, Dart_const_handle adart,
                                      int amark=-1): 
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_orbit_range    
    template<int B1=-1,int B2=-1,int B3=-1,int B4=-1,int B5=-1,
             int B6=-1,int B7=-1,int B8=-1,int B9=-1>    
    struct Dart_of_orbit_range: public CMap_range
    <Self, CMap_dart_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9>,
     CMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9>,
       CMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
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
    <Self, CMap_dart_const_iterator_of_orbit<Self,B1,B2,B3,B4,B5,B6,B7,B8,B9> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_of_orbit
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
                                                     int amark=-1)
    { return Dart_of_orbit_basic_range<>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    Dart_of_orbit_basic_const_range<> darts_of_orbit_basic
    (Dart_const_handle adart,int amark=-1) const
    { return Dart_of_orbit_basic_const_range<>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_basic_range<B1> darts_of_orbit_basic(Dart_handle adart,
                                                       int amark=-1)
    { return Dart_of_orbit_basic_range<B1>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1>
    Dart_of_orbit_basic_const_range<B1> darts_of_orbit_basic
    (Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_basic_range<B1,B2> darts_of_orbit_basic(Dart_handle adart,
                                                          int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2>
    Dart_of_orbit_basic_const_range<B1,B2> darts_of_orbit_basic
    (Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_basic_range<B1,B2,B3> darts_of_orbit_basic(Dart_handle adart,
                                                             int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3>
    Dart_of_orbit_basic_const_range<B1,B2,B3> darts_of_orbit_basic
    (Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_basic_range<B1,B2,B3,B4> darts_of_orbit_basic
    (Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5> darts_of_orbit_basic
    (Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6> darts_of_orbit_basic
    (Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7> darts_of_orbit_basic
    (Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8> darts_of_orbit
    (Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8,B9> 
    darts_of_orbit_basic(Dart_handle adart, int amark=-1)
    { return Dart_of_orbit_basic_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template <unsigned int B1,unsigned int B2,unsigned int B3,unsigned int B4,
              unsigned int B5,unsigned int B6,unsigned int B7,unsigned int B8,
              unsigned int B9>
    Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9> 
    darts_of_orbit_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_orbit_basic_const_range<B1,B2,B3,B4,B5,B6,B7,B8,B9>
        (*this,adart,amark); }
    //**************************************************************************
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    //**************************************************************************
    // Dart_of_cell_basic_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_basic_range: public CMap_range
    <Self, CMap_dart_iterator_basic_of_cell<Self,i,dim>,
     CMap_dart_const_iterator_basic_of_cell<Self,i,dim> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_basic_of_cell<Self,i,dim>,
       CMap_dart_const_iterator_basic_of_cell<Self,i,dim> > Base;
    
      Dart_of_cell_basic_range(Self &amap, Dart_handle adart, int amark=-1) : 
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_basic_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_basic_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_basic_of_cell<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_basic_of_cell<Self,i,dim> > Base;
      
      Dart_of_cell_basic_const_range(const Self &amap, Dart_const_handle adart,
                                     int amark=-1) : 
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_range: public CMap_range
    <Self,CMap_dart_iterator_of_cell<Self,i,dim>,
     CMap_dart_const_iterator_of_cell<Self,i,dim> >
    {
      typedef CMap_range
      <Self,CMap_dart_iterator_of_cell<Self,i,dim>,
       CMap_dart_const_iterator_of_cell<Self,i,dim> > Base;
    
      Dart_of_cell_range(Self &amap, Dart_handle adart) : 
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_cell_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_cell_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_of_cell<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_of_cell<Self,i,dim> > Base;
      
      Dart_of_cell_const_range(const Self &amap, Dart_const_handle adart) : 
        Base(amap, adart)
      {}
    };
    //--------------------------------------------------------------------------
    /// @return a range on all the darts of the given i-cell
    template<unsigned int i, int dim>
    Dart_of_cell_basic_range<i,dim> darts_of_cell_basic(Dart_handle adart,
                                                        int amark=-1)
    { return Dart_of_cell_basic_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i, int dim>
    Dart_of_cell_basic_const_range<i,dim> darts_of_cell_basic
    (Dart_const_handle adart, int amark=-1) const
    { return Dart_of_cell_basic_const_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_basic_range<i,Self::dimension>
    darts_of_cell_basic(Dart_handle adart, int amark=-1)
    { return darts_of_cell_basic<i,Self::dimension>(adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_cell_basic_const_range<i,Self::dimension> 
    darts_of_cell_basic(Dart_const_handle adart, int amark=-1) const
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
    struct Dart_of_involution_basic_range: public CMap_range
    <Self, CMap_dart_iterator_basic_of_involution<Self,i,dim>,
     CMap_dart_const_iterator_basic_of_involution<Self,i,dim> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_basic_of_involution<Self,i,dim>,
       CMap_dart_const_iterator_basic_of_involution<Self,i,dim> > Base;
    
      Dart_of_involution_basic_range(Self &amap, Dart_handle adart,
                                     int amark=-1):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_basic_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_basic_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_basic_of_involution<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_basic_of_involution<Self,i,dim> > Base;
      
      Dart_of_involution_basic_const_range(const Self &amap,
                                           Dart_const_handle adart,
                                           int amark=-1) : 
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    template<unsigned int i,int dim>
    Dart_of_involution_basic_range<i,dim>
    darts_of_involution_basic(Dart_handle adart, int amark=-1)
    { return Dart_of_involution_basic_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i,int dim>
    Dart_of_involution_basic_const_range<i,dim>
    darts_of_involution_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_involution_basic_const_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_basic_range<i,Self::dimension>
    darts_of_involution_basic(Dart_handle adart, int amark=-1)
    { return Dart_of_involution_basic_range<i,Self::dimension>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_basic_const_range<i,Self::dimension>
    darts_of_involution_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_involution_basic_const_range<i,Self::dimension>
        (*this,adart,amark); }
    //**************************************************************************
    // Dart_of_involution_inv_basic_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_inv_basic_range: public CMap_range
    <Self, CMap_dart_iterator_basic_of_involution_inv<Self,i,dim>,
     CMap_dart_const_iterator_basic_of_involution_inv<Self,i,dim> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_basic_of_involution_inv<Self,i,dim>,
       CMap_dart_const_iterator_basic_of_involution_inv<Self,i,dim> > Base;
    
      Dart_of_involution_inv_basic_range(Self &amap, Dart_handle adart,
                                         int amark=-1):
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_inv_basic_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_inv_basic_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_basic_of_involution_inv<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_basic_of_involution_inv<Self,i,dim> >
      Base;
      
      Dart_of_involution_inv_basic_const_range(const Self &amap,
                                               Dart_const_handle adart,
                                               int amark=-1) : 
        Base(amap, adart, amark)
      {}
    };
    //**************************************************************************
    template<unsigned int i,int dim>
    Dart_of_involution_inv_basic_range<i,dim>
    darts_of_involution_inv_basic(Dart_handle adart, int amark=-1)
    { return Dart_of_involution_inv_basic_range<i,dim>(*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i,int dim>
    Dart_of_involution_inv_basic_const_range<i,dim>
    darts_of_involution_inv_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_involution_inv_basic_const_range<i,dim>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_inv_basic_range<i,Self::dimension>
    darts_of_involution_inv_basic(Dart_handle adart, int amark=-1)
    { return Dart_of_involution_inv_basic_range<i,Self::dimension>
        (*this,adart,amark); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_inv_basic_const_range<i,Self::dimension>
    darts_of_involution_inv_basic(Dart_const_handle adart, int amark=-1) const
    { return Dart_of_involution_inv_basic_const_range<i,Self::dimension>
        (*this,adart,amark); }    
    //**************************************************************************
    // Dart_of_involution_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_range: public CMap_range
    <Self, CMap_dart_iterator_of_involution<Self,i,dim>,
     CMap_dart_const_iterator_of_involution<Self,i,dim> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_of_involution<Self,i,dim>,
       CMap_dart_const_iterator_of_involution<Self,i,dim> > Base;
    
      Dart_of_involution_range(Self &amap, Dart_handle adart) :
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_of_involution<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_of_involution<Self,i,dim> > Base;
      
      Dart_of_involution_const_range(const Self &amap, Dart_const_handle adart):
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
    // Dart_of_involution_inv_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_inv_range: public CMap_range
    <Self, CMap_dart_iterator_of_involution_inv<Self,i,dim>,
     CMap_dart_const_iterator_of_involution_inv<Self,i,dim> >
    {
      typedef CMap_range
      <Self, CMap_dart_iterator_of_involution_inv<Self,i,dim>,
       CMap_dart_const_iterator_of_involution_inv<Self,i,dim> > Base;
    
      Dart_of_involution_inv_range(Self &amap, Dart_handle adart) :
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // Dart_of_involution_inv_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct Dart_of_involution_inv_const_range: public CMap_const_range
    <Self, CMap_dart_const_iterator_of_involution_inv<Self,i,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_dart_const_iterator_of_involution_inv<Self,i,dim> > Base;
      
      Dart_of_involution_inv_const_range(const Self &amap,
                                         Dart_const_handle adart):
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    template<unsigned int i,int dim>
    Dart_of_involution_inv_range<i,dim>
    darts_of_involution_inv(Dart_handle adart)
    { return Dart_of_involution_inv_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i,int dim>
    Dart_of_involution_inv_const_range<i,dim>
    darts_of_involution_inv(Dart_const_handle adart) const
    { return Dart_of_involution_inv_const_range<i,dim>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_inv_range<i,Self::dimension>
    darts_of_involution_inv(Dart_handle adart)
    { return Dart_of_involution_inv_range<i,Self::dimension>(*this,adart); }
    //--------------------------------------------------------------------------
    template<unsigned int i>
    Dart_of_involution_inv_const_range<i,Self::dimension>
    darts_of_involution_inv(Dart_const_handle adart) const
    { return Dart_of_involution_inv_const_range<i,Self::dimension>
        (*this,adart); }
    //**************************************************************************
    // Dart_basic_range
    struct Dart_basic_range {
      typedef CMap_dart_iterator_basic_of_all<Self> iterator;
      typedef CMap_dart_const_iterator_basic_of_all<Self> const_iterator;
      Dart_basic_range(Self &amap) : mmap(amap)
      {}
      iterator begin() { return iterator(mmap); }
      iterator end()   { return iterator(mmap,NULL); }
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,NULL); }
      size_type size()
      { return mmap.number_of_darts(); }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      Self & mmap;
    };
    //**************************************************************************
    // Dart_basic_const_range
    struct Dart_basic_const_range {
      typedef CMap_dart_const_iterator_basic_of_all<Self> const_iterator;
      Dart_basic_const_range(Self &amap) : mmap(amap)
      {}
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,NULL); }
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
    struct One_dart_per_incident_cell_range: public CMap_range
    <Self, CMap_one_dart_per_incident_cell_iterator<Self,i,j,dim>,
     CMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
    {
      typedef CMap_range
      <Self, CMap_one_dart_per_incident_cell_iterator<Self,i,j,dim>,
       CMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> > Base;
    
      One_dart_per_incident_cell_range(Self &amap, Dart_handle adart): 
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // One_dart_per_incident_cell_const_range
    template<unsigned int i,unsigned int j,int dim=Self::dimension>
    struct One_dart_per_incident_cell_const_range: public CMap_const_range
    <Self, CMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
    {
      typedef CMap_const_range
      <Self, CMap_one_dart_per_incident_cell_const_iterator<Self,i,j,dim> >
      Base;
      
      One_dart_per_incident_cell_const_range(const Self &amap, 
                                             Dart_const_handle adart) : 
        Base(amap, adart)
      {}
    };
    //**************************************************************************
    // One_dart_per_cell_range
    template<unsigned int i,int dim=Self::dimension>
    struct One_dart_per_cell_range {
      typedef CMap_one_dart_per_cell_iterator<Self,i,dim> iterator;
      typedef CMap_one_dart_per_cell_const_iterator<Self,i,dim> const_iterator;
      One_dart_per_cell_range(Self &amap) : mmap(amap), msize(0)
      {}
      iterator begin() { return iterator(mmap); }
      iterator end()   { return iterator(mmap,NULL); }
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,NULL); }
      size_type size()
      {
        if (msize==0)
          for ( const_iterator it=begin(); it!=end(); ++it)
            ++msize;
        return msize;
      }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      Self & mmap;
      size_type msize;
    };
    //**************************************************************************
    // One_dart_per_cell_const_range
    template<unsigned int i,int dim=Self::dimension>
    struct One_dart_per_cell_const_range {
      typedef CMap_one_dart_per_cell_const_iterator<Self,i,dim> const_iterator;
      One_dart_per_cell_const_range(const Self &amap) : mmap(amap), msize(0)
      {}
      const_iterator begin() const { return const_iterator(mmap); }
      const_iterator end() const   { return const_iterator(mmap,NULL); }
      size_type size()
      {
        if (msize==0)
          for ( const_iterator it=begin(); it!=end(); ++it)
            ++msize;
        return msize;
      }
      bool empty() const
      { return mmap.is_empty(); }
    private:
      const Self & mmap;
      size_type msize;
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

    /** Compute the dual of a Combinatorial_map.
     * @param amap the cmap in which we build the dual of this map.
     * @param adart a dart of the initial map, NULL by default.
     * @return adart of the dual map, the dual of adart if adart!=NULL,
     *         any dart otherwise.
     * As soon as we don't modify this map and amap map, we can iterate
     * simultaneously through all the darts of the two maps and we have
     * each time of the iteration two "dual" darts.
     */
    Dart_handle dual(Self& amap, Dart_handle adart=NULL)
    {
      CGAL_assertion( is_without_boundary(dimension) );

      std::map< Dart_handle, Dart_handle > dual;
      Dart_handle d, d2, res = NULL;
  
      // We clear amap. TODO return a new amap ? (but we need to make
      // a copy contructor and =operator...)
      amap.clear();
  
      // We create a copy of all the dart of the map.
      for (typename Dart_range::iterator it=darts().begin(); it!=darts().end();
           ++it)
      {
        dual[it] = amap.create_dart();
        if ( it==adart && res==NULL ) res = dual[it];
      }
  
      // Then we link the darts by using the dual formula :
      // G(B,b1,b2,...,bn-1,bn) =>
      //    dual(G)=(B, b(n-1)obn, b(n-2)obn,...,b1obn, bn)
      // We suppose darts are run in the same order for both maps.
      typename Dart_range::iterator it2=amap.darts().begin();
      for (typename Dart_range::iterator it=darts().begin(); it!=darts().end();
           ++it, ++it2)
      {
        d = it2; // The supposition on the order allows to avoid d=dual[it];
        CGAL_assertion(it2 == dual[it]);

        // First case outside the loop since we need to use link_beta1
        if ( d->is_free(1) &&
             it->beta(dimension)->beta(dimension-1)!=null_dart_handle )
          amap.link_beta<1>(d, 
                            dual[it->beta(dimension)->beta(dimension-1)]);

        // and during the loop we use link_beta(d1,d2,i)
        for (unsigned int i=dimension-2; i>=1; --i)
        {
          if ( d->is_free(dimension-i) &&
               it->beta(dimension)->beta(i)!=null_dart_handle )
            amap.link_beta(d, dual[it->beta(dimension)->beta(i)], dimension-i);
        }
        if ( d->is_free(dimension) )
        {
          CGAL_assertion ( !it->is_free(dimension) );
          amap.link_beta(d, dual[it->beta(dimension)],dimension);
        }
      }
  
      //  CGAL_postcondition(amap2.is_valid());

      if ( res==NULL ) res = amap.darts().begin();
      return res;
    }

  public:
    /// Void dart. A dart d is i-free if beta_i(d)=null_dart_handle.
    static Dart_handle null_dart_handle;
    
  protected:
    /// Dart container.
    Dart_container mdarts;

    /// Container for the null_dart_handle, static data member.
    static Dart_container mnull_dart_container;
    
    /// Number of times each mark is reserved. 0 if the mark is free.
    mutable size_type mnb_times_reserved_marks[NB_MARKS];

    /// Mask marks to know the value of unmark dart, for each index i.
    mutable std::bitset<NB_MARKS> mmask_marks;

    /// Number of used marks.
    mutable size_type mnb_used_marks;

    /// Index of each mark, in mfree_marks_stack or in mfree_marks_stack.
    mutable size_type mindex_marks[NB_MARKS];

    /// "Stack" of free marks.
    mutable int mfree_marks_stack[NB_MARKS];

    /// "Stack" of used marks.
    mutable int mused_marks_stack[NB_MARKS];

    /// Number of marked darts for each used marks.
    mutable size_type mnb_marked_darts[NB_MARKS];

    /// Tuple of attributes containers
    typename Helper::Attribute_containers mattribute_containers;
  };

  /// Allocation of static data members
  /// mnull_dart_container
  template < unsigned int d_, class Refs, class Items_, class Alloc_ >
  typename Combinatorial_map_base<d_, Refs, Items_, Alloc_>::Dart_container
  Combinatorial_map_base<d_, Refs, Items_, Alloc_>::mnull_dart_container;
  
  /// null_dart_handle
  template < unsigned int d_, class Refs, class Items_, class Alloc_ >
  typename Combinatorial_map_base<d_, Refs, Items_, Alloc_>::Dart_handle
  Combinatorial_map_base<d_, Refs, Items_, Alloc_>::null_dart_handle;
  // =  mnull_dart_container.emplace( std::bitset<NB_MARKS>() );
  // Does not work on windows => segfault
  // Thus we initialize null_dart_handle in the Combinatorial_map constructor
  
  template < unsigned int d_, 
             class Items_=Combinatorial_map_min_items<d_>,
             class Alloc_=CGAL_ALLOCATOR(int) >
  class Combinatorial_map : 
    public Combinatorial_map_base<d_, 
                                  Combinatorial_map<d_,Items_,Alloc_>, 
                                  Items_, Alloc_ >
  {
  public:
    typedef Combinatorial_map<d_, Items_,Alloc_>  Self;
    typedef Combinatorial_map_base<d_, Self, Items_, Alloc_> Base;
    
    typedef typename Base::Dart_handle Dart_handle;
    typedef typename Base::Dart_const_handle Dart_const_handle;
    typedef typename Base::Alloc Alloc;
  };
  
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_H //
// EOF //
