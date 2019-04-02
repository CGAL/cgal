// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_FACE_GRAPH_WRAPPER_H
#define CGAL_FACE_GRAPH_WRAPPER_H 1

#include <CGAL/Functors_for_face_graph_wrapper.h>
#include <CGAL/Combinatorial_map_iterators_base.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL
{
  
////////////////////////////////////////////////////////////////////////////////
/** Class Face_graph_wrapper: to wrap any model of FaceGraph into a Combinatorial map.
 *  For now, only for const models, i.e. does not support modification operators.
 */
template<typename HEG_>
class Face_graph_wrapper
{
public:
  typedef HEG_                    HEG;
  typedef Face_graph_wrapper<HEG> Self;
  typedef std::size_t             size_type;

  struct Dart_container
  {
    typedef typename boost::graph_traits<HEG>::halfedge_iterator iterator;
    typedef typename boost::graph_traits<HEG>::halfedge_iterator const_iterator; // TODO ?
    // typedef My_halfedge_iterator<HEG> iterator;
    // typedef My_halfedge_iterator<HEG> const_iterator; // TODO ?
  };
  
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_handle;
  typedef typename boost::graph_traits<HEG>::halfedge_descriptor Dart_const_handle;

   typedef Dart_handle Null_handle_type;
  // typedef CGAL::Void* Null_handle_type;
   static const Null_handle_type null_handle; //=Dart_handle();
   static const Null_handle_type null_dart_handle; //=Dart_handle();

  /// Number of marks
  static const size_type NB_MARKS = 32;
  static const size_type INVALID_MARK = NB_MARKS;

  /// The dimension of the combinatorial map.
  static const unsigned int dimension=2;

  typedef typename boost::graph_traits<HEG>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<HEG>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<HEG>::face_descriptor   face_descriptor;
  typedef boost::undirected_tag                                directed_category;
  typedef boost::disallow_parallel_edge_tag                    edge_parallel_category; 

  struct SM_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                       public virtual boost::vertex_list_graph_tag,
                                       public virtual boost::edge_list_graph_tag
  {};
  typedef SM_graph_traversal_category traversal_category;

  
  Face_graph_wrapper(const HEG& f) : m_fg(f),
                                     mdarts(*this),
                                     mnb_used_marks(0),
                                     m_nb_darts(0)

  {
    mmask_marks.reset();

    for (size_type i=0; i<NB_MARKS; ++i)
    {
      mfree_marks_stack[i]       =i;
      mindex_marks[i]            =i;
      mnb_marked_darts[i]        =0;
      mnb_times_reserved_marks[i]=0;
      // m_marks[i]                 =nullptr;
    }

    m_nb_darts=darts().size(); // Store locally the number of darts: the HEG must not be modified
  }

  const HEG& get_fg() const
  { return m_fg; }
  
  template<unsigned int i>
  bool is_free(Dart_const_handle dh) const
  { return Is_free<HEG, i>::value(m_fg, dh); }
  bool is_free(Dart_const_handle dh, unsigned int i) const
  {
    if (i==2) { return Is_free<HEG, 2>::value(m_fg, dh); }
    return false; // Not possible to have 0 or 1 free dart with an HEG.
  }

  Dart_const_handle get_beta(Dart_const_handle ADart, int B1) const
  {
    CGAL_assertion(B1>=0 && B1<=(int)dimension);
    if (B1==1) return Get_beta<HEG, 1>::value(m_fg, ADart);
    if (B1==2) return Get_beta<HEG, 2>::value(m_fg, ADart); 
    return Get_beta<HEG, 0>::value(m_fg, ADart);
  }
  template<int B1>
  Dart_const_handle get_beta(Dart_const_handle ADart) const
  {
    CGAL_assertion(B1>=0 && B1<=(int)dimension);
    return Get_beta<HEG, B1>::value(m_fg, ADart);
  }

  bool is_empty() const
  { return number_of_darts()==0; }
  
  /* ??  bool is_dart_used(Dart_const_handle dh) const
      { return true; ?? } */
  
  int highest_nonfree_dimension(Dart_const_handle dh) const
  { return (Is_free<HEG, 2>(m_fg, dh)?1:2); }

  Dart_const_handle previous(Dart_const_handle ADart) const
  { return get_beta<0>(ADart); }  
  Dart_const_handle next(Dart_const_handle ADart) const
  { return get_beta<1>(ADart); }
  Dart_const_handle opposite(Dart_const_handle dh)
  { return get_beta<2>(dh); }
  Dart_const_handle other_extremity(Dart_const_handle dh)
  { return get_beta<1>(dh); }

  template<unsigned int dim>
  Dart_const_handle opposite(Dart_const_handle ADart) const
  { return this->template get_beta<dim>(ADart); }
  Dart_const_handle other_orientation(Dart_const_handle ADart) const
  { return ADart; }  

  bool is_previous_exist(Dart_const_handle) const
  { return true; }
  bool is_next_exist(Dart_const_handle) const
  { return true; }
  template<unsigned int dim>
  bool is_opposite_exist(Dart_const_handle ADart) const
  { return !this->template is_free<dim>(ADart); }

  template<typename ...Betas>
  Dart_handle beta(Dart_handle ADart, Betas... betas)
  { return CGAL::internal::Beta_functor<Self, Dart_handle, Betas...>::
      run(*this, ADart, betas...); }
  template<typename ...Betas>
  Dart_const_handle beta(Dart_const_handle ADart, Betas... betas) const
  { return CGAL::internal::Beta_functor<const Self, Dart_const_handle, Betas...>::
      run(*this, ADart, betas...); }
  template<int... Betas>
  Dart_handle beta(Dart_handle ADart)
  { return CGAL::internal::Beta_functor_static<Self, Dart_handle, Betas...>::
      run(*this, ADart); }
  template<int... Betas>
  Dart_const_handle beta(Dart_const_handle ADart) const
  { return CGAL::internal::Beta_functor_static<const Self, Dart_const_handle, Betas...>::
      run(*this, ADart); }
  
  size_type number_of_darts() const
  { return m_nb_darts; }
  
  size_type number_of_halfedges() const
  { return number_of_darts(); }

  size_type number_of_used_marks() const
  { return mnb_used_marks; }

  bool is_reserved(size_type amark) const
  {
    CGAL_assertion(amark<NB_MARKS);
    return (mnb_times_reserved_marks[amark]!=0);
  }
  
  size_type number_of_marked_darts(size_type amark) const
  {
    CGAL_assertion( is_reserved(amark) );
    return mnb_marked_darts[amark];
  }
  
  size_type number_of_unmarked_darts(size_type amark) const
  {
    CGAL_assertion( is_reserved(amark) );
    return number_of_darts() - number_of_marked_darts(amark);
  }
  
  bool is_whole_map_unmarked(size_type amark) const
  { return number_of_marked_darts(amark)==0; }
  
  bool is_whole_map_marked(size_type amark) const
  {  return number_of_marked_darts(amark)==number_of_darts(); }

  class Exception_no_more_available_mark {};
  
  size_type get_new_mark() const
  {
    if (mnb_used_marks==NB_MARKS)
    {
      std::cerr << "Not enough Boolean marks: "
        "increase NB_MARKS in item class." << std::endl;
      std::cerr << "  (exception launched)" << std::endl;
      throw Exception_no_more_available_mark();
    }
    
    size_type m=mfree_marks_stack[mnb_used_marks];
    mused_marks_stack[mnb_used_marks]=m;
    
    mindex_marks[m]=mnb_used_marks;
    mnb_times_reserved_marks[m]=1;
    
    ++mnb_used_marks;
    CGAL_assertion(is_whole_map_unmarked(m));

    // m_marks[m]=new std::vector<bool>(CGAL::num_halfedges(m_fg),  mmask_marks[m]); // TODO use property map
    // CGAL::num_halfedges is an upper bound; depending on the removed halfedges

    m_marks_new[m]=get(CGAL::dynamic_halfedge_property_t<bool>(),
                       const_cast<HEG&>(m_fg)); // Strange...
    for (typename Dart_range::const_iterator it(darts().begin()),
         itend(darts().end()); it!=itend; ++it)
    { set_dart_mark(*it, m, mmask_marks[m]); }

    return m;
  }
  
  void share_a_mark(size_type amark) const
  {
    CGAL_assertion( is_reserved(amark) );
    ++mnb_times_reserved_marks[amark];
  }

  size_type get_number_of_times_mark_reserved(size_type amark) const
  {
    CGAL_assertion( amark<NB_MARKS );
    return mnb_times_reserved_marks[amark];
  }

  void negate_mark(size_type amark) const
  {
    CGAL_assertion(is_reserved(amark));
    
    mnb_marked_darts[amark]=number_of_darts()-mnb_marked_darts[amark];
    mmask_marks.flip(amark);
  }

  void mark_null_dart( size_type /*amark*/) const
  {}
  
  bool get_dart_mark(Dart_const_handle ADart, size_type amark) const
  {
    // return (*(m_marks[amark]))[ADart];
    return get(m_marks_new[amark], ADart);
  }
  void set_dart_mark(Dart_const_handle ADart, size_type amark, bool avalue) const
  { // (*(m_marks[amark]))[ADart]=avalue;
    put(m_marks_new[amark], ADart, avalue);
  }

  void flip_dart_mark(Dart_const_handle ADart, size_type amark) const
  { set_dart_mark(ADart, amark, !get_dart_mark(ADart, amark)); }

  bool is_marked(Dart_const_handle adart, size_type amark) const
  {
    CGAL_assertion(is_reserved(amark));
    return get_dart_mark(adart, amark)!=mmask_marks[amark];
  }

  void set_mark_to(Dart_const_handle adart, size_type amark,
                   bool astate) const
  {
    CGAL_assertion(is_reserved(amark));
    
    if (is_marked(adart, amark)!=astate)
    {
      if (astate) ++mnb_marked_darts[amark];
      else --mnb_marked_darts[amark];
      
      flip_dart_mark(adart, amark);
    }
  }
  
  void mark(Dart_const_handle adart, size_type amark) const
  {
    CGAL_assertion(is_reserved(amark));
    
    if (is_marked(adart, amark)) return;
    
    ++mnb_marked_darts[amark];
    flip_dart_mark(adart, amark);
  }
  
  void unmark(Dart_const_handle adart, size_type amark) const
  {
    CGAL_assertion( adart!=this->null_dart_handle );
    CGAL_assertion( is_reserved(amark) );
    
    if (!is_marked(adart, amark)) return;
    
    --mnb_marked_darts[amark];
    flip_dart_mark(adart, amark);
  }
  
  void unmark_all(size_type amark) const
  {
    CGAL_assertion( is_reserved(amark) );
    
    if ( is_whole_map_marked(amark) )
    {
      negate_mark(amark);
    }
    else if ( !is_whole_map_unmarked(amark) )
    {
      for (typename Dart_range::const_iterator it(darts().begin()),
           itend(darts().end()); it!=itend; ++it)
        unmark(*it, amark);
    }
    CGAL_assertion(is_whole_map_unmarked(amark));
  }
  
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
    mfree_marks_stack[ mnb_used_marks ]=amark;
    mindex_marks[amark] = mnb_used_marks;
    
    mnb_times_reserved_marks[amark]=0;

    // delete m_marks[amark]; m_marks[amark]=nullptr; // TODO use property map
    m_marks_new[amark]=get(CGAL::dynamic_halfedge_property_t<bool>(),
                           const_cast<HEG&>(m_fg)); // To erase the property map ??
  }
  
  bool is_without_boundary(unsigned int i) const
  {
    CGAL_assertion(1<=i && i<=dimension);
    if (i==1) return true;
    
    for ( typename Dart_range::const_iterator it(darts().begin()),
            itend(darts().end()); it!=itend; ++it)
    { if (is_free<2>(it)) return false; }
    return true;
  }
  
  bool is_without_boundary() const
  { return is_without_boundary(2); }

  //**************************************************************************
  // Dart_of_cell_range
  template<unsigned int i>
  struct Dart_of_cell_range
  {
    typedef CGAL::FGW_cell_iterator<Self, i> iterator;
    typedef CGAL::FGW_cell_iterator<Self, i> const_iterator;
    Dart_of_cell_range(const Self &amap, Dart_handle adart) : mmap(amap),
                                                              m_initdart(adart),
                                                              msize(0)
    {}
    const_iterator begin() const { return const_iterator(mmap, m_initdart); }
    const_iterator end() const   { return const_iterator(mmap, m_initdart, mmap.null_handle); }
    size_type size() const
    {
      if (msize==0)
      {
        for (const_iterator it=begin(), itend=end(); it!=itend; ++it)
        { ++msize; }
      }
      return msize;
    }
    
    bool empty() const
    { return mmap.is_empty(); }
  private:
    const Self & mmap;
    Dart_handle m_initdart;
    mutable typename Self::size_type msize;
  };
  //**************************************************************************
  // Dart_of_cell_const_range
  /*  template<unsigned int i,int dim=Self::dimension>
  struct Dart_of_cell_const_range // TODO REMOVE ??
  {}; */
  //--------------------------------------------------------------------------
  template<unsigned int i>
  Dart_of_cell_range<i> darts_of_cell(Dart_handle adart)
  { return Dart_of_cell_range<i>(*this,adart); }
  //--------------------------------------------------------------------------
  template<unsigned int i>
  Dart_of_cell_range<i> darts_of_cell(Dart_const_handle adart) const
  { return Dart_of_cell_range<i>(*this,adart); } // Before it was Dart_of_cell_const_range<i>
  //**************************************************************************
  // Dart_range
  struct Dart_range {
    typedef CGAL::FGW_dart_iterator_basic_of_all<Self> iterator;
    typedef CGAL::FGW_dart_iterator_basic_of_all<Self> const_iterator;
    Dart_range(const Self &amap) : mmap(amap), msize(0)
    {}
    iterator begin() { return iterator(mmap); }
    iterator end()   { return iterator(mmap,mmap.null_handle); }
    const_iterator begin() const { return const_iterator(mmap); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
    size_type size() const
    {
      if (msize==0)
      {
        for (const_iterator it=begin(), itend=end(); it!=itend; ++it)
        { ++msize; }
      }
      return msize;
    }
    bool empty() const
    { return mmap.is_empty(); }
  private:
    const Self & mmap;
    mutable typename Self::size_type msize;
  };
  //**************************************************************************
  // Dart_const_range // TODO REMOVE ?
  /*  struct Dart_const_range {
    typedef CGAL::FGW_dart_iterator_basic_of_all<Self, true> const_iterator;
    Dart_const_range(const Self &amap) : mmap(amap), msize(0)
    {}
    const_iterator begin() const { return const_iterator(mmap); }
    const_iterator end() const   { return const_iterator(mmap,mmap.null_handle); }
    size_type size() const
    {
      if (msize==0)
      {
        for (const_iterator it=begin(), itend=end(); it!=itend; ++it)
        { ++msize; }
      }
      return msize;
    }
    
    bool empty() const
    { return mmap.is_empty(); }
  private:
    const Self & mmap;
    mutable typename Self::size_type msize;
    };*/
  //**************************************************************************
  Dart_range darts()
  { return Dart_range(*this); }
  //--------------------------------------------------------------------------
  Dart_range darts() const
  { return Dart_range(*this); } // Before it was Dart_const_range(*this)
  //**************************************************************************

  template <unsigned int i>
  bool belong_to_same_cell(Dart_const_handle adart1,
                           Dart_const_handle adart2) const
  {
    for (typename Dart_of_cell_range<i>::iterator it=darts_of_cell<i>(adart1).begin(),
           itend=darts_of_cell<i>(adart1).begin(); it!=itend; ++it)
    { if (*it==adart2) { return true; } }
    return false;
  }
  
  template <unsigned int i>
  bool is_whole_cell_unmarked(Dart_const_handle adart, size_type amark) const
  {
    for (typename Dart_of_cell_range<i>::iterator it=darts_of_cell<i>(adart).begin(),
           itend=darts_of_cell<i>(adart).begin(); it!=itend; ++it)
    { if (is_marked(*it, amark)) { return false; } }
    return true;
  }

  template <unsigned int i>
  bool is_whole_cell_marked(Dart_const_handle adart, size_type amark) const
  {
    for (typename Dart_of_cell_range<i>::iterator it=darts_of_cell<i>(adart).begin(),
           itend=darts_of_cell<i>(adart).end(); it!=itend; ++it)
    { if (!is_marked(*it, amark)) { return false; } }
    return true;
  }

  template <unsigned int i>
  size_type mark_cell(Dart_const_handle adart, size_type amark) const
  {
    size_type res=0;
    for (typename Dart_of_cell_range<i>::iterator it=darts_of_cell<i>(adart).begin(),
           itend=darts_of_cell<i>(adart).end(); it!=itend; ++it)
    { mark(*it, amark); ++res; }
    return res;
}

  size_type mark_cell(Dart_const_handle adart, unsigned int i, size_type amark) const
  {
    if (i==0) { return mark_cell<0>(adart, amark); }
    else if (i==1) { return mark_cell<1>(adart, amark); }
    return mark_cell<2>(adart, amark);
  } 

  template <unsigned int i>
  size_type unmark_cell(Dart_const_handle adart, size_type amark) const
  {
    size_type res=0;
    for (typename Dart_of_cell_range<i>::iterator it=darts_of_cell<i>(adart).begin(),
           itend=darts_of_cell<i>(adart).end(); it!=itend; ++it)
    { unmark(*it, amark); ++res; }
    return res;
  }

  std::vector<unsigned int>
  count_marked_cells(size_type amark, const std::vector<unsigned int>& acells) const
  {
    std::vector<unsigned int> res(dimension+2);
    std::vector<size_type> marks(dimension+2);

    // Initialization of the result
    for (unsigned int i=0; i<dimension+2; ++i)
    {
      res[i]=0;
      marks[i]=INVALID_MARK;
    }
    
    // Mark reservation
    for (unsigned int i=0; i<acells.size(); ++i)
    {
      CGAL_assertion(acells[i]<=dimension+1);
      if (marks[acells[i]]==INVALID_MARK )
      {
        marks[acells[i]]=get_new_mark();
        assert(is_whole_map_unmarked(marks[acells[i]]));
      }
    }

    // Counting and marking cells
    for (typename Dart_range::const_iterator it(darts().begin()),
           itend(darts().end()); it!=itend; ++it)
    {
      if (is_marked(*it, amark))
      {
        for (unsigned int i=0; i<acells.size(); ++i)
        {
          if (!is_marked(*it, marks[acells[i]]))
          {
            mark_cell(*it, acells[i], marks[acells[i]]);
            ++res[acells[i]];
          }
        }
      }
    }

    // Unmarking darts
    std::vector<size_type> tounmark;
    for (unsigned int i=0; i<acells.size(); ++i)
    {
      if (is_whole_map_marked(marks[acells[i]]) ||
          is_whole_map_unmarked(marks[acells[i]]))
      { free_mark(marks[acells[i]]); }
      else
      { tounmark.push_back(marks[acells[i]]); }
    }

    if (tounmark.size()>0)
    {
      for (typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        for (unsigned int i=0; i<tounmark.size(); ++i)
        { unmark(*it, tounmark[i]); }
      }
      for (unsigned int i=0; i<tounmark.size(); ++i)
      {
        CGAL_assertion(is_whole_map_unmarked(tounmark[i]));
        free_mark(tounmark[i]);
      }
    }

    return res;
  }

  std::vector<unsigned int>
  count_cells(const std::vector<unsigned int>& acells) const
  {
    std::vector<unsigned int> res;
    size_type m=get_new_mark();
    negate_mark(m); // We mark all the cells.
    
    res=count_marked_cells(m, acells);
    
    negate_mark(m); // We unmark the cells
    free_mark(m);
    
    return res;
  }

  std::vector<unsigned int> count_all_cells() const
  {
    std::vector<unsigned int> dim(dimension+1);
    
    for ( unsigned int i=0; i<=dimension; ++i)
    { dim[i]=i; }
    
    return count_cells(dim);
  }

  std::ostream& display_characteristics(std::ostream & os) const
  {
    std::vector<unsigned int> cells(dimension+1);
    for ( unsigned int i=0; i<=dimension; ++i)
    { cells[i]=i; }
    
    std::vector<unsigned int> res=count_cells(cells);
    
    os<<"#Darts="<<number_of_darts();
    for (unsigned int i=0; i<=dimension; ++i)
      os<<", #"<<i<<"-cells="<<res[i];
    //  os<<", #ccs="<<res[dimension+1];

    return os;
  }

protected:
  const HEG& m_fg;
  Dart_range mdarts;
  std::size_t m_nb_darts;
  
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

  /// Array of property maps; one for each reserved mark.
  // mutable std::vector<bool>* m_marks[NB_MARKS];
  typedef typename boost::property_map
              <HEG, CGAL::dynamic_halfedge_property_t<bool> >::type MarkPMap;
  mutable MarkPMap m_marks_new[NB_MARKS];
};

  /// null_handle
  // template <typename HEG>
  // const typename Face_graph_wrapper<HEG>::Null_handle_type
  // Face_graph_wrapper<HEG>::null_handle=nullptr;
  template <typename HEG>
  const typename Face_graph_wrapper<HEG>::Null_handle_type
  Face_graph_wrapper<HEG>::null_handle=typename Face_graph_wrapper<HEG>::Dart_handle();

  /// null_dart_handle
  // template <typename HEG>
  // const typename Face_graph_wrapper<HEG>::Null_handle_type
  // Face_graph_wrapper<HEG>::null_dart_handle=nullptr;
  template <typename HEG>
  const typename Face_graph_wrapper<HEG>::Null_handle_type
  Face_graph_wrapper<HEG>::null_dart_handle=typename Face_graph_wrapper<HEG>::Dart_handle();

} // Namespace CGAL

#endif // CGAL_FACE_GRAPH_WRAPPER_H //
// EOF //
