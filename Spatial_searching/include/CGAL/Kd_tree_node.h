// Copyright (c) 2002,2011  Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H

#include <CGAL/license/Spatial_searching.h>



#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Has_member.h>
#include <CGAL/internal/Search_helpers.h>
#include <boost/cstdint.hpp>

namespace CGAL {

  CGAL_GENERATE_MEMBER_DETECTOR(contains_point_given_as_coordinates);

  template <class SearchTraits, class Splitter, class UseExtendedNode, class EnablePointsCache>
  class Kd_tree;

  template < class TreeTraits, class Splitter, class UseExtendedNode, class EnablePointsCache >
  class Kd_tree_node {

    friend class Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache>;

    typedef Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache> Kdt;

    typedef typename Kdt::Node_handle Node_handle;
    typedef typename Kdt::Node_const_handle Node_const_handle;
    typedef typename Kdt::Internal_node_handle Internal_node_handle;
    typedef typename Kdt::Internal_node_const_handle Internal_node_const_handle;
    typedef typename Kdt::Leaf_node_handle Leaf_node_handle;
    typedef typename Kdt::Leaf_node_const_handle Leaf_node_const_handle;
    typedef typename TreeTraits::Point_d Point_d;

    typedef typename TreeTraits::FT FT;
    typedef typename Kdt::Separator Separator;
    typedef typename Kdt::Point_d_iterator Point_d_iterator;
    typedef typename Kdt::iterator iterator;
    typedef typename Kdt::D D;

    bool leaf;

  public :
    Kd_tree_node(bool leaf) : leaf(leaf) { }

    bool is_leaf() const { return leaf; }

    std::size_t
    num_items() const
    {
      if (is_leaf()){
        Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
        return node->size();
      }
      else {
        Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        return node->lower()->num_items() + node->upper()->num_items();
      }
    }

    std::size_t
    num_nodes() const
    {
      if (is_leaf()) return 1;
      else {
        Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        return node->lower()->num_nodes() + node->upper()->num_nodes();
      }
    }

    int
    depth(const int current_max_depth) const
    {
      if (is_leaf()){
        return current_max_depth;
      }
      else {
        Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        return
          (std::max)( node->lower()->depth(current_max_depth + 1),
                      node->upper()->depth(current_max_depth + 1));
      }
    }

    int
    depth() const
    {
      return depth(1);
    }

    template <class OutputIterator>
    OutputIterator
    tree_items(OutputIterator it) const {
      if (is_leaf()) {
        Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
        if (node->size()>0)
          for (iterator i=node->begin(); i != node->end(); i++)
          {*it=*i; ++it;}
      }
      else {
        Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        it=node->lower()->tree_items(it);
        it=node->upper()->tree_items(it);
      }
      return it;
    }


    boost::optional<Point_d>
    any_tree_item() const {
      boost::optional<Point_d> result = boost::none;
      if (is_leaf()) {
         Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
         if (node->size()>0){
           return boost::make_optional(*(node->begin()));
         }
        }
      else {
         Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
          result = node->lower()->any_tree_item();
          if(! result){
            result = node->upper()->any_tree_item();
          }
      }
      return result;
    }


     void
    indent(int d) const
    {
      for(int i = 0; i < d; i++){
        std::cout << " ";
      }
    }


    void
    print(int d = 0) const
    {
      if (is_leaf()) {
        Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
        indent(d);
        std::cout << "leaf" << std::endl;
        if (node->size()>0)
          for (iterator i=node->begin(); i != node->end(); i++)
          {indent(d);std::cout << *i << std::endl;}
      }
      else {
        Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        indent(d);
        std::cout << "lower tree" << std::endl;
        node->lower()->print(d+1);
        indent(d);
        std::cout << "separator: dim = " << node->cutting_dimension() << "  val = " << node->cutting_value() << std::endl;
        indent(d);
        std::cout << "upper tree" << std::endl;
        node->upper()->print(d+1);
      }
    }


    template <class OutputIterator, class FuzzyQueryItem>
    OutputIterator
    search(OutputIterator it, const FuzzyQueryItem& q,
           Kd_tree_rectangle<FT,D>& b,
           typename Kdt::const_iterator tree_points_begin,
           typename std::vector<FT>::const_iterator cache_begin,
           int dim) const
    {
      if (is_leaf()) {
        Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
        if (node->size() > 0)
        {
          typename internal::Has_points_cache<Kdt, internal::has_Enable_points_cache<Kdt>::type::value>::type dummy;
          it = search_in_leaf(node, q, tree_points_begin, cache_begin, dim, it, dummy);
        }
      }
      else {
         Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
         // after splitting b denotes the lower part of b
         Kd_tree_rectangle<FT,D> b_upper(b);
         node->split_bbox(b, b_upper);

         if (q.outer_range_contains(b))
           it=node->lower()->tree_items(it);
         else
           if (q.inner_range_intersects(b))
             it=node->lower()->search(it,q,b,tree_points_begin,cache_begin,dim);
         if  (q.outer_range_contains(b_upper))
           it=node->upper()->tree_items(it);
         else
           if (q.inner_range_intersects(b_upper))
             it=node->upper()->search(it,q,b_upper,tree_points_begin,cache_begin,dim);
      };
      return it;
    }


    template <class FuzzyQueryItem>
    boost::optional<Point_d>
    search_any_point(const FuzzyQueryItem& q,
                     Kd_tree_rectangle<FT,D>& b,
                     typename Kdt::const_iterator tree_points_begin,
                     typename std::vector<FT>::const_iterator cache_begin,
                     int dim) const
    {
      boost::optional<Point_d> result = boost::none;
      if (is_leaf()) {
        Leaf_node_const_handle node =
          static_cast<Leaf_node_const_handle>(this);
        if (node->size()>0)
        {
          typename internal::Has_points_cache<Kdt, internal::has_Enable_points_cache<Kdt>::type::value>::type dummy;
          result = search_any_point_in_leaf(node, q, tree_points_begin, cache_begin, dim, dummy);
        }
      }
      else {
         Internal_node_const_handle node =
          static_cast<Internal_node_const_handle>(this);
        // after splitting b denotes the lower part of b
        Kd_tree_rectangle<FT,D> b_upper(b);
        node->split_bbox(b, b_upper);

        if (q.inner_range_intersects(b)) {
          result = node->lower()->search_any_point(q,b,tree_points_begin,cache_begin,dim);
          if(result)
            return result;
        }
        if (q.inner_range_intersects(b_upper))
          result = node->upper()->search_any_point(q,b_upper,tree_points_begin,cache_begin,dim);
      }
      return result;
    }

  private:

    // If contains_point_given_as_coordinates does not exist in `FuzzyQueryItem`
    template <typename FuzzyQueryItem>
    bool contains(
      const FuzzyQueryItem& q,
      Point_d const& p,
      typename std::vector<FT>::const_iterator /*it_coord_begin*/,
      typename std::vector<FT>::const_iterator /*it_coord_end*/,
      Tag_false /*has_contains_point_given_as_coordinates*/) const
    {
      return q.contains(p);
    }
    // ... or if it exists
    template <typename FuzzyQueryItem>
    bool contains(
      const FuzzyQueryItem& q,
      Point_d const& /*p*/,
      typename std::vector<FT>::const_iterator it_coord_begin,
      typename std::vector<FT>::const_iterator it_coord_end,
      Tag_true /*has_contains_point_given_as_coordinates*/) const
    {
      return q.contains_point_given_as_coordinates(it_coord_begin, it_coord_end);
    }

    // With cache
    template<class FuzzyQueryItem, class OutputIterator>
    OutputIterator search_in_leaf(
      Leaf_node_const_handle node,
      const FuzzyQueryItem &q,
      typename Kdt::const_iterator tree_points_begin,
      typename std::vector<FT>::const_iterator cache_begin,
      int dim,
      OutputIterator oit,
      Tag_true /*has_points_cache*/) const
    {
      typename Kdt::iterator it_node_point = node->begin(), it_node_point_end = node->end();
      typename std::vector<FT>::const_iterator cache_point_it = cache_begin + dim*(it_node_point - tree_points_begin);
      for (; it_node_point != it_node_point_end; ++it_node_point, cache_point_it += dim)
      {
        Boolean_tag<has_contains_point_given_as_coordinates<FuzzyQueryItem>::value> dummy;
        if (contains(q, *it_node_point, cache_point_it, cache_point_it + dim, dummy))
          *oit++ = *it_node_point;
      }
      return oit;
    }

    // Without cache
    template<class FuzzyQueryItem, class OutputIterator>
    OutputIterator search_in_leaf(
      Leaf_node_const_handle node,
      const FuzzyQueryItem &q,
      typename Kdt::const_iterator /*tree_points_begin*/,
      typename std::vector<FT>::const_iterator /*cache_begin*/,
      int /*dim*/,
      OutputIterator oit,
      Tag_false /*has_points_cache*/) const
    {
      for (iterator i = node->begin(); i != node->end(); ++i)
      {
        if (q.contains(*i))
          *oit++ = *i;
      }
      return oit;
    }

    // With cache
    template<class FuzzyQueryItem>
    boost::optional<Point_d> search_any_point_in_leaf(
      Leaf_node_const_handle node,
      const FuzzyQueryItem &q,
      typename Kdt::const_iterator tree_points_begin,
      typename std::vector<FT>::const_iterator cache_begin,
      int dim,
      Tag_true /*has_points_cache*/) const
    {
      boost::optional<Point_d> result = boost::none;
      typename Kdt::iterator it_node_point = node->begin(), it_node_point_end = node->end();
      typename std::vector<FT>::const_iterator cache_point_it = cache_begin + dim*(it_node_point - tree_points_begin);
      for (; it_node_point != it_node_point_end; ++it_node_point, cache_point_it += dim)
      {
        Boolean_tag<has_contains_point_given_as_coordinates<FuzzyQueryItem>::value> dummy;
        if (contains(q, *it_node_point, cache_point_it, cache_point_it + dim, dummy))
        {
          result = *it_node_point;
          break;
        }
      }
      return result;
    }

    // Without cache
    template<class FuzzyQueryItem>
    boost::optional<Point_d> search_any_point_in_leaf(
      Leaf_node_const_handle node,
      const FuzzyQueryItem &q,
      typename Kdt::const_iterator /*tree_points_begin*/,
      typename std::vector<FT>::const_iterator /*cache_begin*/,
      int /*dim*/,
      Tag_false /*has_points_cache*/) const
    {
      boost::optional<Point_d> result = boost::none;
      for (iterator i = node->begin(); i != node->end(); ++i)
      {
        if (q.contains(*i))
        {
          result = *i;
          break;
        }
      }
      return result;
    }
  };


  template < class TreeTraits, class Splitter, class UseExtendedNode, class EnablePointsCache >
  class Kd_tree_leaf_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode, EnablePointsCache >{

    friend class Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache>;

    typedef typename Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache>::iterator iterator;
    typedef Kd_tree_node< TreeTraits, Splitter, UseExtendedNode, EnablePointsCache> Base;
    typedef typename TreeTraits::Point_d Point_d;

  private:

    // private variables for leaf nodes
    boost::int32_t n; // denotes number of items in a leaf node
    iterator data; // iterator to data in leaf node


  public:

    // default constructor
    Kd_tree_leaf_node()
      : Base(true)
    {}

    Kd_tree_leaf_node(unsigned int n_ )
      : Base(true), n(n_)
    {}

    // members for all nodes

    // members for leaf nodes only
    inline
    unsigned int
    size() const
    {
      return n;
    }

    inline
    iterator
    begin() const
    {
      return data;
    }

    inline
    iterator
    end() const
    {
      return data + n;
    }

    inline
    void
    drop_last_point()
    {
      --n;
    }

  }; //leaf node



  template < class TreeTraits, class Splitter, class UseExtendedNode, class EnablePointsCache>
  class Kd_tree_internal_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode, EnablePointsCache >{

    friend class Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache>;

    typedef Kd_tree<TreeTraits, Splitter, UseExtendedNode, EnablePointsCache> Kdt;

    typedef Kd_tree_node< TreeTraits, Splitter, UseExtendedNode, EnablePointsCache> Base;
    typedef typename Kdt::Node_handle Node_handle;
    typedef typename Kdt::Node_const_handle Node_const_handle;

    typedef typename TreeTraits::FT FT;
    typedef typename Kdt::Separator Separator;
    typedef typename Kdt::D D;

  private:

       // private variables for internal nodes
    boost::int32_t cut_dim;
    FT cut_val;
    Node_handle lower_ch, upper_ch;


    // private variables for extended internal nodes
    FT upper_low_val;
    FT upper_high_val;
    FT lower_low_val;
    FT lower_high_val;


  public:

    // default constructor
    Kd_tree_internal_node()
      : Base(false), cut_dim(-1), cut_val(0)
      , lower_ch (nullptr), upper_ch (nullptr)
      , upper_low_val(0), upper_high_val(0)
      , lower_low_val(0), lower_high_val(0)
    {}

    // members for internal node and extended internal node

    inline
    Node_const_handle
    lower() const
    {
      return lower_ch;
    }

    inline
    Node_const_handle
    upper() const
    {
      return upper_ch;
    }

    inline
    Node_handle
    lower()
    {
      return lower_ch;
    }

    inline
    Node_handle
    upper()
    {
      return upper_ch;
    }

    inline
    void
    set_lower(Node_handle nh)
    {
      lower_ch = nh;
    }

    inline
    void
    set_upper(Node_handle nh)
    {
      upper_ch = nh;
    }

    // inline Separator& separator() {return sep; }
    // use instead
    inline
    void set_separator(Separator& sep){
      cut_dim = sep.cutting_dimension();
      cut_val = sep.cutting_value();
    }

    inline
    FT
    cutting_value() const
    {
      return cut_val;
    }

    inline
    int
    cutting_dimension() const
    {
      return cut_dim;
    }

    // members for extended internal node only
    inline
    FT
    upper_low_value() const
    {
      return upper_low_val;
    }

    inline
    FT
    upper_high_value() const
    {
      return upper_high_val;
    }

    inline
    FT
    lower_low_value() const
    {
      return lower_low_val;
    }

    inline
    FT
    lower_high_value() const
    {
      return lower_high_val;
    }

    /*Separator&
    separator()
    {
      return Separator(cutting_dimension,cutting_value);
    }*/

    void split_bbox(Kd_tree_rectangle<FT,D>& l, Kd_tree_rectangle<FT,D>& u) const {
      l.lower()[cut_dim]=lower_low_val;
      l.upper()[cut_dim]=lower_high_val;
      u.lower()[cut_dim]=upper_low_val;
      u.upper()[cut_dim]=upper_high_val;
    }
  };//internal node

 template < class TreeTraits, class Splitter, class EnablePointsCache>
 class Kd_tree_internal_node<TreeTraits,Splitter,Tag_false,EnablePointsCache>
   : public Kd_tree_node< TreeTraits, Splitter, Tag_false, EnablePointsCache >
 {
    friend class Kd_tree<TreeTraits, Splitter, Tag_false, EnablePointsCache>;

    typedef Kd_tree<TreeTraits, Splitter, Tag_false, EnablePointsCache> Kdt;

    typedef Kd_tree_node< TreeTraits, Splitter, Tag_false, EnablePointsCache> Base;
    typedef typename Kdt::Node_handle Node_handle;
    typedef typename Kdt::Node_const_handle Node_const_handle;

    typedef typename TreeTraits::FT FT;
    typedef typename Kdt::Separator Separator;
    typedef typename Kdt::D D;

  private:

       // private variables for internal nodes
    boost::uint8_t cut_dim;
    FT cut_val;

    Node_handle lower_ch, upper_ch;

  public:

    // default constructor
    Kd_tree_internal_node()
      : Base(false)
    {}

    // members for internal node and extended internal node

    inline
    Node_const_handle
    lower() const
    {
      return lower_ch;
    }

    inline
    Node_const_handle
    upper() const
    {
      return upper_ch;
    }

    inline
    Node_handle
    lower()
    {
      return lower_ch;
    }

    inline
    Node_handle
    upper()
    {
      return upper_ch;
    }

    inline
    void
    set_lower(Node_handle nh)
    {
      lower_ch = nh;
    }

    inline
    void
    set_upper(Node_handle nh)
    {
      upper_ch = nh;
    }

    // inline Separator& separator() {return sep; }
    // use instead

    inline
    void set_separator(Separator& sep){
      cut_dim = static_cast<boost::uint8_t>(sep.cutting_dimension());
      cut_val = sep.cutting_value();
    }

    inline
    FT
    cutting_value() const
    {
      return cut_val;
    }

    inline
    int
    cutting_dimension() const
    {
      return cut_dim;
    }

   /* Separator&
    separator()
    {
      return Separator(cutting_dimension,cutting_value);
    }*/

    void split_bbox(Kd_tree_rectangle<FT,D>& l, Kd_tree_rectangle<FT,D>& u) const {
      l.upper()[cut_dim]=cut_val;
      u.lower()[cut_dim]=cut_val;
    }
  };//internal node



} // namespace CGAL
#endif // CGAL_KDTREE_NODE_H
