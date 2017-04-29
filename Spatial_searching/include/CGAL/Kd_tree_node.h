// Copyright (c) 2002,2011  Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H

#include <CGAL/license/Spatial_searching.h>



#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>
#include <boost/cstdint.hpp>

namespace CGAL {

  template <class SearchTraits, class Splitter, class UseExtendedNode> 
  class Kd_tree;

  template < class TreeTraits, class Splitter, class UseExtendedNode > 
  class Kd_tree_node {

     friend class Kd_tree<TreeTraits,Splitter,UseExtendedNode>;

    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Node_handle Node_handle;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Node_const_handle Node_const_handle;
     typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Internal_node_handle Internal_node_handle;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Internal_node_const_handle Internal_node_const_handle;
     typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Leaf_node_handle Leaf_node_handle;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Leaf_node_const_handle Leaf_node_const_handle;
    typedef typename TreeTraits::Point_d Point_d;

    typedef typename TreeTraits::FT FT;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Separator Separator;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Point_d_iterator Point_d_iterator;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::iterator iterator;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::D D;

    bool leaf;

  public : 
    Kd_tree_node(bool leaf_)
      :leaf(leaf_){}

    bool is_leaf() const{
      return leaf;
    }

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
	   Kd_tree_rectangle<FT,D>& b) const
    {
      if (is_leaf()) { 
        Leaf_node_const_handle node = 
          static_cast<Leaf_node_const_handle>(this);
	if (node->size()>0) 
	  for (iterator i=node->begin(); i != node->end(); i++) 
	    if (q.contains(*i)) 
	      {*it++=*i;}
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
	    it=node->lower()->search(it,q,b);
	if  (q.outer_range_contains(b_upper))     
	  it=node->upper()->tree_items(it);
	else
	  if (q.inner_range_intersects(b_upper)) 
	    it=node->upper()->search(it,q,b_upper);
      };
      return it;				
    }


    template <class FuzzyQueryItem>
    boost::optional<Point_d>
    search_any_point(const FuzzyQueryItem& q,
                     Kd_tree_rectangle<FT,D>& b) const
    {
      boost::optional<Point_d> result = boost::none;
      if (is_leaf()) { 
        Leaf_node_const_handle node = 
          static_cast<Leaf_node_const_handle>(this);
	if (node->size()>0) 
	  for (iterator i=node->begin(); i != node->end(); i++) 
	    if (q.contains(*i)) 
	      { result = *i; break; }
      }
      else {
         Internal_node_const_handle node = 
          static_cast<Internal_node_const_handle>(this);
	// after splitting b denotes the lower part of b
	Kd_tree_rectangle<FT,D> b_upper(b);
	node->split_bbox(b, b_upper);
                             
	if (q.outer_range_contains(b)){ 	
          result = node->lower()->any_tree_item();
	}else{
	  if (q.inner_range_intersects(b)){ 
	    result = node->lower()->search_any_point(q,b);
          }
        }
        if(result){
          return result;
        }
	if  (q.outer_range_contains(b_upper)){     
	  result = node->upper()->any_tree_item();
	}else{
	  if (q.inner_range_intersects(b_upper)) 
	    result = node->upper()->search_any_point(q,b_upper);
        }
      }
      return result;				
    }

  };


  template < class TreeTraits, class Splitter, class UseExtendedNode > 
  class Kd_tree_leaf_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode >{

    friend class Kd_tree<TreeTraits,Splitter,UseExtendedNode>;
    
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::iterator iterator;
    typedef Kd_tree_node< TreeTraits, Splitter, UseExtendedNode> Base;
    typedef typename TreeTraits::Point_d Point_d;

  private:
    
    // private variables for leaf nodes
    boost::int32_t n; // denotes number of items in a leaf node
    iterator data; // iterator to data in leaf node

                    
  public:
		
    // default constructor
    Kd_tree_leaf_node() 
    {}

    Kd_tree_leaf_node(bool leaf_ ) 
      : Base(leaf_)
    {}

    Kd_tree_leaf_node(bool leaf_,unsigned int n_ ) 
      : Base(leaf_), n(n_)
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



  template < class TreeTraits, class Splitter, class UseExtendedNode> 
  class Kd_tree_internal_node : public Kd_tree_node< TreeTraits, Splitter, UseExtendedNode >{

    friend class Kd_tree<TreeTraits,Splitter,UseExtendedNode>;

    typedef Kd_tree_node< TreeTraits, Splitter, UseExtendedNode> Base;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Node_handle Node_handle;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Node_const_handle Node_const_handle;

    typedef typename TreeTraits::FT FT;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Separator Separator;
    typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::D D;

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
    {}

    Kd_tree_internal_node(bool leaf_) 
      : Base(leaf_)
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

 template < class TreeTraits, class Splitter> 
 class Kd_tree_internal_node<TreeTraits,Splitter,Tag_false> : public Kd_tree_node< TreeTraits, Splitter, Tag_false >{

    friend class Kd_tree<TreeTraits,Splitter,Tag_false>;

    typedef Kd_tree_node< TreeTraits, Splitter, Tag_false> Base;
    typedef typename Kd_tree<TreeTraits,Splitter,Tag_false>::Node_handle Node_handle;
    typedef typename Kd_tree<TreeTraits,Splitter,Tag_false>::Node_const_handle Node_const_handle;

    typedef typename TreeTraits::FT FT;
    typedef typename Kd_tree<TreeTraits,Splitter,Tag_false>::Separator Separator;
    typedef typename Kd_tree<TreeTraits,Splitter,Tag_false>::D D;

  private:
    
       // private variables for internal nodes
    boost::uint8_t cut_dim;
    FT cut_val;

    Node_handle lower_ch, upper_ch;
                
  public:

    // default constructor
    Kd_tree_internal_node() 
    {}

    Kd_tree_internal_node(bool leaf_) 
      : Base(leaf_)
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
