// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// Copyrigth (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     :  Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_AABB_TREE_H
#define CGAL_AABB_TREE_H

#include <vector>
#include <iterator>
#include <CGAL/AABB_tree/AABB_node.h>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/contains.hpp>

namespace CGAL {

/**
 * @class AABB_tree
 *
 *
 */
template<typename AABBTraits>
class AABB_tree
{
public:
  /// Traits types
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Projection_query Projection_query;
  typedef typename AABBTraits::Projection_return Projection_return;
  typedef typename AABBTraits::Intersection_type Intersection_type;
  typedef typename AABBTraits::Bounding_box Bounding_box;


  /**
   * @brief Constructor
   * @param first the first primitive to insert
   * @param last the last primitive to insert
   *
   * Builds the datastructure. Type ConstPrimitiveIterator can be any const
   * iterator on a container of Primitive::id_type such that Primitive has
   * a constructor taking a ConstPrimitiveIterator as argument.
   */
  template<typename ConstPrimitiveIterator>
  AABB_tree(ConstPrimitiveIterator first, ConstPrimitiveIterator last);

  /// Default constructor
  //AABB_tree() : data_(), p_root_(NULL) {};

  /// Non virtual destructor
  ~AABB_tree();

  template<typename Query>
  bool do_intersect(const Query& q) const;

  template<typename Query>
  int number_of_intersections(const Query& q) const;

  template<typename Query, typename OutputIterator>
  OutputIterator intersected_primitives(const Query& q,
                                        OutputIterator out) const;

  template<typename Query, typename OutputIterator>
  OutputIterator all_intersections(const Query& q,
                                   OutputIterator out) const;

  template<typename Query>
  bool any_intersection(const Query& q,
                        Intersection_type& intersection) const;

  Projection_return closest_point(const Projection_query& q,
                                  const Projection_return& hint) const;


  //////////////////////////////////////////////
  //TODO: document this
  Bounding_box root_bbox() const { return p_root_->bounding_box(); }
  bool isEmpty() const { return data_.empty(); }
  size_t size() const { return data_.size(); }

  /// generic traversal of tree
  template <class Query, class Traversal_traits>
  void traversal(const Query& q, Traversal_traits& traits) const
  {
    p_root_->template traversal<Traversal_traits,Query>(q, traits, data_.size());
  }
  //////////////////////////////////////////////

private:
  typedef AABB_node<AABBTraits> Node;
  typedef typename AABBTraits::Sphere Sphere;

  //-------------------------------------------------------
  // Traits classes for traversal computation
  //-------------------------------------------------------
  /**
   * @class First_intersection_traits
   */
  template<typename Query>
  class First_intersection_traits
  {
  public:
    First_intersection_traits()
      : is_found_(false)
      , result_() {}

    bool go_further() const { return !is_found_; }

    void intersection(const Query& q, const Primitive& primitive)
    {
      is_found_ = AABBTraits().intersection(q, primitive, result_);
    }

    bool do_intersect(const Query& q, const Node& node) const
    {
      return AABBTraits().do_intersect(q, node.bounding_box());
    }

    Intersection_type result() const { return result_; }
    bool is_intersection_found() const { return is_found_; }

  private:
    bool is_found_;
    Intersection_type result_;
  };


  /**
   * @class Counting_traits
   */
  template<typename Query>
  class Counting_traits
  {
  public:
    Counting_traits()
      : intersection_()
      , intersection_nb_(0) {}

    bool go_further() const { return true; }

    void intersection(const Query& q, const Primitive& primitive)
    {
      if( AABBTraits().intersection(q, primitive, intersection_) )
      {
        ++intersection_nb_;
      }
    }

    bool do_intersect(const Query& q, const Node& node) const
    {
      return AABBTraits().do_intersect(q, node.bounding_box());
    }

    int intersection_number() const { return intersection_nb_; }

  private:
    Intersection_type intersection_;
    int intersection_nb_;
  };


  /**
   * @class Listing_intersection_traits
   */
  template<typename Query, typename Output_iterator>
  class Listing_intersection_traits
  {
  public:
    Listing_intersection_traits(Output_iterator out_it)
      : intersection_()
      , out_it_(out_it) {}

    bool go_further() const { return true; }

    void intersection(const Query& q, const Primitive& primitive)
    {
      if( AABBTraits().intersection(q, primitive, intersection_) )
      {
        *out_it_++ = intersection_;
      }
    }

    bool do_intersect(const Query& q, const Node& node) const
    {
      return AABBTraits().do_intersect(q, node.bounding_box());
    }

  private:
    Intersection_type intersection_;
    Output_iterator out_it_;
  };


  /**
   * @class Listing_primitive_traits
   */
  template<typename Query, typename Output_iterator>
  class Listing_primitive_traits
  {
  public:
    Listing_primitive_traits(Output_iterator out_it)
      : intersection_()
      , out_it_(out_it) {}

    bool go_further() const { return true; }

    void intersection(const Query& q, const Primitive& primitive)
    {
      if( AABBTraits().intersection(q, primitive, intersection_) )
      {
        *out_it_++ = primitive;
      }
    }

    bool do_intersect(const Query& q, const Node& node) const
    {
      return AABBTraits().do_intersect(q, node.bounding_box());
    }

  private:
    Intersection_type intersection_;
    Output_iterator out_it_;
  };

  /**
   * @class Projection_traits
   */
  class Projecting_traits
  {
  public:
    Projecting_traits(const Projection_query& query,
                      const Projection_return& hint)
      : projection_(hint)
      , center_(query)
      , sphere_(AABBTraits().sphere(query,hint))         { }

    bool go_further() const { return true; }

    void intersection(const Projection_query& q, const Primitive& primitive)
    {
      // We don't use q here because it is embedded in sphere_ and we don't
      // want to compute sphere everytime

      Projection_return projection;
      if ( AABBTraits().intersection(sphere_, primitive, projection) )
      {
        const Sphere sphere = AABBTraits().sphere(center_, projection);
        if ( AABBTraits().is_smaller(sphere, sphere_) )
        {
          projection_ = projection;
          sphere_ = sphere;
        }
      }
    }

    bool do_intersect(const Projection_query& q, const Node& node) const
    {
      return AABBTraits().do_intersect(sphere_, node.bounding_box());
    }

    Projection_return projection() const { return projection_; }

  private:
    Projection_return projection_;
    Projection_query center_;
    Sphere sphere_;
  };


private:



private:
  // set of input primitives (halfedge or face handles)
  std::vector<Primitive> data_;
  // single root node
  Node* p_root_;


private:
  // Disabled copy constructor & assignment operator
  typedef AABB_tree<AABBTraits> Self;
  AABB_tree(const Self& src);
  Self& operator=(const Self& src);

};  // end class AABB_tree



template<typename Tr>
template<typename ConstPrimitiveIterator>
AABB_tree<Tr>::AABB_tree(ConstPrimitiveIterator first,
                         ConstPrimitiveIterator last)
: data_()
, p_root_(NULL)
{
  // Insert each primitive into tree
  // TODO: get number of elements to reserve space ?
  while ( first != last )
  {
    data_.push_back(Primitive(first));
    ++first;
  }

  p_root_ = new Node[data_.size()-1]();
  p_root_->expand(data_.begin(), data_.end(), data_.size());
}


template<typename Tr>
AABB_tree<Tr>::~AABB_tree()
{
  delete[] p_root_;
}


template<typename Tr>
template<typename Query>
bool
AABB_tree<Tr>::do_intersect(const Query& query) const
{
  typedef First_intersection_traits<Query> Traversal_traits;
  Traversal_traits traversal_traits;

  this->traversal(query, traversal_traits);
  return traversal_traits.is_intersection_found();
}



template<typename Tr>
template<typename Query>
int
AABB_tree<Tr>::number_of_intersections(const Query& query) const
{
  typedef Counting_traits<Query> Traversal_traits;
  Traversal_traits traversal_traits;

  this->traversal(query, traversal_traits);
  return traversal_traits.intersection_number();
}


template<typename Tr>
template<typename Query, typename OutputIterator>
OutputIterator
AABB_tree<Tr>::intersected_primitives(const Query& query,
                                      OutputIterator out) const
{
  typedef Listing_primitive_traits<Query, OutputIterator> Traversal_traits;
  Traversal_traits traversal_traits(out);

  this->traversal(query, traversal_traits);
  return out;
}



template<typename Tr>
template<typename Query, typename OutputIterator>
OutputIterator
AABB_tree<Tr>::all_intersections(const Query& query,
                                 OutputIterator out) const
{
  typedef Listing_intersection_traits<Query, OutputIterator> Traversal_traits;
  Traversal_traits traversal_traits(out);

  this->traversal(query, traversal_traits);
  return out;
}


template<typename Tr>
template<typename Query>
bool
AABB_tree<Tr>::any_intersection(const Query& query,
                                Intersection_type& intersection) const
{
  typedef First_intersection_traits<Query> Traversal_traits;
  Traversal_traits traversal_traits;

  this->traversal(query, traversal_traits);

  intersection = traversal_traits.result();
  return traversal_traits.is_intersection_found();
}


template<typename Tr>
typename AABB_tree<Tr>::Projection_return
AABB_tree<Tr>::closest_point(const Projection_query& query,
                             const Projection_return& hint) const
{
  Projecting_traits traversal_traits(query,hint);

  this->traversal(query, traversal_traits);
  return traversal_traits.projection();
}


} // end namespace CGAL

#endif // CGAL_AABB_TREE_H
