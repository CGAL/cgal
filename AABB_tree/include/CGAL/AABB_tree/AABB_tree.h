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

//template <class Kernel, class Input, class PSC>
//class AABB_tree
//{
//public:
//
//  // basic kernel object types
//  typedef typename Kernel::FT FT;
//  typedef typename CGAL::Bbox_3 Bbox;
//  typedef typename Kernel::Ray_3 Ray;
//  typedef typename Kernel::Line_3 Line;
//  typedef typename Kernel::Plane_3 Plane;
//  typedef typename Kernel::Point_3 Point;
//  typedef typename Kernel::Vector_3 Vector;
//  typedef typename Kernel::Segment_3 Segment;
//  typedef typename Kernel::Triangle_3 Triangle;
//  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;
//
//  // Nodes of the tree
//  typedef AABB_node<Kernel,Input,PSC> Node;
//  typedef typename Node::Point_with_input Point_with_input;
//
//private:
//
//  // set of input primitives (halfedge or face handles)
//  std::vector<Input> m_data;
//
//  // single root node
//  Node *m_root;
//
//public:
//  // life cycle
//  AABB_tree() : m_root(NULL) {}
//  ~AABB_tree()
//  {
//    cleanup();
//  }
//
//  void cleanup()
//  {
//    m_data.clear();
//    delete [] m_root; m_root = NULL;
//  }
//
//  // build tree when input = face_handle
//  bool build_faces(const PSC& psc)
//  {
//    cleanup();
//    set_face_data(psc);
//    if(!empty())
//    {
//      m_root = new Node[m_data.size()-1]();
//      m_root->expand(m_data.begin(),m_data.end(), m_data.size());
//      return true;
//    }
//    return false;
//  }
//
//private:
//
//  void set_face_data(const PSC& psc)
//  {
//    unsigned int nbf = psc.size_of_facets();
//    m_data.reserve(nbf);
//    typename PSC::Facet_const_iterator f;
//    for(f = psc.facets_begin(); f != psc.facets_end(); ++f)
//      m_data.push_back(f);
//  }
//
//public:
//  Bbox_3 bbox() const
//  {
//    return m_root->bbox();
//  }
//
//  size_t size() const {
//    return m_data.size();
//  }
//
//  double max_bbox_length() const
//  {
//    return m_root->max_length();
//  }
//
//  bool empty()
//  {
//    return m_data.size() < 2; // TODO: change this requirement to < 1
//  }
//
//  // --------------------QUERY FUNCTIONS----------------------//
//
//  // generic traversal
//  template <class T, class Traits>
//  void traversal(const T& x, Traits& traits) const
//  {
//    m_root->template traversal<Traits,T>(x, traits, m_data.size());
//  }
//
//
//  template<class QueryType, class ResultType>
//  class First_intersection_traits
//  {
//  private:
//    ResultType& r;
//  public:
//    bool go_further() const
//    {
//      return !r.first;
//    }
//    First_intersection_traits(ResultType& result) : r(result) {}
//
//    bool intersection(const QueryType& q, const Input& i)
//    {
//      ResultType result;
//      if(Node::intersection(q, i, result.second))
//      {
//	r.first = true;
//	r.second = result.second;
//	return true;
//      }
//      return false;
//    }
//    bool do_intersect(const QueryType& q, const Node& node) const
//    {
//      return Node::do_intersect(q, node);
//    }
//  };
//
//  typedef boost::mpl::vector<Plane, Ray, Line, Segment> Allowed_query_types;
//
//  // The following function template is restricted to that T can only be in
//  // {Ray, Line, Segment}. It return type is bool.
//  // The trick uses enable_if and the Boost MPL.
//  template <class T>
//  typename boost::enable_if<
//    typename boost::mpl::contains<Allowed_query_types,
//                                  T>::type,
//    bool>::type
//  first_intersection(const T& x,
//                     Point_with_input& pwh)
//  {
//    typedef std::pair<bool,Point_with_input> Result_type;
//    typedef First_intersection_traits<T, Result_type> Traits;
//
//    Result_type result;
//    Traits traits(result);
//    m_root->template traversal<Traits,T>(x, traits, m_data.size());
//    if(result.first)
//    {
//      pwh = result.second;
//      return true;
//    }
//    return false;
//  }
//
//  template<class QueryType, class Output_iterator, class Value_type>
//  class Listing_traits
//  {
////     typedef typename std::iterator_traits<Output_iterator>::value_type
////     Pt;
//    typedef Value_type Pt;
//  private:
//    Output_iterator& out_it;
//  public:
//    bool go_further()
//    {
//      return true;
//    }
//    Listing_traits(Output_iterator& out_it_) : out_it(out_it_) {}
//    bool intersection(const QueryType& q, const Input& i)
//    {
//      Pt p;
//      if(Node::intersection(q, i, p))
//      {
//        *out_it++ = p;
//        return true;
//      }
//      return false;
//    }
//    bool do_intersect(const QueryType& q, const Node& node)
//    {
//      return Node::do_intersect(q, node);
//    }
//  }; // end class Listing_traits<QueryType,Container>
//
//  // The following function template is restricted to that T can only be in
//  // {Ray, Line, Segment}. It return type is Output_iterator&.
//  // The trick uses enable_if and the Boost MPL.
//  template <class Value_type, class T, typename Output_iterator>
//  typename boost::enable_if<
//    typename boost::mpl::contains<Allowed_query_types,
//                                  T>::type,
//    Output_iterator>::type
//  all_intersection(const T& x,
//                   Output_iterator out)
//  {
//    typedef Listing_traits<T, Output_iterator, Value_type> Traits;
//
//    Traits traits(out);
//    m_root->template traversal<Traits,T>(x, traits, m_data.size());
//    return out;
//  }
//
//  template <class Value_type>
//  class Counting_iterator {
//    typedef Counting_iterator<Value_type> Self;
//    int& i;
//  public:
//    Counting_iterator(int& i_) : i(i_) {};
//
//    struct Proxy {
//      Proxy& operator=(const Value_type&) { return *this; };
//    };
//
//    Proxy operator*() {
//      return Proxy();
//    }
//
//    Self& operator++() {
//      ++i;
//      return *this;
//    }
//
//    Self& operator++(int) {
//      ++i;
//      return *this;
//    }
//  };
//
//  // The following function template is restricted to that T can only be in
//  // {Ray, Line, Segment}. It return type is int.
//  // The trick uses enable_if and the Boost MPL.
//  template <class T>
//  typename boost::enable_if<
//    typename boost::mpl::contains<Allowed_query_types,
//                                  T>::type,
//    int>::type
//  count_intersections(const T& x)
//  {
//    typedef Listing_traits<T, Counting_iterator<Point>, Point> Traits;
//    int result = 0;
//    Counting_iterator<Point> counting_it(result);
//    Traits traits(counting_it);
//    m_root->template traversal<Traits,T>(x, traits, m_data.size());
//    return result ;
//  }
//}; // end class AABB_tree
//




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
   * iterator on a container of Primitive::Data_type such that Primitive has
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
