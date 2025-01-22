// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Olivier Billet, Mariette Yvinec

#ifndef CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
#define CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/basic.h>
#include <utility>
#include <map>
#include <set>
#include <list>

#include <boost/stl_interfaces/iterator_interface.hpp>

#include <CGAL/unordered_flat_map.h>
#include <CGAL/Skiplist.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/assertions.h>
#include <CGAL/Has_timestamp.h>
#include <CGAL/IO/io.h>

#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
#  define CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2 CGAL_CDT_2_DEBUG_INTERSECTIONS
#endif
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
#  include <CGAL/IO/io.h>
#  include <CGAL/Constrained_triangulation_2.h>
#  include <iostream>
#endif

namespace CGAL {

// T               is expected to be Vertex_handle
// Compare         is a comparison operator for type T
// Point           the point type of vertices
template <class T, class Compare, class Point>
class Polyline_constraint_hierarchy_2
{
public:
  typedef T                                       Vertex_handle;
  typedef std::pair<T, T>                         Subconstraint;

  using size_type = std::size_t;

private:
  class Node {
  public:
    explicit Node(Vertex_handle vh, bool input = false)
      : vertex_(vh), point_(vertex_->point()), input_(input)
    {}
    const Point& point() const { return point_; }
    Vertex_handle vertex() const { return vertex_; }
    bool& input() { return input_; }
    const bool& input() const { return input_; }
  private:
    Vertex_handle vertex_;
    Point point_;
    bool input_;
  };

  typedef CGAL::Skiplist<Node>  Vertex_list;
  typedef Vertex_list* Vertex_list_ptr;

public:
  // the base line is always
  class Point_it
    : public boost::iterator_adaptor<
    Point_it
    , typename Vertex_list::all_iterator
    , const Point&
    >
  {
  public:
    Point_it() : Point_it::iterator_adaptor_() {}
    Point_it(typename Vertex_list::all_iterator it) : Point_it::iterator_adaptor_(it) {}
  private:
    friend class boost::iterator_core_access;
    const Point& dereference() const { return this->base()->point(); }
  };

  // only nodes with a vertex_handle that is still in the triangulation
  class Vertex_it
    : public boost::iterator_adaptor<
    Vertex_it
    , typename Vertex_list::skip_iterator
    , Vertex_handle
    , std::bidirectional_iterator_tag
    , Vertex_handle>
  {
    typedef typename Vertex_list::skip_iterator Base_it;
  public:
    Vertex_it() : Vertex_it::iterator_adaptor_() {}
    Vertex_it(Base_it it) : Vertex_it::iterator_adaptor_(it) {}
    operator Point_it() const { return Point_it(this->base()); }
    bool& input() { return this->base()->input(); }
    const bool& input() const { return this->base()->input(); }

    friend bool operator==(const Vertex_it& a, const Base_it& it) { return a.base() == it; }
    friend bool operator==(const Base_it& it, const Vertex_it& a) { return a.base() == it; }
    friend bool operator!=(const Vertex_it& a, const Base_it& it) { return a.base() != it; }
    friend bool operator!=(const Base_it& it, const Vertex_it& a) { return a.base() != it; }

  private:
    friend class boost::iterator_core_access;
    Vertex_handle dereference() const { return this->base()->vertex(); }
  };

  struct Vertex_list_with_info {
    const Polyline_constraint_hierarchy_2* hierarchy_ptr = nullptr;
    Vertex_list vl{};
    bool may_share_subconstraint_with_others = false;
  };

  struct Constraint_id
  {
    Vertex_list_with_info* vl_with_info_ptr = nullptr;
    size_type id = (std::numeric_limits<size_type>::max)();

    Constraint_id(std::nullptr_t = nullptr) {}
    Constraint_id(Vertex_list_with_info* ptr, size_type id) : vl_with_info_ptr(ptr), id(id) {}

    void destroy() {
      if(vl_with_info_ptr != nullptr) {
        delete vl_with_info_ptr;
        (*this) = nullptr;
      }
    }

    Vertex_list_ptr vl_ptr() const {
      return vl_with_info_ptr == nullptr ? nullptr : std::addressof(vl_with_info_ptr->vl);
    }

    operator std::pair<Subconstraint, Vertex_list_ptr>() const {
      Subconstraint subconstraint = vl_with_info_ptr == nullptr
                                        ? Subconstraint()
                                        : Subconstraint(vl_ptr()->front().vertex(), vl_ptr()->back().vertex());
      return {subconstraint, vl_ptr()};
    }

    bool same_hierarchy(const Constraint_id& other) const {
      return vl_with_info_ptr == nullptr || other.vl_with_info_ptr == nullptr ||
             vl_with_info_ptr->hierarchy_ptr == other.vl_with_info_ptr->hierarchy_ptr;
    }

    Constraint_id& operator=(std::nullptr_t) {
      vl_with_info_ptr = nullptr;
      id = (std::numeric_limits<size_type>::max)();
      return *this;
    }
    bool operator==(std::nullptr_t n) const { return vl_with_info_ptr == n; }
    bool operator!=(std::nullptr_t n) const { return vl_with_info_ptr != n; }

    bool operator==(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (id == other.id));
      return vl_ptr() == other.vl_ptr();
    }

    bool operator!=(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (id == other.id));
      return vl_ptr() != other.vl_ptr();
    }

    bool operator<(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (id == other.id));
      return id < other.id;
    }

    // forward a new Vertex_list operations
    decltype(auto) begin() const { return vl_ptr()->skip_begin(); }
    decltype(auto) end() const { return vl_ptr()->skip_end(); }
    decltype(auto) elements() const { return vl_ptr()->skip_elements(); }
    decltype(auto) clear() const { return vl_ptr()->clear(); }
    decltype(auto) size() const { return vl_ptr()->skip_size(); }
    decltype(auto) front() const { return vl_ptr()->front(); }
    decltype(auto) back() const { return vl_ptr()->back(); }

  }; // end Constraint_id

  class Pair_compare {
    Compare comp;

  public:
    Pair_compare(const Compare& comp) : comp(comp) {}

    bool operator()(const Subconstraint& sc1, const Subconstraint& sc2) const {
      if(comp(sc1.first, sc2.first)) {
        return true;
      } else if((! comp(sc2.first, sc1.first)) && //  !less(sc1,sc2) && !less(sc2,sc1) == equal
                comp(sc1.second, sc2.second)) {
        return true;
      } else {
        return false;
      }
    }
  };

  class Context {
    friend class Polyline_constraint_hierarchy_2<T,Compare,Point>;
  private:
    Constraint_id enclosing;
    Vertex_it     pos;
  public:
    Context() : enclosing(nullptr) {}

    Vertex_it    vertices_begin()const { return enclosing.begin();}
    Vertex_it    current()const {return pos;}
    Vertex_it    vertices_end()const {return enclosing.end();}
    Constraint_id  id()const { return enclosing; }
    size_type    number_of_vertices() const {return enclosing.size(); }
  };

  typedef std::list<Context>              Context_list;
  typedef typename Context_list::iterator Context_iterator;

  static void fix_contexts(Context_list& context_list) {
    const bool multiple_contexts = context_list.size() > 1;
    if(false == multiple_contexts) return;
    for(auto& context : context_list) {
      CGAL_assertion(context.enclosing != nullptr);
      context.enclosing.vl_with_info_ptr->may_share_subconstraint_with_others = true;
    }
  }

  typedef std::set<Constraint_id>           Constraints_set;
#if CGAL_USE_BARE_STD_MAP
  typedef std::map<Subconstraint, Context_list*,
                   Pair_compare>            Sc_to_c_map;
#else
  typedef CGAL::unordered_flat_map<Subconstraint, Context_list*,
                                   boost::hash<Subconstraint>> Sc_to_c_map;
#endif
  typedef typename Constraints_set::iterator Constraint_iterator;
  typedef const Constraints_set& Constraints;
  typedef typename Sc_to_c_map::const_iterator    Sc_iterator;
  typedef typename Sc_to_c_map::iterator Sc_it;
  typedef Sc_iterator Subconstraint_and_contexts_iterator;
  typedef const Sc_to_c_map& Subconstraints_and_contexts;

  class Subconstraint_iterator : public boost::stl_interfaces::proxy_iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
                            Subconstraint_iterator,
#endif
                            std::bidirectional_iterator_tag,
                            Subconstraint>
  {
    using base_type = boost::stl_interfaces::proxy_iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
        Subconstraint_iterator,
#endif
        std::bidirectional_iterator_tag,
        Subconstraint>;

    const Polyline_constraint_hierarchy_2* hierarchy = nullptr;
    Constraint_iterator constraint_it{};
    Vertex_it vertex_it{};

  public:
    // - The object is singular if and only if `hierarchy==nullptr`.
    //
    // - The end value is when `constraint_it` is the end iterator of `constraints_set`.
    //   In that case `vertex_it` must be singular.
    //
    // - Otherwise all members must be valid pointers or dereferencable iterators.

    bool is_singular() const {
      CGAL_assertion((hierarchy == nullptr) == (constraint_it == Constraint_iterator{}));
      CGAL_assertion((hierarchy != nullptr) || (vertex_it == Vertex_it{}));
      return hierarchy == nullptr;
    }

    bool is_valid() const {
      if(hierarchy == nullptr) {
        CGAL_assertion(is_singular());
        return false;
      }
      if(constraint_it == Constraint_iterator{}) {
        return false;
      }
      if(constraint_it == hierarchy->constraints_set.end()) {
        return vertex_it == Vertex_it{};
      }
      if(vertex_it == Vertex_it{}) {
        return false;
      }
      return true;
    }

    bool is_end() const {
      return constraint_it == hierarchy->constraints_set.end();
    }

    bool is_begin() const {
      return constraint_it == hierarchy->constraints_set.begin() &&
             vertex_it == begin_or_null(constraint_it);
    }

    bool is_dereferenceable() const {
      return is_valid() &&
          constraint_it != hierarchy->constraints_set.end() &&
          vertex_it != constraint_it->end() &&
          std::next(vertex_it) != constraint_it->end();
    }

    bool equal(const Subconstraint_iterator& other) const {
      if(hierarchy != other.hierarchy) return false;
      if(is_singular()) return true;
      return (constraint_it == other.constraint_it &&
                  vertex_it == other.vertex_it);
    }

    Vertex_it begin_or_null(Constraint_iterator constraint_it) const {
      if(constraint_it == hierarchy->constraints_set.end()) {
        return Vertex_it();
      }
      return constraint_it->begin();
    }

    bool already_seen() const {
      if(is_end()) {
        return false;
      }
      if(false == constraint_it->vl_with_info_ptr->may_share_subconstraint_with_others) {
        return false;
      }
      auto [va, vb] = this->operator*();
      auto it = hierarchy->sc_to_c_map.find(hierarchy->sorted_pair(va, vb));
      CGAL_assertion(it != hierarchy->sc_to_c_map.end());

      const Context_list& ctx_list = *it->second;

      const Context& ctx = ctx_list.front();
      // if this context does not correspond to *this, return true

      if(ctx.enclosing != *constraint_it) {
        return true;
      }

      return (ctx.pos != vertex_it && ctx.pos != std::next(vertex_it));
    }

    static auto non_null(const Polyline_constraint_hierarchy_2* hierarchy) {
      CGAL_precondition(hierarchy != nullptr);
      return hierarchy;
    }

  public:
    Subconstraint_iterator() = default;

    // Constructors for begin and end. The constructors are public, but only the
    // hierarchy can create an iterator of this class, through its friendship of
    // the nested class `Construction_access`: Construction_access::begin_tag() and
    // Construction_access::end_tag().
    class Construction_access
    {
    private:
      friend class Subconstraint_iterator;
      friend class Polyline_constraint_hierarchy_2<T, Compare, Point>;

      static auto begin_tag() { return Begin_tag(); }
      static auto end_tag() { return End_tag(); }

      struct Begin_tag
      {};
      struct End_tag
      {};
    };
    //
    // constructor for the begin iterator
    explicit Subconstraint_iterator(typename Construction_access::Begin_tag,
                                    const Polyline_constraint_hierarchy_2* hierarchy)
        : hierarchy(non_null(hierarchy))
        , constraint_it(hierarchy->constraints_set.begin())
        , vertex_it(begin_or_null(hierarchy->constraints_set.begin()))
    {
      if(already_seen()) {
        ++(*this);
      }
    }
    //
    // constructor for the end iterator
    explicit Subconstraint_iterator(typename Construction_access::End_tag,
                                    const Polyline_constraint_hierarchy_2* hierarchy)
        : hierarchy(non_null(hierarchy))
        , constraint_it(hierarchy->constraints_set.end())
        , vertex_it() {}

    Subconstraint operator*() const {
      CGAL_precondition(is_dereferenceable());
      return Subconstraint(*vertex_it, *std::next(vertex_it));
    }

    friend bool operator==(const Subconstraint_iterator& lhs, const Subconstraint_iterator& rhs) {
      return lhs.equal(rhs);
    }

    using base_type::operator++;
    Subconstraint_iterator& operator++() {
      CGAL_precondition(is_valid() && false == is_end());

      do {
        ++vertex_it;
        CGAL_assertion(vertex_it != constraint_it->end());

        if(std::next(vertex_it) == constraint_it->end()) {
          ++constraint_it;
          vertex_it = begin_or_null(constraint_it);
        }
      } while(already_seen());
      return *this;
    }

    using base_type::operator--;
    Subconstraint_iterator& operator--() {
      CGAL_precondition(is_valid() && false == is_begin());

      do {
        if(constraint_it == hierarchy->constraints_set.end() || vertex_it == constraint_it->begin()) {
          --constraint_it;
          vertex_it = std::prev(constraint_it->end(), 2);
        } else {
          --vertex_it;
        }
      } while(already_seen());
      return *this;
    }
  }; // end class Subconstraint_iterator
  typedef Iterator_range<Subconstraint_iterator> Subconstraints;

private:
  Compare          comp;
  Constraints_set  constraints_set;
  Sc_to_c_map      sc_to_c_map;

public:
  Polyline_constraint_hierarchy_2(const Compare& comp)
    : comp(comp)
#if CGAL_USE_BARE_STD_MAP
    , sc_to_c_map(Pair_compare(comp))
#else
    , sc_to_c_map()
#endif
  { }
  Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch);
  Polyline_constraint_hierarchy_2(Polyline_constraint_hierarchy_2&&) = default;
  ~Polyline_constraint_hierarchy_2(){ clear();}
  void clear();
  Polyline_constraint_hierarchy_2& operator=(const Polyline_constraint_hierarchy_2& ch);
  Polyline_constraint_hierarchy_2& operator=(Polyline_constraint_hierarchy_2&& ch) = default;

  // Query
  bool is_subconstraint(T va, T vb) const;

  Vertex_it vertices_in_constraint_begin(Constraint_id cid) const
  { return cid.begin(); }
  Vertex_it vertices_in_constraint_end(Constraint_id cid) const
  { return cid.end(); }
  auto vertices_in_constraint(Constraint_id cid) const
  { return Iterator_range<Vertex_it>(cid.begin(), cid.end()); }


  Point_it points_in_constraint_begin(Constraint_id cid) const
  { return cid.vl_ptr()->all_begin(); }
  Point_it points_in_constraint_end(Constraint_id cid) const
  { return cid.vl_ptr()->all_end(); }

  bool enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const;
  bool next_along_sc(T va, T vb, T& w) const;
  void oriented_end(T va, T vb, T& vc) const;

  Context context(T va, T vb);
  size_type number_of_enclosing_constraints(T va, T vb) const;
  Context_iterator contexts_begin(T va, T vb) const;
  Context_iterator contexts_end(T va, T vb) const;
  Iterator_range<Context_iterator> contexts_range(T va, T vb) const;
  size_type number_of_constraints() const  { return constraints_set.size();}
  size_type number_of_subconstraints()const {return sc_to_c_map.size();}


  // insert/remove
  void add_Steiner(const T va, const T vb, const T vx);
  Constraint_id insert_constraint_old_API(T va, T vb);
  Constraint_id insert_constraint(T va, T vb);
  void append_constraint(Constraint_id cid, T va, T vb);
  void swap(Constraint_id first, Constraint_id second);
  void remove_constraint(Constraint_id cid);

  void split_constraint(T va, T vb, T vc);

  void simplify(Vertex_it u,
                Vertex_it v,
                Vertex_it w);

  size_type remove_points_without_corresponding_vertex(Constraint_id);
  size_type remove_points_without_corresponding_vertex();

  Constraint_id concatenate(Constraint_id first, Constraint_id second);
  Constraint_id concatenate2(Constraint_id first, Constraint_id second);
  Constraint_id split(Constraint_id first, Vertex_it vcit);
  Constraint_id split2(Constraint_id first, Vertex_it vcit);

  // iterators

  Subconstraint_and_contexts_iterator subconstraints_and_contexts_begin() const
  {
    return sc_to_c_map.begin();
  }

  Subconstraint_and_contexts_iterator subconstraints_and_contexts_end() const
  {
    return sc_to_c_map.end();
  }

  Subconstraint_iterator subconstraints_begin() const {
    BOOST_STL_INTERFACES_STATIC_ASSERT_CONCEPT(Subconstraint_iterator, std::bidirectional_iterator);
    BOOST_STL_INTERFACES_STATIC_ASSERT_ITERATOR_TRAITS(
      Subconstraint_iterator, std::bidirectional_iterator_tag, std::bidirectional_iterator,
      Subconstraint, Subconstraint, typename Subconstraint_iterator::pointer, std::ptrdiff_t);
    return Subconstraint_iterator(Subconstraint_iterator::Construction_access::begin_tag(),
                                  this);
  }

  Subconstraint_iterator subconstraints_end() const {
    return Subconstraint_iterator(Subconstraint_iterator::Construction_access::end_tag(),
                                  this);
  }

  Sc_iterator sc_begin() const{ return sc_to_c_map.begin(); }
  Sc_iterator sc_end()   const{ return sc_to_c_map.end();   }
  Constraint_iterator  constraints_begin()  const{ return constraints_set.begin(); }
  Constraint_iterator  constraints_end()    const{ return constraints_set.end();   }

  // Ranges
  const auto& constraints() const { return constraints_set; }
  const auto& subconstraints_and_contexts() const { return sc_to_c_map; }
  auto subconstraints() const {
    return Iterator_range<Subconstraint_iterator>(subconstraints_begin(), subconstraints_end());
  }

  // Helper functions
  void copy(const Polyline_constraint_hierarchy_2& ch);
  void copy(const Polyline_constraint_hierarchy_2& ch, std::map<Vertex_handle,Vertex_handle>& vmap);
  void swap(Polyline_constraint_hierarchy_2& ch);

private:
  template <typename F>
  void for_context_lists_of_all_subconstraints(Constraint_id cid, const F& f)
  {
    auto vl = cid.vl_ptr();
    for(Vertex_it it = vl->skip_begin(), succ = it, end = vl->skip_end(); ++succ != end; ++it) {
      auto scit = sc_to_c_map.find(sorted_pair(*it, *succ));
      CGAL_assertion(scit != sc_to_c_map.end());
      Context_list* hcl = scit->second;
      f(hcl, it, scit);
    }
  }

  static void replace_first_in_context_list(Context_list* hcl, Constraint_id old_id, Constraint_id new_id)
  {
    // std::find_if is a sort of std::for_each with a break
    [[maybe_unused]] auto it = std::find_if(hcl->begin(), hcl->end(), [&](Context& ctxt) {
      if(ctxt.enclosing == old_id) {
        ctxt.enclosing = new_id;
        return true;
      }
      return false;
    });
  }

  static void update_first_context_position(Context_list* hcl, Constraint_id id, Vertex_it new_pos)
  {
    [[maybe_unused]] auto it = std::find_if(hcl->begin(), hcl->end(), [&](Context& ctxt) {
      if(ctxt.enclosing == id) {
        ctxt.pos = new_pos;
        return true;
      }
      return false;
    });
  }

  static void remove_first_in_context_list(Context_list* hcl, Constraint_id id)
  {
    auto it = std::find_if(hcl->begin(), hcl->end(), [&](Context& ctxt) {
      return ctxt.enclosing == id;
    });
    if(it != hcl->end()) {
      hcl->erase(it);
    }
  }

  Constraint_id new_constraint_id() const {
    // TODO: handle ids
    auto id = number_of_constraints() == 0 ? 0 : constraints_set.rbegin()->id + 1;
    return Constraint_id(new Vertex_list_with_info{this}, id);
  }
  Subconstraint sorted_pair(T va, T vb) const;
  Subconstraint sorted_pair(Subconstraint sc) {
    const auto& [va, vb] = sc; return sorted_pair(va, vb);
  }
  Vertex_it get_pos(T va, T vb) const;
  bool      get_contexts(T va, T vb,
                         Context_iterator& ctxt,
                         Context_iterator& past) const;

  Context_list* get_context_list(T va, T vb) const;

  //to_debug
public:
  void   print(std::ostream& os = std::cout) const;
};

template <class T, class Compare, class Point>
Polyline_constraint_hierarchy_2<T,Compare,Point>::
Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch)
  : comp(ch.comp)
  , sc_to_c_map()
{
  copy(ch);
}

template <class T, class Compare, class Point>
Polyline_constraint_hierarchy_2<T,Compare,Point>&
Polyline_constraint_hierarchy_2<T,Compare,Point>::
operator=(const Polyline_constraint_hierarchy_2& ch){
  copy(ch);
  return *this;
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
copy(const Polyline_constraint_hierarchy_2& other)
{
  // create a identity transfer vertex map
  std::map<Vertex_handle, Vertex_handle>  vmap;
  for(const auto& cid : other.constraints()) {
    for(const auto& node : cid.elements()) {
      auto v = node.vertex();
      vmap[v] = v;
    }
  }
  copy(other, vmap);
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
copy(const Polyline_constraint_hierarchy_2& other, std::map<Vertex_handle,Vertex_handle>& vmap)
  // copy with a transfer vertex map
{
  std::map<Constraint_id, Constraint_id> cstr_map;
  clear();
  // copy constraints_set
  for(const auto& cid1: other.constraints()) {
    Constraint_id cid2 = new_constraint_id();
    cstr_map[cid1] = cid2;
    for(const auto& node : cid1.elements()) {
      cid2.vl_ptr()->push_back(Node(vmap[node.vertex()], node.input()));
    }
    cid2.vl_with_info_ptr->may_share_subconstraint_with_others =
        cid1.vl_with_info_ptr->may_share_subconstraint_with_others;
    constraints_set.insert(cid2);
  }
  // copy sc_to_c_map
  for(const auto& [sc1, hcl1] : other.subconstraints_and_contexts()) {
    Context_list* hcl2 = new Context_list;
    Vertex_handle uu2 = vmap[sc1.first];
    Vertex_handle vv2 = vmap[sc1.second];
    Subconstraint sc2 = sorted_pair(uu2, vv2);
    sc_to_c_map[sc2] = hcl2;
    for(const Context& ctxt1 : *hcl1) {
      // vertices of the enclosing constraints
      Context ctxt2;
      ctxt2.enclosing = cstr_map[ctxt1.enclosing];
      ctxt2.pos = ctxt2.enclosing.begin();
      Vertex_it aux = ctxt1.enclosing.begin();
      while(aux != ctxt1.pos) {
        ++aux;
        ++ctxt2.pos;
      }
      hcl2->push_back(ctxt2);
    }
  }

  comp = other.comp;
  return;
}


template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
swap(Polyline_constraint_hierarchy_2& ch)
{
  using std::swap;
  swap(comp, ch.comp);
  constraints_set.swap(ch.constraints_set);
  sc_to_c_map.swap(ch.sc_to_c_map);
}


template <class T, class Compare, class Point>
bool Polyline_constraint_hierarchy_2<T,Compare,Point>::
is_subconstraint(T va, T vb) const
{
  return( sc_to_c_map.find(sorted_pair(va, vb)) != sc_to_c_map.end() );
}


// used by Constrained_triangulation_plus_2::intersect with Exact_intersection_tag
template <class T, class Compare, class Point>
bool Polyline_constraint_hierarchy_2<T,Compare,Point>::
enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const
{
  Context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  // va = hcit->enclosing.front().vertex();
  // vb = hcit->enclosing.back().vertex();
  // Vertex_list_ptr vl = hcit->enclosing;
  Vertex_it pos = hcit->pos;
  if(vaa != *pos){
    std::swap(vaa,vbb);
  }
  while(!pos.input()){
    --pos;
  }
  va = *pos;
  pos = hcit->pos;
  ++pos;
  CGAL_assertion(vbb == *pos);
  while(!pos.input()){
    ++pos;
  }
  vb = *pos;
  return true;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Context
Polyline_constraint_hierarchy_2<T,Compare,Point>::
context(T va, T vb)
{
  Context_iterator hcit, past;
  if(!get_contexts(va,vb, hcit ,past)) CGAL_assertion(false);
  return *hcit;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::size_type
Polyline_constraint_hierarchy_2<T,Compare,Point>::
number_of_enclosing_constraints(T va, T vb) const
{
  Context_list* hcl = get_context_list(va,vb);
  CGAL_assertion(nullptr != hcl);
  return hcl->size();
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Context_iterator
Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts_begin(T va, T vb) const
{
   Context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_assertion(false);
   return first;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Context_iterator
Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts_end(T va, T vb) const
{
   Context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_assertion(false);
   return last;
}

template <class T, class Compare, class Point>
auto
Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts_range(T va, T vb) const -> Iterator_range<Context_iterator> {
  Context_iterator first, last;
  if( !get_contexts(va,vb,first,last)) return { first, first };
  else return { first, last };
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
swap(Constraint_id constr_a, Constraint_id constr_b) {
  auto substitute_enclosing_in_vertex_list = [this](Constraint_id cid, Constraint_id old_id, Constraint_id new_id) {
    for_context_lists_of_all_subconstraints(cid, [&](Context_list* hcl, Vertex_it, Sc_it) {
      replace_first_in_context_list(hcl, old_id, new_id);
    });
  };

  Vertex_list_ptr constr_a_vl = constr_a.vl_ptr();
  Vertex_list_ptr constr_b_vl = constr_b.vl_ptr();

  substitute_enclosing_in_vertex_list(constr_a, constr_a, nullptr);
  substitute_enclosing_in_vertex_list(constr_b, constr_b, constr_a);
  substitute_enclosing_in_vertex_list(constr_a, nullptr, constr_b);

  constr_a_vl->swap(*constr_b_vl);
  std::swap(constr_a.vl_with_info_ptr->may_share_subconstraint_with_others,
            constr_b.vl_with_info_ptr->may_share_subconstraint_with_others);
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
remove_constraint(Constraint_id cid)
{
  constraints_set.erase(cid);

  for_context_lists_of_all_subconstraints(cid, [&](Context_list* hcl, Vertex_it, Sc_it scit) {
    remove_first_in_context_list(hcl, cid);

    // If the constraint passes several times through the same subconstraint,
    // the above call to `remove_first_in_context_list` may remove them in the wrong order.

    // If this was the only context in the list, delete the context list
    if(hcl->empty()) {
      sc_to_c_map.erase(scit);
      delete hcl;
    }
  });
  cid.destroy();
}


// This function removes vertex v from the polyline constraint
// It only works for one polyline passing through v
// and for the case that the constrained edge u,w has no intersections
template <class T, class Compare, class Point>
void Polyline_constraint_hierarchy_2<T,Compare,Point>::simplify(Vertex_it uc,
                                                                Vertex_it vc,
                                                                Vertex_it wc)
{
  // TODO: How do we (want to) deal with u == w ???
  Vertex_handle u = *uc, v = *vc, w = *wc;
  typename Sc_to_c_map::iterator uv_sc_iter = sc_to_c_map.find(sorted_pair(u, v));
  typename Sc_to_c_map::iterator vw_sc_iter = sc_to_c_map.find(sorted_pair(v, w));
  Context_list*  uv_hcl = uv_sc_iter->second;
  Context_list*  vw_hcl = vw_sc_iter->second;
  // AF:  what is input() about???
  if(vc.input()){
    uc.input() = true;
    wc.input() = true;
  }

  // Take contexts from the two context lists depending on the orientation of the constraints
  // These are the contexts where current is either u or w
  // remove from uv_hcl the contexts where current is not u
  // remove from vw_hcl the contexts where current is not w
  // splice into uv_hcl
  typename Context_list::iterator it = uv_hcl->begin();
  while(it != uv_hcl->end()){
    if((*it->current()) != u){
      it = uv_hcl->erase(it);
    }else{
      // Remove the list item which points to v
      Vertex_list_ptr vertex_list = it->id().vl_ptr();
      Vertex_it vc_in_context = it->current();
      vc_in_context = std::next(vc_in_context);
      vertex_list->skip(vc_in_context.base());
      ++it;
    }
  }
  it = vw_hcl->begin();
  while(it != vw_hcl->end()){
    if((*it->current()) != w){
      it = vw_hcl->erase(it);
    }else{
      // Remove the list item which points to v
      Vertex_list_ptr vertex_list = it->id().vl_ptr();
      Vertex_it vc_in_context = it->current();
      vc_in_context = std::next(vc_in_context);
      vertex_list->skip(vc_in_context.base());
      ++it;
    }
  }

  uv_hcl->splice(uv_hcl->end(),*vw_hcl);
  delete vw_hcl;

  sc_to_c_map.erase(uv_sc_iter);
  sc_to_c_map.erase(vw_sc_iter);

  // reuse other context list
  sc_to_c_map[sorted_pair(u,w)] = uv_hcl;
}


template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::size_type
Polyline_constraint_hierarchy_2<T,Compare,Point>::remove_points_without_corresponding_vertex(Constraint_id cid)
{
  size_type n = 0;
  for(Point_it it = points_in_constraint_begin(cid);
      it != points_in_constraint_end(cid);) {
    if(cid.vl_ptr()->is_skipped(it.base())) {
      it = cid.vl_ptr()->erase(it.base());
      ++n;
    }else{
      ++it;
    }
  }
  return n;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::size_type
Polyline_constraint_hierarchy_2<T,Compare,Point>::remove_points_without_corresponding_vertex()
{
  size_type n = 0;
  for(const auto& cid : constraints_set){
    n+= remove_points_without_corresponding_vertex(cid);
  }
  return n;
}


template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::concatenate(Constraint_id constr_a, Constraint_id constr_b)
{
  // constr_a is [A, ..., M]
  // constr_b is [M, ..., B]
  // we want:
  //   constr_a = [A, ..., M, ..., B]
  //   constr_b = []

  // std::cerr << std::format("concatenate({}, {}) ", constr_a.id, constr_b.id) << std::endl;
  Vertex_list_ptr constr_a_vl = constr_a.vl_ptr();
  Vertex_list_ptr constr_b_vl = constr_b.vl_ptr();
  constraints_set.erase(constr_b);

  for_context_lists_of_all_subconstraints(constr_b, [&](Context_list* hcl, Vertex_it, Sc_it) {
    replace_first_in_context_list(hcl, constr_b, constr_a);
  });

  // now we really concatenate the vertex lists
  // Note that all iterators pointing into constr_a remain valid.
  // This concerns user code, as well as  the data member "pos" of the Context class
  CGAL_assertion(constr_a_vl->back().vertex() == constr_b_vl->front().vertex());
  constr_a_vl->pop_back(); // because it is the same as constr_b_vl.front()
  constr_a_vl->splice(constr_a_vl->skip_end(), *constr_b_vl, constr_b_vl->skip_begin(), constr_b_vl->skip_end());

  for_context_lists_of_all_subconstraints(constr_a, [&](Context_list* hcl, Vertex_it it, Sc_it) {
    update_first_context_position(hcl, constr_a, it);
  });

  if(constr_b.vl_with_info_ptr->may_share_subconstraint_with_others) {
    constr_a.vl_with_info_ptr->may_share_subconstraint_with_others = true;
  }
  constr_b.destroy();
  return constr_a;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::concatenate2(Constraint_id constr_a, Constraint_id constr_b)
{
  // constr_a is [A, ..., M]
  // constr_b is [M, ..., B]
  // we want:
  //   constr_a =
  //   constr_b = [A, ..., M, ..., B]

  Vertex_list_ptr constr_a_vl = constr_a.vl_ptr();
  Vertex_list_ptr constr_b_vl = constr_b.vl_ptr();
  constraints_set.erase(constr_a); // DIFF

  for_context_lists_of_all_subconstraints(constr_a, [&](Context_list* hcl, Vertex_it, Sc_it) { // DIFF
    replace_first_in_context_list(hcl, constr_a, constr_b); // DIFF
  });

  // now we really concatenate the vertex lists
  // Note that all iterators pointing into constr_b remain valid.
  CGAL_assertion(constr_a_vl->back().vertex() == constr_b_vl->front().vertex());
  constr_a_vl->pop_back(); // because it is the same as constr_b_vl.front()
  constr_b_vl->splice(constr_b_vl->skip_begin(), *constr_a_vl, constr_a_vl->skip_begin(), constr_a_vl->skip_end()); // DIFF

  // Note that for VC8 with iterator debugging the iterators pointing into constr_a
  // are NOT valid      So we have to update them
  for_context_lists_of_all_subconstraints(constr_b /*DIFF*/, [&](Context_list* hcl, Vertex_it it, Sc_it) {
    update_first_context_position(hcl, constr_b, it); // DIFF
  });

  if(constr_a.vl_with_info_ptr->may_share_subconstraint_with_others) {
    constr_b.vl_with_info_ptr->may_share_subconstraint_with_others = true;
  }
  constr_a.destroy(); // DIFF
  return constr_b; // DIFF
}


  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint
template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::split(Constraint_id constr, Vertex_it vcit)
{
  // constr is [A, ..., B], vcit points to M in [A, ..., B]

  Constraint_id new_constr = new_constraint_id();
  constraints_set.insert(new_constr);

  Vertex_list_ptr constr_vl = constr.vl_ptr();
  Vertex_list_ptr new_vl = new_constr.vl_ptr();

  vcit.input() = true;
  auto middle_node = Node(*vcit, true);

  // Let's split, that way:
  //       constr = [A, ..., M]
  //   new_constr = [M, ..., B]

  // The splice does:
  //       constr = [A, ..., M[ (without M)
  //   new_constr = [M, ..., B]
  new_vl->splice(new_vl->skip_end(), *constr_vl, vcit.base(), constr_vl->skip_end());
  constr_vl->push_back(middle_node); // add M to the end of constr

  CGAL_assertion(vcit.base() == new_vl->skip_begin());
  CGAL_assertion(new_vl->front().input() == true);
  CGAL_assertion(constr_vl->back().input() == true);

  for_context_lists_of_all_subconstraints(new_constr, [&](Context_list* hcl, Vertex_it, Sc_it) {
    replace_first_in_context_list(hcl, constr, new_constr);
  });

  new_constr.vl_with_info_ptr->may_share_subconstraint_with_others =
      constr.vl_with_info_ptr->may_share_subconstraint_with_others;
  return new_constr;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::split2(Constraint_id constr, Vertex_it vcit)
{
  // constr is [A, ..., B], vcit points to M in [A, ..., B]

  Constraint_id new_constr = new_constraint_id();
  constraints_set.insert(new_constr);

  Vertex_list_ptr constr_vl = constr.vl_ptr();
  Vertex_list_ptr new_vl = new_constr.vl_ptr();

  vcit.input() = true;
  auto middle_node = Node(*vcit, true);

  // Let's split, that way:
  //       constr = [M, ..., B]
  //   new_constr = [A, ..., M]

  // The splice does:
  //       constr = [M, ..., B]
  //   new_constr = [A, ..., M[  (without M)
  new_vl->splice(new_vl->skip_end(), *constr_vl, constr_vl->skip_begin(), vcit.base());
  new_vl->push_back(middle_node); // add M to the end of new_constr

  CGAL_assertion(vcit.base() == constr_vl->skip_begin());
  CGAL_assertion(constr_vl->front().input() == true);
  CGAL_assertion(new_vl->back().input() == true);

  for_context_lists_of_all_subconstraints(new_constr, [&](Context_list* hcl, Vertex_it, Sc_it) {
    replace_first_in_context_list(hcl, constr, new_constr);
  });

  new_constr.vl_with_info_ptr->may_share_subconstraint_with_others =
      constr.vl_with_info_ptr->may_share_subconstraint_with_others;
  return new_constr;
}


/*
when a constraint is inserted,
it is, at first, both  a constraint and a subconstraint
 */
template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::
insert_constraint(T va, T vb){
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  using CGAL::IO::oformat;
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "C_hierachy.insert_constraint( "
            << IO::oformat(va) << ", " << IO::oformat(vb) << ")\n";
#endif // CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  Subconstraint sc = sorted_pair(va, vb);
  Constraint_id cid = new_constraint_id();
  auto children = cid.vl_ptr();
  auto& fathers = sc_to_c_map[sc];
  if(fathers == nullptr){
    fathers = new Context_list;
  }

  children->push_front(Node(va, true));  // was sc.first
  children->push_back(Node(vb, true));   // was sc.second
  constraints_set.insert(cid);
  Context ctxt;
  ctxt.enclosing = cid;
  ctxt.pos       = children->skip_begin();
  fathers->push_front(ctxt);
  fix_contexts(*fathers);

  return cid;
}


template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::
insert_constraint_old_API(T va, T vb){
  return insert_constraint(va, vb);
}


template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
append_constraint(Constraint_id cid, T va, T vb){
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  using CGAL::IO::oformat;
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "C_hierachy.append_constraint( ..., "
            << IO::oformat(va) << ", " << IO::oformat(vb) << ")\n";
#endif // CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  Subconstraint sc = sorted_pair(va, vb);
  auto& fathers = sc_to_c_map[sc];
  if(fathers == nullptr){
    fathers = new Context_list;
  }

  typename Vertex_list::skip_iterator last_pos = std::prev(cid.end());
  CGAL_assertion(last_pos->vertex() == va);
  cid.vl_ptr()->push_back(Node(vb, true));
  Context ctxt;
  ctxt.enclosing = cid;
  ctxt.pos       = last_pos;
  fathers->push_front(ctxt);
  fix_contexts(*fathers);
}


template <class T, class Compare, class Point>
void Polyline_constraint_hierarchy_2<T,Compare,Point>::
clear()
{
  // clean and delete vertices lists
  for(auto cid : constraints()) {
    cid.clear();
    cid.destroy();
  }
  // clean and delete context lists
  for(auto& [_, cl_ptr] : sc_to_c_map) {
    cl_ptr->clear();
    delete cl_ptr;
  }
  sc_to_c_map.clear();
  constraints_set.clear();
}


template <class T, class Compare, class Point>
bool Polyline_constraint_hierarchy_2<T,Compare,Point>::
next_along_sc(T va, T vb, T& w) const
{
  // find the next vertex after vb along any enclosing constrained
  // return false if there is no ....
  Context_iterator  ctxtit, past;
  if(!get_contexts(va, vb, ctxtit, past)) CGAL_assertion(false);

  for(; ctxtit != past; ctxtit++) {
    Vertex_it pos = ctxtit->pos;
    if((*pos) == va) {
      ++pos;
      ++pos;
      if(pos != ctxtit->enclosing.end()) {
        w = *pos;
        return true;
      }
    } else {
      if(pos != ctxtit->enclosing.begin()) {
        --pos;
        w = *pos;
        return true;
      }
    }
  }
  return false;
}


/*
  same as add_Steiner
  precondition : va,vb est une souscontrainte.
*/
template <class T, class Compare, class Point>
void Polyline_constraint_hierarchy_2<T,Compare,Point>::
split_constraint(T va, T vb, T vc){
  add_Steiner(va, vb, vc);
}


template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
add_Steiner(const T va, const T vb, const T vc){
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  using CGAL::IO::oformat;
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "C_hierachy.add_Steinter( "
            << IO::oformat(va) << ", " << IO::oformat(vb) << ", " << IO::oformat(vc)
            << ")\n";
#endif // CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
  Sc_iterator sc_iter_va_vb = sc_to_c_map.find(sorted_pair(va, vb));
  if(sc_iter_va_vb == sc_to_c_map.end()) {
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
      std::cerr << CGAL::internal::cdt_2_indent_level
                << "  -> the constraint is already split\n";
#endif // CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
    return;
  }

  Context_list* va_vb_cl = sc_iter_va_vb->second;
  sc_to_c_map.erase(sc_iter_va_vb);
  Context_list* va_vc_cl = get_context_list(va,vc);
  Context_list* vc_vb_cl = get_context_list(vc,vb);

  bool vc_vb_was_already_a_subconstraint = vc_vb_cl != nullptr;

  if(vc_vb_cl == nullptr) {
    vc_vb_cl = new Context_list;
  }

  for(Context& ctxt : *va_vb_cl) {
    Vertex_it pos = ctxt.current();
    Vertex_it next_pos = std::next(pos);

    // [pos, next_pos] is the subconstraint (va,vb) or (vb,va)
    CGAL_assertion( (va == *pos && vb == *next_pos) || (vb == *pos && va == *next_pos) );
    // insert vc in enclosing constraint, just before `next_pos`
    Vertex_it pos_vc = ctxt.enclosing.vl_ptr()->insert(next_pos.base(), Node(vc));
    pos = std::prev(pos_vc);

    // now (pos, pos_vc, next_pos) is (va,vc,vb) or (vb,vc,va)

    // change ctxt in va_vb_cl to the context of (va,vc)
    Context& va_vc_ctxt = ctxt;
    Context vc_vb_ctxt{ctxt};
    if(*pos == va) { // (pos, pos_vc, next_pos) is (va,vc,vb)
      va_vc_ctxt.pos = pos;
      vc_vb_ctxt.pos = pos_vc;
    }
    else { // (pos, pos_vc, next_pos) is (vb,vc,va)
      vc_vb_ctxt.pos = pos;
      va_vc_ctxt.pos = pos_vc;
    }
    vc_vb_cl->push_back(vc_vb_ctxt);
  }

  if (va_vc_cl != nullptr) { // (va,vc) was already a subconstraint
    va_vc_cl->splice(va_vc_cl->end(), *va_vb_cl);
    delete va_vb_cl;
  } else {
    va_vc_cl = va_vb_cl;
    sc_to_c_map.emplace(sorted_pair(va,vc), va_vc_cl);
  }

  if (false == vc_vb_was_already_a_subconstraint) {
    sc_to_c_map.emplace(sorted_pair(vc,vb), vc_vb_cl);
  }
  fix_contexts(*va_vc_cl);
  fix_contexts(*vc_vb_cl);
}


template <class T, class Compare, class Point>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Subconstraint
Polyline_constraint_hierarchy_2<T,Compare,Point>::
sorted_pair(T va, T vb) const
{
  return comp(va, vb) ? Subconstraint(va,vb) : Subconstraint(vb,va);
}

template <class T, class Compare, class Point>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Context_list*
Polyline_constraint_hierarchy_2<T,Compare,Point>::
get_context_list(T va, T vb) const
{
  Sc_iterator sc_iter = sc_to_c_map.find(sorted_pair(va, vb));
  if(sc_iter == sc_to_c_map.end())
    return nullptr;
  else
    return sc_iter->second;
}

template <class T, class Compare, class Point>
inline
bool
Polyline_constraint_hierarchy_2<T,Compare,Point>::
get_contexts(T va, T vb,
             Context_iterator& ctxt,
             Context_iterator& past) const
{
  Context_list* hcl = get_context_list(va,vb);
  if (nullptr == hcl) return false;
  ctxt = hcl->begin();
  past = hcl->end();
  return true;
}



template <class T, class Compare, class Point>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Vertex_it
Polyline_constraint_hierarchy_2<T,Compare,Point>::
get_pos(T va, T vb) const
  //return pos in the first context
{
    return (*sc_to_c_map.find(sorted_pair(va,vb))).second->begin().pos;
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
oriented_end(T va, T vb, T& vc) const
{
  Context_iterator ctxt, past;
  if(!get_contexts(va,vb, ctxt, past) ) CGAL_assertion(false);
  if(*(ctxt->pos) == va)
    vc = ctxt->enclosing.back();
  else
    vc = ctxt->enclosing.front();
}


template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
print(std::ostream& os) const
{
  std::map<Vertex_handle,int>  vertex_num_map;
  int num = 0;
  for(const auto& cid : constraints()) {
    for (const auto& node : cid.elements()){
      vertex_num_map.emplace(node.vertex(), ++num);
    }
  }

  auto disp_vertex = [&vertex_num_map](Vertex_handle v) {
    return CGAL::IO::oformat(
        [v, &vertex_num_map](auto& os) -> decltype(os)& {
          constexpr bool star_v_has_timestamp =
              internal::has_timestamp_v<CGAL::cpp20::remove_cvref_t<decltype(*v)>>;
          if constexpr(star_v_has_timestamp) {
            return os << '#' << v->time_stamp();
          } else {
            auto it = vertex_num_map.find(v);
            auto n = (it == vertex_num_map.end()) ? -1 : it->second;
            return os << n;
          }
        },
        IO_manip_tag{});
  };

  os << "# number of constraints: " << number_of_constraints() << std::endl;
  for(const auto& cid : constraints()) {
    os << "constraint(" << cid.id << ") ";
    os << cid.vl_ptr();
    os << "\n  vertex list ";
    for(const auto& node : cid.elements()) {
      os << disp_vertex(node.vertex()) << " ";
    }
    os << "\n  corresponding points:  ";
    for(const auto& node : cid.elements()) {
      os << node.point() << "   ";
    }
    if(cid.vl_with_info_ptr->may_share_subconstraint_with_others) {
      os << "\n  (may have non-simple context lists)";
    }
    os << std::endl;
  }
  os << std::endl;
  os << "# number of subconstraints: " << sc_to_c_map.size() << std::endl;
  for(const auto& subconstraint : subconstraints()) {
    os << "subconstraint (";
    os << disp_vertex(subconstraint.first) << ", "
              << disp_vertex(subconstraint.second) << ")";
    Context_iterator cb, ce;
    get_contexts(subconstraint.first, subconstraint.second, cb, ce);

    os << "  enclosing:  ";
    for(; cb != ce; cb++) {
      os << "(cid " << cb->id().id << ") " << cb->id().vl_ptr();
      os << ", pos: " << std::distance(cb->vertices_begin(), cb->pos) << "   ";
    }
    os << std::endl;
  }
  return;
}

} //namespace CGAL
#endif // CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
