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
#include <array>
#include <queue>
#include <iterator>

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
  using T_point_ref = decltype(std::declval<T>()->point());
  static_assert(std::is_same_v<CGAL::cpp20::remove_cvref_t<Point>, CGAL::cpp20::remove_cvref_t<T_point_ref>>,
                "The point type of the vertex handle must be the same as the point type of the hierarchy.");
public:
  using Vertex_handle = T;
  using Vertex_handle_compare = Compare;
  using Subconstraint = std::pair<T, T>;

  using size_type = typename Vertex_handle::size_type;

private:
  class Node {
  public:
    explicit Node(Vertex_handle vh, bool input = false)
      : vertex_(vh), point_(vertex_->point()), input_(input)
    {}
    const Point& point() const { return point_; }
    const Vertex_handle& vertex() const { return vertex_; }
    bool& input() { return input_; }
    const bool& input() const { return input_; }
  private:
    Vertex_handle vertex_;
    Point point_;
    bool input_;
  };

  using Vertex_list = CGAL::Skiplist<Node>;
  using Vertex_list_ptr = Vertex_list*;

public:
  // the base line is always
  class Point_it : public boost::stl_interfaces::iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
      Point_it,
#endif
      std::bidirectional_iterator_tag,
      const Point>
  {
public:
    using Self = Point_it;
    using Base = boost::stl_interfaces::iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
        Point_it,
#endif
        std::bidirectional_iterator_tag,
        const Point>;
    using Base_it = typename Vertex_list::all_iterator;

    Base_it base() const { return it; }
    Base_it& base_reference() { return it; }
    const Base_it& base_reference() const { return it; }

    Point_it() = default;
    Point_it(Base_it it) : it(it) {}

    const Point& operator*() const { return base()->point(); }

    Self& operator++() { ++it; return *this; }
    Self operator++(int i) { return static_cast<Self>(Base::operator++(i)); }

    Self& operator--() { --it; return *this; }
    Self operator--(int i) { return static_cast<Self>(Base::operator--(i)); }

  private:
    friend bool operator==(const Self& lhs, const Self& rhs) {
      return lhs.base() == rhs.base();
    }
    friend bool operator==(const Self& lhs, const Base_it& rhs) {
      return lhs.base() == rhs;
    }
    friend bool operator==(const Base_it& lhs, const Self& rhs) {
      return lhs == rhs.base();
    }
    Base_it it;
  };

  BOOST_STL_INTERFACES_STATIC_ASSERT_CONCEPT(Point_it, std::bidirectional_iterator);
#if BOOST_VERSION >= 108300
  BOOST_STL_INTERFACES_STATIC_ASSERT_ITERATOR_TRAITS(Point_it, std::bidirectional_iterator_tag,
    std::bidirectional_iterator, Point, const Point&, const Point*, std::ptrdiff_t);
#endif

  // only nodes with a vertex_handle that is still in the triangulation
  class Vertex_it : public boost::stl_interfaces::iterator_interface<
  #if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
        Vertex_it,
  #endif
        std::bidirectional_iterator_tag,
        const Vertex_handle>
  {
  public:
    using Self = Vertex_it;
    using Base = boost::stl_interfaces::iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
          Vertex_it,
#endif
          std::bidirectional_iterator_tag,
          const Vertex_handle>;
    using Base_it = typename Vertex_list::skip_iterator;

    Base_it base() const { return it; }

    Base_it& base_reference() { return it; }
    const Base_it& base_reference() const { return it; }

    Vertex_it() = default;
    Vertex_it(Base_it it) : it(it) {
      [[maybe_unused]] Vertex_it itt;
      static_assert(std::is_same_v<decltype(++itt), Vertex_it&>);
    }

    const Vertex_handle& operator*() const { return this->base()->vertex(); }

    Self& operator++() { ++it; return *this; }
    Self operator++(int i) { return static_cast<Self>(Base::operator++(i)); }

    Self& operator--() { --it; return *this; }
    Self operator--(int i) { return static_cast<Self>(Base::operator--(i)); }

    operator Point_it() const { return Point_it(this->base()); }

    bool& input() { return this->base()->input(); }
    const bool& input() const { return this->base()->input(); }
  private:
    friend bool operator==(const Self& lhs, const Self& rhs) {
      return lhs.base() == rhs.base();
    }
    friend bool operator==(const Self& lhs, const Base_it& rhs) {
      return lhs.base() == rhs;
    }
    friend bool operator==(const Base_it& lhs, const Self& rhs) {
      return lhs == rhs.base();
    }
    Base_it it;
  };
  BOOST_STL_INTERFACES_STATIC_ASSERT_CONCEPT(Vertex_it, std::bidirectional_iterator);
#if BOOST_VERSION >= 108300
  BOOST_STL_INTERFACES_STATIC_ASSERT_ITERATOR_TRAITS(Vertex_it, std::bidirectional_iterator_tag,
    std::bidirectional_iterator, Vertex_handle, const Vertex_handle&, const Vertex_handle*, std::ptrdiff_t);
#endif
  struct Vertex_list_with_info {
    const Polyline_constraint_hierarchy_2* hierarchy_ptr = nullptr;
    Vertex_list vl{};
    bool may_share_subconstraint_with_others = false;
  };

  class Constraint_id
  {
    Vertex_list_with_info* vl_with_info_ptr = nullptr;
    size_type id = (std::numeric_limits<size_type>::max)();
  public:
    Constraint_id(std::nullptr_t = nullptr) {}
    Constraint_id(Vertex_list_with_info* ptr, size_type id) : vl_with_info_ptr(ptr), id(id) {}

    void destroy() {
      if(vl_with_info_ptr != nullptr) {
        delete vl_with_info_ptr;
        (*this) = nullptr;
      }
    }

    auto index() const { return id; }

    Vertex_list_with_info* vl_with_info_pointer() const { return vl_with_info_ptr; }

    Vertex_list_ptr vl_ptr() const {
      return vl_with_info_ptr == nullptr ? nullptr : std::addressof(vl_with_info_ptr->vl);
    }

    bool  may_share() const { return vl_with_info_ptr->may_share_subconstraint_with_others; }
    bool& may_share()       { return vl_with_info_ptr->may_share_subconstraint_with_others; }

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
      *this = Constraint_id{};
      return *this;
    }
    bool operator==(std::nullptr_t n) const { return vl_with_info_ptr == n; }
    bool operator!=(std::nullptr_t n) const { return vl_with_info_ptr != n; }

    bool operator==(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (index() == other.index()));
      return vl_ptr() == other.vl_ptr();
    }

    bool operator!=(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (index() == other.index()));
      return vl_ptr() != other.vl_ptr();
    }

    bool operator<(const Constraint_id& other) const
    {
      CGAL_assertion(same_hierarchy(other));
      CGAL_assertion((vl_ptr() == other.vl_ptr()) == (index() == other.index()));
      return index() < other.index();
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
    Constraint_id enclosing = nullptr;
    Vertex_it     pos{};
  public:
    Context() = default;
    Context(Constraint_id enclosing, Vertex_it pos) : enclosing(enclosing), pos(pos) {}
    Context(Constraint_id enclosing, typename Vertex_list::skip_iterator pos) : enclosing(enclosing), pos(pos) {}

    Vertex_it    vertices_begin()const { return enclosing.begin();}
    Vertex_it    current()const {return pos;}
    Vertex_it    vertices_end()const {return enclosing.end();}
    Constraint_id  id()const { return enclosing; }
    size_type    number_of_vertices() const {return enclosing.size(); }
  };

  using Context_list = std::list<Context>;
  using Context_iterator = typename Context_list::iterator;

  static void fix_may_share_in_contexts_constraints(Context_list& context_list) {
    const auto cl_size = context_list.size();
    switch(cl_size) {
    case 0:
      CGAL_unreachable();
    case 1:
      return;
    default:
      for(auto& context : context_list) {
        CGAL_assertion(context.enclosing != nullptr);
        context.enclosing.may_share() = true;
      }
    }
  }

  using Constraints_set = std::set<Constraint_id>;
#if CGAL_USE_BARE_STD_MAP
  using Sc_to_c_map = std::map<Subconstraint, Context_list*, Pair_compare>;
#else
  using Sc_to_c_map = CGAL::unordered_flat_map<Subconstraint, Context_list*, boost::hash<Subconstraint>>;
#endif
  using Constraint_iterator = typename Constraints_set::iterator;
  using Constraints = const Constraints_set&;
  using Sc_iterator = typename Sc_to_c_map::const_iterator;
  using Sc_it = typename Sc_to_c_map::iterator;
  using Subconstraint_and_contexts_iterator = Sc_iterator;
  using Subconstraints_and_contexts = const Sc_to_c_map&;

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
      return hierarchy == nullptr;
    }

    bool is_valid() const {
      return (hierarchy != nullptr);
    }

    bool is_end() const {
      return constraint_it == hierarchy->constraints_end();
    }

    bool is_begin() const {
      return constraint_it == hierarchy->constraints_begin() &&
             vertex_it == begin_or_null(constraint_it);
    }

    bool is_dereferenceable() const {
      return is_valid() &&
          constraint_it != hierarchy->constraints_end() &&
          vertex_it != constraint_it->end() &&
          std::next(vertex_it) != constraint_it->end();
    }

    bool equal(const Subconstraint_iterator& other) const {
      if(hierarchy != other.hierarchy) return false;
      if(is_singular()) return true;

      return (constraint_it == other.constraint_it) && (this->is_end() == other.is_end());
    }

    Vertex_it begin_or_null(Constraint_iterator constraint_it) const {
      if(constraint_it == hierarchy->constraints_end()) {
        return Vertex_it();
      }
      return constraint_it->begin();
    }

    bool already_seen() const {
      if(is_end()) {
        return false;
      }
      if(false == constraint_it->may_share()) {
        return false;
      }
      auto [va, vb] = this->operator*();
      auto it = hierarchy->find_contexts(va, vb);
      CGAL_assertion(it != hierarchy->contexts_not_found());

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
        , constraint_it(hierarchy->constraints_begin())
        , vertex_it(begin_or_null(hierarchy->constraints_begin()))
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
        , constraint_it(hierarchy->constraints_end())
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
        if(constraint_it == hierarchy->constraints_end() || vertex_it == constraint_it->begin()) {
          --constraint_it;
          vertex_it = std::prev(constraint_it->end(), 2);
        } else {
          --vertex_it;
        }
      } while(already_seen());
      return *this;
    }
  }; // end class Subconstraint_iterator
  using Subconstraints = Iterator_range<Subconstraint_iterator>;

private:
  struct Priv { // encapsulate the private members in a struct, to detect direct access to them
    Priv(Compare comp)
    : comp(comp)
    #if CGAL_USE_BARE_STD_MAP
    , sc_to_c_map(Pair_compare(comp))
#else
    , sc_to_c_map()
#endif
    {}

    Compare               comp;
    Sc_to_c_map           sc_to_c_map;
    std::queue<size_type> free_ids;
    Constraints_set       constraints_set;
  } priv;
public:
  Polyline_constraint_hierarchy_2(const Compare& comp) : priv(comp)  {}
  Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch) : priv(ch.priv.comp) {}
  Polyline_constraint_hierarchy_2(Polyline_constraint_hierarchy_2&&) = default;

  ~Polyline_constraint_hierarchy_2(){ clear();}
  void clear();

  Polyline_constraint_hierarchy_2& operator=(const Polyline_constraint_hierarchy_2& ch) { return copy(ch); }
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

  std::array<Vertex_handle, 2> enclosing_constraint(T vaa, T vbb) const;

  Context context(T va, T vb);
  size_type number_of_enclosing_constraints(T va, T vb) const;
  Context_iterator contexts_begin(T va, T vb) const;
  Context_iterator contexts_end(T va, T vb) const;
  using Contexts = Iterator_range<Context_iterator>;
  Contexts contexts(T va, T vb) const;
  Context_list* get_context_list(T va, T vb) const;

  size_type number_of_constraints() const { return priv.constraints_set.size(); }
  size_type number_of_subconstraints() const { return priv.sc_to_c_map.size(); }

  // insert/remove
  void add_Steiner(const T va, const T vb, const T vx);
  Constraint_id insert_constraint_old_API(T va, T vb);
  Constraint_id insert_constraint(T va, T vb);
  void append_constraint(Constraint_id cid, T va, T vb);
  void swap(Constraint_id first, Constraint_id second);
  void remove_constraint(Constraint_id cid);

  void split_constraint(T va, T vb, T vc);

  void simplify(Vertex_it u, Vertex_it v, Vertex_it w);

  size_type remove_points_without_corresponding_vertex(Constraint_id);
  size_type remove_points_without_corresponding_vertex();

  Constraint_id concatenate(Constraint_id first, Constraint_id&& second);
  Constraint_id prepend(Constraint_id&& first, Constraint_id second);
  Constraint_id split_tail(Constraint_id first, Vertex_it vcit);
  Constraint_id split_head(Constraint_id first, Vertex_it vcit);

  // iterators

  Subconstraint_and_contexts_iterator subconstraints_and_contexts_begin() const
  {
    return priv.sc_to_c_map.begin();
  }

  Subconstraint_and_contexts_iterator subconstraints_and_contexts_end() const
  {
    return priv.sc_to_c_map.end();
  }

  Subconstraint_iterator subconstraints_begin() const {
    BOOST_STL_INTERFACES_STATIC_ASSERT_CONCEPT(Subconstraint_iterator, std::bidirectional_iterator);
#if BOOST_VERSION >= 108300
    BOOST_STL_INTERFACES_STATIC_ASSERT_ITERATOR_TRAITS(
      Subconstraint_iterator, std::bidirectional_iterator_tag, std::bidirectional_iterator,
      Subconstraint, Subconstraint, typename Subconstraint_iterator::pointer, std::ptrdiff_t);
#endif
    return Subconstraint_iterator(Subconstraint_iterator::Construction_access::begin_tag(),
                                  this);
  }

  Subconstraint_iterator subconstraints_end() const {
    return Subconstraint_iterator(Subconstraint_iterator::Construction_access::end_tag(),
                                  this);
  }

  Constraint_iterator  constraints_begin()  const{ return priv.constraints_set.begin(); }
  Constraint_iterator  constraints_end()    const{ return priv.constraints_set.end();   }

  // Ranges
  const auto& constraints() const { return priv.constraints_set; }
  const auto& subconstraints_and_contexts() const { return priv.sc_to_c_map; }
  auto subconstraints() const {
    return Iterator_range<Subconstraint_iterator>(subconstraints_begin(), subconstraints_end());
  }

  // Helper functions
  Polyline_constraint_hierarchy_2& copy(const Polyline_constraint_hierarchy_2& ch);
  Polyline_constraint_hierarchy_2& copy(const Polyline_constraint_hierarchy_2& ch,
                                        std::map<Vertex_handle,Vertex_handle>& vmap);
  void swap(Polyline_constraint_hierarchy_2& ch);

private:
  // a few member functions to encapsulate more of the uses of `sc_to_c_map`
  auto find_contexts(Vertex_handle va, Vertex_handle vb) { return priv.sc_to_c_map.find(sorted_pair(va, vb)); }
  auto find_contexts(Vertex_handle va, Vertex_handle vb) const { return priv.sc_to_c_map.find(sorted_pair(va, vb)); }
  auto contexts_not_found() { return priv.sc_to_c_map.end(); }
  auto contexts_not_found() const { return priv.sc_to_c_map.end(); }
  void erase_context(Sc_iterator it) { priv.sc_to_c_map.erase(it); }
  auto& contexts_of(Vertex_handle va, Vertex_handle vb) { return priv.sc_to_c_map[sorted_pair(va, vb)]; }
  //
  // then the uses of `constraints_set`
  Constraint_id create_new_constraint() {
    size_type id; // uninitialized
    if(priv.free_ids.empty()) {
      id = priv.constraints_set.size();
    } else {
      id = priv.free_ids.front();
      priv.free_ids.pop();
    }
    Constraint_id cid{new Vertex_list_with_info{this}, id};
    priv.constraints_set.insert(cid);
    return cid;
  }

  void erase_constraint(Constraint_id cid) {
    priv.free_ids.push(cid.index());
    priv.constraints_set.erase(cid);
    cid.destroy();
  }

  //
  // functions to traverse and act on the context lists
  //
  template <typename F>
  void for_context_lists_of_all_subconstraints(Constraint_id cid, const F& f)
  {
    auto vl = cid.vl_ptr();
    for(Vertex_it it = vl->skip_begin(), succ = it, end = vl->skip_end(); ++succ != end; ++it) {
      auto scit = find_contexts(*it, *succ);
      CGAL_assertion(scit != contexts_not_found());
      Context_list* hcl = scit->second;
      f(hcl, it, scit);
    }
  }
  //
  static void replace_first_in_context_list(Context_list* hcl, Constraint_id old_id, Constraint_id new_id)
  {
    // std::find_if is a sort of std::for_each with a break when the lambda returns true
    [[maybe_unused]] auto it = std::find_if(hcl->begin(), hcl->end(), [&](Context& ctxt) {
      if(ctxt.enclosing == old_id) {
        ctxt.enclosing = new_id;
        return true;
      }
      return false;
    });
  }
  //
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
  //
  static void remove_first_in_context_list(Context_list* hcl, Constraint_id id)
  {
    auto it = std::find_if(hcl->begin(), hcl->end(), [&](Context& ctxt) {
      return ctxt.enclosing == id;
    });
    if(it != hcl->end()) {
      hcl->erase(it);
    }
  }

  //
  // other utilities as private member functions
  //
  Subconstraint sorted_pair(T va, T vb) const {
    return priv.comp(va, vb) ? Subconstraint(va,vb) : Subconstraint(vb,va);
  }
  Subconstraint sorted_pair(Subconstraint sc) {
    const auto& [va, vb] = sc; return sorted_pair(va, vb);
  }

public:
  //to_debug
  void print(std::ostream& os = std::cout) const;
};

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
copy(const Polyline_constraint_hierarchy_2& other) -> Polyline_constraint_hierarchy_2&
{
  // create a identity transfer vertex map
  std::map<Vertex_handle, Vertex_handle>  vmap;
  for(const auto& cid : other.constraints()) {
    for(const auto& node : cid.elements()) {
      auto v = node.vertex();
      vmap[v] = v;
    }
  }
  return copy(other, vmap);
}

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T, Compare, Point>::
copy(const Polyline_constraint_hierarchy_2& other, std::map<Vertex_handle, Vertex_handle>& vmap)
    -> Polyline_constraint_hierarchy_2&
// copy with a transfer vertex map
{
  std::map<Constraint_id, Constraint_id> cstr_map;
  clear();
  // copy constraints_set
  for(const auto& cid1: other.constraints()) {
    Constraint_id cid2 = create_new_constraint();
    cstr_map[cid1] = cid2;
    for(const auto& node : cid1.elements()) {
      cid2.vl_ptr()->push_back(Node(vmap[node.vertex()], node.input()));
    }
    cid2.may_share() = cid1.may_share();
  }
  // copy sc_to_c_map
  for(const auto& [sc1, hcl1] : other.subconstraints_and_contexts()) {
    Vertex_handle uu2 = vmap[sc1.first];
    Vertex_handle vv2 = vmap[sc1.second];
    Context_list* hcl2 = new Context_list;
    contexts_of(uu2, vv2) = hcl2;
    for(const auto& ctxt : *hcl1) {
      const auto cid1 = ctxt.id();
      const auto pos1 = ctxt.current();
      // vertices of the enclosing constraints
      Constraint_id cid2 = cstr_map[cid1];
      auto pos2 = std::next(Vertex_it(cid2.begin()),
                            std::distance(Vertex_it(cid1.begin()), pos1));
      hcl2->emplace_back(cid2, pos2);
    }
  }

  priv.comp = other.priv.comp;
  return *this;
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
swap(Polyline_constraint_hierarchy_2& ch)
{
  using std::swap;
  swap(priv.comp, ch.priv.comp);
  priv.free_ids.swap(ch.priv.free_ids);
  priv.constraints_set.swap(ch.priv.constraints_set);
  priv.sc_to_c_map.swap(ch.priv.sc_to_c_map);
}


template <class T, class Compare, class Point>
bool Polyline_constraint_hierarchy_2<T,Compare,Point>::
is_subconstraint(T va, T vb) const
{
  return( find_contexts(va, vb) != contexts_not_found() );
}


// used by Constrained_triangulation_plus_2::intersect with Exact_intersection_tag
template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
enclosing_constraint(T vaa, T vbb) const -> std::array<Vertex_handle, 2>
{
  std::array<Vertex_handle, 2> result;
  auto [hcit, past] = contexts(vaa, vbb);
  if (hcit == past) return result;
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
  result[0] = *pos;
  pos = std::next(hcit->pos);
  CGAL_assertion(vbb == *pos);
  while(!pos.input()){
    ++pos;
  }
  result[1] = *pos;
  return result;
}

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
context(T va, T vb) -> Context
{
  auto [hcit, past] = contexts(va, vb);
  CGAL_assertion(hcit != past); CGAL_USE(past);
  return *hcit;
}

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
number_of_enclosing_constraints(T va, T vb) const -> size_type
{
  Context_list* hcl = get_context_list(va,vb);
  CGAL_assertion(nullptr != hcl);
  return hcl->size();
}

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts_begin(T va, T vb) const -> Context_iterator
{
  return contexts(va, vb).begin();
}

template <class T, class Compare, class Point>
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts_end(T va, T vb) const -> Context_iterator
{
  return contexts(va, vb).end();
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

  substitute_enclosing_in_vertex_list(constr_a, constr_a, nullptr);
  substitute_enclosing_in_vertex_list(constr_b, constr_b, constr_a);
  substitute_enclosing_in_vertex_list(constr_a, nullptr, constr_b);

  constr_a.vl_ptr()->swap(*constr_b.vl_ptr());
  std::swap(constr_a.may_share(), constr_b.may_share());
}

template <class T, class Compare, class Point>
void
Polyline_constraint_hierarchy_2<T,Compare,Point>::
remove_constraint(Constraint_id cid)
{
  for_context_lists_of_all_subconstraints(cid, [&](Context_list* hcl, Vertex_it, Sc_it scit) {
    remove_first_in_context_list(hcl, cid);

    // If the constraint passes several times through the same subconstraint,
    // the above call to `remove_first_in_context_list` may remove them in the wrong order.

    // If this was the only context in the list, delete the context list
    if(hcl->empty()) {
      erase_context(scit);
      delete hcl;
    }
  });

  erase_constraint(cid);
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
  typename Sc_to_c_map::iterator uv_sc_iter = find_contexts(u, v);
  typename Sc_to_c_map::iterator vw_sc_iter = find_contexts(v, w);
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

  erase_context(uv_sc_iter);
  erase_context(vw_sc_iter);

  // reuse other context list
  contexts_of(u,w) = uv_hcl;
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
  for(const auto& cid : constraints()){
    n+= remove_points_without_corresponding_vertex(cid);
  }
  return n;
}


template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::concatenate(Constraint_id constr_a, Constraint_id&& constr_b)
{
  if(constr_a == nullptr) {
    swap(constr_a, constr_b);
  };
  if(constr_b == nullptr) return constr_a;

  // constr_a is [A, ..., M]
  // constr_b is [M, ..., B]
  // we want:
  //   constr_a = [A, ..., M, ..., B]
  //   constr_b = []

  // std::cerr << std::format("concatenate({}, {}) ", constr_a.id, constr_b.id) << std::endl;
  Vertex_list_ptr constr_a_vl = constr_a.vl_ptr();
  Vertex_list_ptr constr_b_vl = constr_b.vl_ptr();

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

  if(constr_b.may_share()) {
    constr_a.may_share() = true;
  }
  erase_constraint(constr_b);
  return constr_a;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::prepend(Constraint_id&& constr_a, Constraint_id constr_b)
{
  if(constr_b == nullptr) {
    swap(constr_a, constr_b);
  };
  if(constr_a == nullptr) return constr_b;

  // constr_a is [A, ..., M]
  // constr_b is [M, ..., B]
  // we want:
  //   constr_a =
  //   constr_b = [A, ..., M, ..., B]

  Vertex_list_ptr constr_a_vl = constr_a.vl_ptr();
  Vertex_list_ptr constr_b_vl = constr_b.vl_ptr();

  for_context_lists_of_all_subconstraints(constr_a, [&](Context_list* hcl, Vertex_it, Sc_it) { // DIFF
    replace_first_in_context_list(hcl, constr_a, constr_b); // DIFF
  });

  // now we really concatenate the vertex lists
  // Note that all iterators pointing into constr_b remain valid.
  CGAL_assertion(constr_a_vl->back().vertex() == constr_b_vl->front().vertex());
  constr_a_vl->pop_back(); // because it is the same as constr_b_vl.front()
  constr_b_vl->splice(constr_b_vl->skip_begin(), *constr_a_vl, constr_a_vl->skip_begin(), constr_a_vl->skip_end()); // DIFF

  for_context_lists_of_all_subconstraints(constr_b /*DIFF*/, [&](Context_list* hcl, Vertex_it it, Sc_it) {
    update_first_context_position(hcl, constr_b, it); // DIFF
  });

  if(constr_a.may_share()) {
    constr_b.may_share() = true;
  }
  erase_constraint(constr_a); // DIFF
  return constr_b; // DIFF
}


  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint
template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::split_tail(Constraint_id constr, Vertex_it vcit)
{
  // constr is [A, ..., B], vcit points to M in [A, ..., B]

  Constraint_id new_constr = create_new_constraint();

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

  bool new_constr_share = false;
  for_context_lists_of_all_subconstraints(new_constr, [&](Context_list* hcl, Vertex_it, Sc_it) {
    if(constr.may_share() && hcl->size() > 1) {
      new_constr_share = true;
    }
    replace_first_in_context_list(hcl, constr, new_constr);
  });

  new_constr.may_share() = new_constr_share;
  return new_constr;
}

template <class T, class Compare, class Point>
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Point>::split_head(Constraint_id constr, Vertex_it vcit)
{
  // constr is [A, ..., B], vcit points to M in [A, ..., B]

  Constraint_id new_constr = create_new_constraint();

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

  bool new_constr_share = false;
  for_context_lists_of_all_subconstraints(new_constr, [&](Context_list* hcl, Vertex_it, Sc_it) {
    if(constr.may_share() && hcl->size() > 1) {
      new_constr_share = true;
    }
    replace_first_in_context_list(hcl, constr, new_constr);
  });

  new_constr.may_share() = new_constr_share;
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
  Constraint_id cid = create_new_constraint();
  auto& context_list_ptr = contexts_of(va, vb);
  if(context_list_ptr == nullptr){
    context_list_ptr = new Context_list;
  }

  auto children = cid.vl_ptr();
  children->push_front(Node(va, true));
  children->push_back(Node(vb, true));
  context_list_ptr->emplace_front(cid, cid.begin());
  fix_may_share_in_contexts_constraints(*context_list_ptr);

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
  auto& context_list_ptr = contexts_of(va, vb);
  if(context_list_ptr == nullptr){
    context_list_ptr = new Context_list;
  }

  auto pos_va = std::prev(cid.end());
  CGAL_assertion(pos_va->vertex() == va);
  cid.vl_ptr()->push_back(Node(vb, true));
  context_list_ptr->emplace_front(cid, pos_va);
  fix_may_share_in_contexts_constraints(*context_list_ptr);
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
  for(auto& [_, cl_ptr] : priv.sc_to_c_map) {
    cl_ptr->clear();
    delete cl_ptr;
  }
  priv = Priv(priv.comp);
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
  Sc_iterator sc_iter_va_vb = find_contexts(va, vb);
  if(sc_iter_va_vb == contexts_not_found()) {
#ifdef CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
      std::cerr << CGAL::internal::cdt_2_indent_level
                << "  -> the constraint is already split\n";
#endif // CGAL_DEBUG_POLYLINE_CONSTRAINT_HIERARCHY_2
    return;
  }

  Context_list* va_vb_cl = sc_iter_va_vb->second;
  erase_context(sc_iter_va_vb);
  Context_list*& vc_vb_cl_ref = contexts_of(vc,vb);

  if(vc_vb_cl_ref == nullptr) {
    vc_vb_cl_ref = new Context_list;
  }

  Context_list* vc_vb_cl = vc_vb_cl_ref;

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

  Context_list*& va_vc_cl = contexts_of(va,vc);
  if (va_vc_cl != nullptr) { // (va,vc) was already a subconstraint
    va_vc_cl->splice(va_vc_cl->end(), *va_vb_cl);
    delete va_vb_cl;
  } else {
    va_vc_cl = va_vb_cl;
  }

  fix_may_share_in_contexts_constraints(*va_vc_cl);
  fix_may_share_in_contexts_constraints(*vc_vb_cl);
}


template <class T, class Compare, class Point>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Point>::Context_list*
Polyline_constraint_hierarchy_2<T,Compare,Point>::
get_context_list(T va, T vb) const
{
  Sc_iterator sc_iter = find_contexts(va, vb);
  if(sc_iter == contexts_not_found())
    return nullptr;
  else
    return sc_iter->second;
}

template <class T, class Compare, class Point>
inline
auto Polyline_constraint_hierarchy_2<T,Compare,Point>::
contexts(T va, T vb) const -> Iterator_range<Context_iterator>
{
  Context_list* hcl = get_context_list(va,vb);
  if (nullptr == hcl) return {};
  else return {hcl->begin(), hcl->end()};
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
            CGAL_USE(vertex_num_map);
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
    os << "constraint(" << cid.index() << ") ";
    os << cid.vl_ptr();
    os << "\n  vertex list ";
    for(const auto& node : cid.elements()) {
      os << disp_vertex(node.vertex()) << " ";
    }
    os << "\n  corresponding points:  ";
    for(const auto& node : cid.elements()) {
      os << node.point() << "   ";
    }
    if(cid.may_share()) {
      os << "\n  (may have non-simple context lists)";
    }
    os << std::endl;
  }
  os << std::endl;
  os << "# number of subconstraints: " << number_of_subconstraints() << std::endl;
  for(const auto& subconstraint : subconstraints()) {
    os << "subconstraint (";
    os << disp_vertex(subconstraint.first) << ", "
       << disp_vertex(subconstraint.second) << ")";

    os << "  enclosing:  ";
    for(const auto& ctxt : contexts(subconstraint.first, subconstraint.second)) {
      os << "(cid " << ctxt.id().index() << ") " << ctxt.id().vl_ptr();
      os << ", pos: " << std::distance(ctxt.vertices_begin(), ctxt.pos) << "   ";
    }
    os << std::endl;
  }
  return;
}

} //namespace CGAL
#endif // CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
