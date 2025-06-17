// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    Andreas Fabri, Laurent Rineau

#ifndef CGAL_TDS_INTERNAL_INDEXED_STORAGE_3_H
#define CGAL_TDS_INTERNAL_INDEXED_STORAGE_3_H

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/license/TDS_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/TDS_3/internal/Dummy_tds_3.h>
#include <CGAL/TDS_3/internal/Triangulation_ds_circulators_3.h>
#include <CGAL/TDS_3/internal/Triangulation_ds_iterators_3.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_data_structure_3_fwd.h>

#include <CGAL/assertions.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Time_stamper.h>
#include <CGAL/utility.h>

#include <boost/stl_interfaces/iterator_interface.hpp>

#include <array>
#include <iostream>
#include <list>
#include <type_traits>
#include <vector>
#include <optional>

#ifdef CGAL_LINKED_WITH_TBB
#  include <tbb/enumerable_thread_specific.h>
#endif

namespace CGAL {

//utilities for copy_tds
namespace internal { namespace TDS_3{
  template <class Vertex_src,class Vertex_tgt>
  struct Default_index_vertex_converter;

  template <class Cell_src,class Cell_tgt>
  struct Default_index_cell_converter;

  template <class Vertex>
  struct Default_index_vertex_converter<Vertex,Vertex>
  {
    const Vertex& operator()(const Vertex& src) const {
      return src;
    }

    void operator()(const Vertex&,Vertex) const {}
  };

  template <class Cell>
  struct Default_index_cell_converter<Cell,Cell>{
    const Cell& operator()(const Cell& src) const {
      return src;
    }

    void operator()(const Cell&,Cell) const {}
  };
} } //namespace internal::TDS_3

  template <typename TDS_3 = void>
  struct Cell {

    using Triangulation_data_structure = TDS_3;
    using Vertex_handle = typename TDS_3::Vertex_handle;
    using Cell_handle   = typename TDS_3::Cell_handle;
    using Vertex_index  = typename TDS_3::Vertex_index;
    using Cell_index    = typename TDS_3::Cell_index;
    using Index         = Cell_index;
    using TDS_data      = typename TDS_3::Cell_data;

    struct Storage {
      std::array<Vertex_index,4> ivertices{ {{}, {}, {}, {}} };
      std::array<Cell_index,4>   ineighbors{ {{}, {}, {}, {}} };
    };

    Cell()
      : tds_(nullptr), index_(Cell_handle().index())
    {}

    Cell(TDS_3* tds_, Cell_index index)
      : tds_(tds_), index_(index)
    {}

    Cell& operator=(const Cell& other)
    {
      if (this != &other) {
        tds_ = other.tds();
        index_ = other.index_;
      }
      return *this;
    }

    Cell(const Cell& other)
      : tds_(other.tds()), index_(other.index_)
    {}

    bool operator==(const Cell& other) const
    {
      return tds_ == other.tds() && index_ == other.index_;
    }

    bool operator!=(const Cell& other) const
    {
      return !(*this == other);
    }

    bool operator<(const Cell& other) const
    {
      return (tds_ == other.tds()) ? (index_ < other.index_) : (tds_ < other.tds());
    }

    Cell_index index() const
    {
      return index_;
    }

    decltype(auto) storage() {
      return tds()->cell_storage()[index()];
    }

    decltype(auto) storage() const {
      return tds()->cell_storage()[index()];
    }

    Vertex_handle vertex(int i) const
    {
      return Vertex_handle(tds(), storage().ivertices[i]);
    }

    int index(Vertex_handle v) const
    {
      for (int i = 0; i < 4; ++i) {
        if (v.index() == storage().ivertices[i]) {
          return i;
        }
      }
      return -1; // Not found
    }

    bool has_vertex(Vertex_handle v) const
    {
      for (int i = 0; i < 4; ++i) {
        if (v.index() == storage().ivertices[i]) {
          return true;
        }
      }
      return false;
    }

    bool has_vertex(Vertex_handle v, int & i) const
    {
      for (i = 0; i < 4; ++i) {
        if (v.index() == storage().ivertices[i]) {
          return true;
        }
      }
      return false;
    }

    Cell_handle neighbor(int i) const
    {
      return Cell_handle(tds(), storage().ineighbors[i]);
    }

    int index(Cell_handle n) const
    {
      for (int i = 0; i < 4; ++i) {
        if (n.index() == storage().ineighbors[i]) {
          return i;
        }
      }
      return -1; // Not found
    }

    bool has_neighbor(Cell_handle n) const
    {
      for (int i = 0; i < 4; ++i) {
        if (n.index() == storage().ineighbors[i]) {
          return true;
        }
      }
      return false;
    }

    bool has_neighbor(Cell_handle n, int & i) const
    {
      for (i = 0; i < 4; ++i) {
        if (n.index() == storage().ineighbors[i]) {
          return true;
        }
      }
      return false;
    }

    void set_vertex(int i, Vertex_handle v)
    {
      storage().ivertices[i] = v.index();
    }

    void set_neighbor(int i, Cell_handle n)
    {
      storage().ineighbors[i] = n.index();
    }

    void set_vertices(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
    {
      storage().ivertices[0] = v0.index();
      storage().ivertices[1] = v1.index();
      storage().ivertices[2] = v2.index();
      storage().ivertices[3] = v3.index();
    }

    void set_neighbors(Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3)
    {
      storage().ineighbors[0] = n0.index();
      storage().ineighbors[1] = n1.index();
      storage().ineighbors[2] = n2.index();
      storage().ineighbors[3] = n3.index();
    }

    // CHECKING
    bool is_valid(bool  = false,
                  int   = 0) const
    {
      return true;
    }

    TDS_data& tds_data() {
      return tds()->cell_tds_data_pmap()[index()];
    }

    const TDS_data& tds_data() const {
      return tds()->cell_tds_data_pmap()[index()];
    }

    template <typename TDS>
    struct Rebind_TDS {
      using Other = Cell<TDS>;
    };

    TDS_3* tds() const // AF: constness issue
    {
      return tds_;
    }

    friend std::ostream& operator<<(std::ostream& os, const Cell&)
    {
      // Non combinatorial information. Default = nothing.
      return os;
    }

    friend std::istream& operator>>(std::istream& is, Cell)
    {
      return is;
    }

    TDS_3* tds_;
    Cell_index index_;
  };


  // Specialization for void.
  template <>
  class Cell<void>
  {
  public:
    using Triangulation_data_structure = internal::Dummy_tds_3;
    using Vertex_handle = Triangulation_data_structure::Vertex_handle;
    using Cell_handle = Triangulation_data_structure::Cell_handle;

    template <typename TDS2>
    struct Rebind_TDS { using Other = Cell<TDS2>; };
  };



  template <typename TDS_3 = void>
  struct Vertex{

    using Triangulation_data_structure = TDS_3;
    using Cell_handle   = typename TDS_3::Cell_handle;
    using Cell_index    = typename TDS_3::Cell_index;
    using Vertex_index  = typename TDS_3::Vertex_index;
    using Index         = Vertex_index;
    using Vertex_handle [[maybe_unused]] = typename TDS_3::Vertex_handle;

    struct Storage {
      Cell_index icell;
    };

    Vertex()
      : tds_(nullptr), index_()
    {}

    Vertex(TDS_3* tds, Vertex_index index)
      : tds_(tds), index_(index)
    {}

    bool operator==(const Vertex& other) const
    {
      return tds() == other.tds() && index_ == other.index_;
    }

    bool operator!=(const Vertex& other) const
    {
      return !(*this == other);
    }

    bool operator<(const Vertex& other) const
    {
      return (tds() == other.tds()) ? (index_ < other.index_) : (tds() < other.tds());
    }

    Vertex_index index() const
    {
      return index_;
    }

    decltype(auto) storage() {
      return tds()->vertex_storage()[index()];
    }

    decltype(auto) storage() const {
      return tds()->vertex_storage()[index()];
    }

    Cell_handle cell() const
    {
      return Cell_handle(tds(), storage().icell);
    }

    void set_cell(Cell_handle c)
    {
      storage().icell = c.index();
    }

    bool is_valid(bool = false, int = 0) const
    {
      return true;
    }

    template <typename TDS>
    struct Rebind_TDS {
      using Other = Vertex<TDS>;
    };

    TDS_3* tds() const // AF: constness issue
    {
      return tds_;
    }

    TDS_3* tds_;
    Vertex_index index_;

    friend std::ostream& operator<<(std::ostream& os, const Vertex&)
    {
      return os;
    }

    friend std::istream& operator>>(std::istream& is, Vertex)
    {
      return is;
    }
  };

  template <typename GT, typename Vb = Vertex<>>
  struct VertexWithPoint
  : public Vb
  {
  public:
    using Vb::Vb;
    using Point = typename GT::Point_3;
    struct Storage : public Vb::Storage {
      typename GT::Point_3 point;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
      using Other = VertexWithPoint<GT,Vb2>;
    };

    decltype(auto) storage() {
      return this->tds()->vertex_storage()[this->index()];
    }

    decltype(auto) storage() const {
      return this->tds()->vertex_storage()[this->index()];
    }

    const Point& point() const
    {
      return storage().point;
    }

    void set_point(const Point& p)
    {
      storage().point = p;
    }

    friend std::ostream& operator<<(std::ostream& os, const VertexWithPoint& v)
    {
      os << v.point();
      return os;
    }

    friend std::istream& operator>>(std::istream& is, VertexWithPoint v)
    {
      typename GT::Point_3 p;
      if(is >> p) {
        v.set_point(p);
      }
      return is;
    }
  };


  template <typename GT, typename Vb = Vertex<>>
  struct VertexWithWeightedPoint
  : public Vb
  {
  public:
    using Vb::Vb;
    using Point = typename GT::Weighted_point_3;
    struct Storage : public Vb::Storage {
      typename GT::Weighted_point_3 point;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
      using Other = VertexWithWeightedPoint<GT,Vb2>;
    };

    decltype(auto) storage() {
      return this->tds()->vertex_storage()[this->index()];
    }

    decltype(auto) storage() const {
      return this->tds()->vertex_storage()[this->index()];
    }

    const Point& point() const
    {
      return storage().point;
    }

    void set_point(const Point& p)
    {
      storage().point = p;
    }

    friend std::ostream& operator<<(std::ostream& os, const VertexWithWeightedPoint& v)
    {
      os << v.point();
      return os;
    }

    friend std::istream& operator>>(std::istream& is, VertexWithWeightedPoint v)
    {
      typename GT::Weighted_point_3 p;
      if(is >> p) {
        v.set_point(p);
      }
      return is;
    }
  };


  template <typename Vb>
  class Vertex4Hierarchy
    : public Vb
  {
  public:
    using Vb::Vb; // inherit constructors
    using TDS = typename Vb::Triangulation_data_structure;
    using Vertex_handle = typename TDS::Vertex_handle;

    struct Storage : public Vb::Storage {
      Vertex_handle up;
      Vertex_handle down;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Vb2 = typename Vb::template Rebind_TDS<TDS2>::Other;
      using Other = Vertex4Hierarchy<Vb2>;
    };

    decltype(auto) storage() {
      return this->tds()->vertex_storage()[this->index()];
    }

    decltype(auto) storage() const {
      return this->tds()->vertex_storage()[this->index()];
    }

    Vertex_handle up()   const { return storage().up; }
    Vertex_handle down() const { return storage().down; }
    void set_up(Vertex_handle u)   { storage().up=u; }
    void set_down(Vertex_handle d) { storage().down=d; }

  };

  // Specialization for void.
  template <>
  class Vertex<void>
  {
  public:
    using Triangulation_data_structure = internal::Dummy_tds_3;
    using Cell_handle = Triangulation_data_structure::Cell_handle;
    using Vertex_handle = Triangulation_data_structure::Vertex_handle;

    template <typename TDS2>
    struct Rebind_TDS { using Other = Vertex<TDS2>; };
  };

  // AF : factorize as also in Surface_mesh and Point_set_3
  template<typename Derived_id>
  class Index
  {
  public:

    using size_type = std::uint32_t;
    static constexpr size_type invalid_index = (std::numeric_limits<size_type>::max)();

    /// Constructor. %Default construction creates an invalid index.
    /// We write <a href="https://en.cppreference.com/w/cpp/types/numeric_limits">
    /// <tt>(std::numeric_limits<size_type>::max)()</tt></a>.
    explicit Index(size_type _idx=invalid_index) : idx_(_idx) {}

    /// Get the underlying index of this index
    operator size_type() const { return idx_; }

    /// reset index to be invalid (index=(std::numeric_limits<size_type>::max)())
    void reset() { idx_ = invalid_index; }

    /// return whether the index is valid, i.e., the index is not equal to `%std::numeric_limits<size_type>::max()`.
    bool is_valid() const {
      return idx_ != invalid_index;
    }

    // Compatibility with OpenMesh handle
    size_type idx() const {
      return idx_;
    }
    size_type& idx() {
      return idx_;
    }
    // For convenience
    size_type id() const {
      return idx_;
    }
    size_type& id() {
      return idx_;
    }

    /// increments the internal index. This operation does not
    /// guarantee that the index is valid or undeleted after the
    /// increment.
    Index& operator++() { ++idx_; return *this; }
    /// decrements the internal index. This operation does not
    /// guarantee that the index is valid or undeleted after the
    /// decrement.
    Index& operator--() { --idx_; return *this; }

    /// increments the internal index. This operation does not
    /// guarantee that the index is valid or undeleted after the
    /// increment.
    Index operator++(int) { Index tmp(*this); ++idx_; return tmp; }
    /// decrements the internal index. This operation does not
    /// guarantee that the index is valid or undeleted after the
    /// decrement.
    Index operator--(int) { Index tmp(*this); --idx_; return tmp; }

    Index operator+=(std::ptrdiff_t n) { idx_ = size_type(std::ptrdiff_t(idx_) + n); return *this; }

    template<class T> friend bool operator==(const T&, const Derived_id&) = delete;
    template<class T> friend bool operator!=(const T&, const Derived_id&) = delete;
    template<class T> friend bool operator< (const T&, const Derived_id&) = delete;
    template<class T> friend bool operator<=(const T&, const Derived_id&) = delete;
    template<class T> friend bool operator> (const T&, const Derived_id&) = delete;
    template<class T> friend bool operator>=(const T&, const Derived_id&) = delete;

    template<class T> friend bool operator==(const Derived_id&, const T&) = delete;
    template<class T> friend bool operator!=(const Derived_id&, const T&) = delete;
    template<class T> friend bool operator< (const Derived_id&, const T&) = delete;
    template<class T> friend bool operator<=(const Derived_id&, const T&) = delete;
    template<class T> friend bool operator> (const Derived_id&, const T&) = delete;
    template<class T> friend bool operator>=(const Derived_id&, const T&) = delete;

    friend bool operator==(const Derived_id& _lhs, const Derived_id& _rhs) { return _lhs.id() == _rhs.id(); }
    friend bool operator!=(const Derived_id& _lhs, const Derived_id& _rhs) { return !(_lhs == _rhs); }
    friend bool operator< (const Derived_id& _lhs, const Derived_id& _rhs) { return _lhs.id() < _rhs.id(); }
    friend bool operator<=(const Derived_id& _lhs, const Derived_id& _rhs) { return _lhs.id() <= _rhs.id(); }
    friend bool operator> (const Derived_id& _lhs, const Derived_id& _rhs) { return _lhs.id() > _rhs.id(); }
    friend bool operator>=(const Derived_id& _lhs, const Derived_id& _rhs) { return _lhs.id() >= _rhs.id(); }
    friend bool operator==(const Derived_id& _lhs, const std::nullptr_t&) { return _lhs.id() == invalid_index; }
    friend bool operator!=(const Derived_id& _lhs, std::nullptr_t&) { return _lhs.id() != invalid_index; }

    friend std::ostream& operator<<(std::ostream& os, const Derived_id& idx)
    {
      return (os << idx.output_prefix() << idx.id() );
    }

  protected:
    size_type idx_ = invalid_index;
  }; // end class Index

  template <class T>
  std::size_t hash_value(const Index<T>&  i)
  {
    return i.id();
  }

  template <typename T, typename Index_type_, typename Element_container>
  class Index_handle {
    using Element = T;
    using size_type = typename Element_container::size_type;
    using Proxy = boost::stl_interfaces::proxy_arrow_result<Element>;
  public:
    using Index_type = Index_type_;
    using value_type = Element;
    using reference = Element;
    using pointer = Proxy;

    Index_handle() = default;

    Index_handle(Element_container* container, size_type idx)
      : cont_(container), idx_(idx) {}

    Element operator*() const {
      return Element(cont_, Index_type(idx_));
    }

    Proxy operator->() const {
      return Proxy{this->operator*()};
    }

    Index_type index() const {
      return Index_type{idx_};
    }

    auto container() const {
      return cont_;
    }

    bool operator==(const Index_handle& other) const {
      CGAL_assertion(container() == nullptr || other.container() == nullptr ||
                     container() == other.container());
      return index() == other.index();
    }

    bool operator!=(const Index_handle& other) const {
      return !(*this == other);
    }

    bool operator<(const Index_handle& other) const {
      return (container() == other.container()) ? (index() < other.index()) : (container() < other.container());
    }

    bool operator==( std::nullptr_t ) const {
      return cont_ == nullptr && idx_ == Index_type::invalid_index;
    }

    bool operator!=( std::nullptr_t ) const {
      return !(*this == nullptr);
    }

    friend std::ostream& operator<<(std::ostream& os, const Index_handle& h)
    {
      return os << "#" << h.index();
    }

    friend std::size_t hash_value(const Index_handle& h) {
      return static_cast<std::size_t>(h.index().id());
    }

  private:
    Element_container* cont_ = nullptr;
    size_type idx_ = Index_type::invalid_index;
  }; // end class Index_handle

} // end namespace CGAL

namespace std {
  template <typename T, typename Index_type, typename Element_container>
  struct hash<CGAL::Index_handle<T, Index_type, Element_container>> {
    using Handle = CGAL::Index_handle<T, Index_type, Element_container>;
    std::size_t operator()(const Handle& h) const {
      return hash_value(h);
    }
  };
}

namespace CGAL {

  template <typename Handle_, typename Element_container>
  class Index_iterator // make it derive from Handle_
      : public boost::stl_interfaces::v1::proxy_iterator_interface<Index_iterator<Handle_, Element_container>,
                                                                  std::random_access_iterator_tag,
                                                                  typename Handle_::value_type>
  {
    using Index = typename Handle_::Index_type;

    using Facade = boost::stl_interfaces::v1::proxy_iterator_interface<Index_iterator<Handle_, Element_container>,
                                                                      std::random_access_iterator_tag,
                                                                      typename Handle_::value_type>;

  public:
    using value_type = typename Handle_::value_type;
    using reference = value_type;
    using Concurrency_tag = typename Element_container::Concurrency_tag;

    Index_iterator()
        : idx_()
        , cont_(nullptr) {}

    Index_iterator(const Index& h, const Element_container* m)
        : idx_(h)
        , cont_(const_cast<Element_container*>(m)) { // AF: todo make const_cast safe
      if(cont_ && cont_->has_garbage()) {
        while(cont_->has_valid_index(idx_) && cont_->is_removed(idx_))
          ++idx_;
      }
    }

    auto handle() const {
      static_assert(std::is_base_of_v<Facade, Index_iterator<Handle_, Element_container>>);

      CGAL_assertion(cont_ != nullptr);
      CGAL_assertion(cont_->has_valid_index(idx_));
      return Handle_(cont_, idx_.id());
    }

    operator Handle_() const { return handle(); }

    reference operator*() const { return value_type{cont_, idx_}; }

    typename Handle_::pointer operator->() const { return typename Handle_::pointer{value_type{cont_, idx_}}; }

    using Facade::operator++;
    Index_iterator& operator++() {
      ++idx_;
      CGAL_assertion(cont_ != nullptr);

      if(cont_->has_garbage())
        while(cont_->has_valid_index(idx_) && cont_->is_removed(idx_))
          ++idx_;
      return *this;
    }

    using Facade::operator--;
    Index_iterator& operator--() {
      --idx_;
      CGAL_assertion(cont_ != nullptr);
      if(cont_->has_garbage())
        while(cont_->has_valid_index(idx_) && cont_->is_removed(idx_))
          --idx_;
      return *this;
    }

    Index_iterator& operator+=(std::ptrdiff_t n) {
      CGAL_assertion(cont_ != nullptr);

      if(cont_->has_garbage()) {
        if(n > 0)
          for(std::ptrdiff_t i = 0; i < n; ++i)
            this->operator++();
        else
          for(std::ptrdiff_t i = 0; i < -n; ++i)
            this->operator--();
      } else
        idx_ += n;
      return *this;
    }

    std::ptrdiff_t operator-(const Index_iterator& other) const {
      if(cont_->has_garbage()) {
        bool forward = (other.idx_ > idx_);

        std::ptrdiff_t out = 0;
        Index_iterator it = *this;
        while(!(it == other)) {
          if(forward) {
            ++it;
            ++out;
          } else {
            --it;
            --out;
          }
        }
        return out;
      }

      // else
      return std::ptrdiff_t(other.idx_) - std::ptrdiff_t(this->idx_);
    }

    bool operator==(const Index_iterator& other) const { return this->idx_ == other.idx_; }

  private:
    Index idx_;
    Element_container* cont_;
  }; // end class Index_iterator

  template <class I, class T, class ConcurrencyTag = Sequential_tag>
  struct Property_map : Properties::Property_map_base<I, T, Property_map<I, T>, ConcurrencyTag>
  {
    using Base = Properties::Property_map_base<I, T, Property_map<I, T>, ConcurrencyTag>;

    using Base::Base;
  };

  template <typename Index_type,
            typename Element_type,
            typename Storage_type,
            typename Container,
            typename ConcurrencyTag,
            typename Container::size_type& (*free_list_next_function_)(Storage_type&),
            char prefix>
  struct Indexed_container
  {
    using Self = Indexed_container<Index_type,
                                   Element_type,
                                   Storage_type,
                                   Container,
                                   ConcurrencyTag,
                                   free_list_next_function_,
                                   prefix>;
    using Handle = Index_handle<Element_type, Index_type, Container>;
    using size_type = typename Container::size_type;
    using Concurrency_tag = ConcurrencyTag;

    static constexpr bool is_parallel = std::is_convertible_v<Concurrency_tag, Parallel_tag>;
#ifdef CGAL_LINKED_WITH_TBB
    struct Freelist_handler {
      Freelist_handler(size_type invalid_index = Index_type::invalid_index) : freelist(invalid_index) {}

      size_type freelist;
      size_type number_of_removed_elements = 0;
    };
    using free_list_type = std::conditional_t<is_parallel,
                                              tbb::enumerable_thread_specific<Freelist_handler>,
                                              size_type>;
#else
    static_assert(!is_parallel,
                  "In CGAL triangulations, `Parallel_tag` can only be used with the Intel TBB library. "
                  "Make TBB available in the build system and then define the macro `CGAL_LINKED_WITH_TBB`.");
    using free_list_type = size_type;
#endif

    Properties::Property_container<Self, Index_type, Concurrency_tag> properties_;
    size_type nb_of_removed_elements_ = 0;
    free_list_type freelist_{Index_type::invalid_index};
    size_type anonymous_property_nb = 0;
    static constexpr bool recycle_ = true;
    bool garbage_ = false;
    Property_map<Index_type, Storage_type, Concurrency_tag> storage_ =
        add_property_map<Storage_type>(prefix + std::string(":storage"), Storage_type()).first;
    Property_map<Index_type, bool, Concurrency_tag> removed_ =
        add_property_map<bool>(prefix + std::string(":removed"), false).first;

    template <typename Key, typename T>
    struct Get_property_map {
      using type = Property_map<Key, T, Concurrency_tag>;
    };

    using Time_stamper = CGAL::Time_stamper_impl<Element_type>;
    template <typename U> using EraseCounterStrategy =
      internal::Erase_counter_strategy<internal::has_increment_erase_counter<U>::value>;

    template<class T>
    std::pair<Property_map<Index_type, T, ConcurrencyTag>, bool>
    add_property_map(std::string name=std::string(), const T t=T()) {
      if(name.empty()){
        std::ostringstream oss;
        oss << "anonymous-property-" << anonymous_property_nb++;
        name = std::string(oss.str());
      }
      return properties_.template add<T>(name, t);
    }

    size_type size() const
    {
      return static_cast<size_type>(properties_.size());
    }

    void reserve(size_type n)
    {
      properties_.reserve(n);
    }

    void clear() {
      properties_.resize(0);
      properties_.shrink_to_fit();

      freelist_ = free_list_type{Index_type::invalid_index};
      nb_of_removed_elements_ = 0;
      garbage_ = false;
    }

    size_type& free_list()
    {
      if constexpr(is_parallel) {
        return freelist_.local().freelist; // TBB
      } else {
        return freelist_; // Sequential
      }
    }

    size_type& local_number_of_removed_elements()
    {
      if constexpr(is_parallel) {
        return freelist_.local().number_of_removed_elements; // TBB
      } else {
        return nb_of_removed_elements_; // Sequential
      }
    }

    size_type number_of_removed_elements() const {
      if constexpr(is_parallel) {
        size_type result = 0;
        for(const auto& f : freelist_) {
          result += f.number_of_removed_elements;
        }
        return result;
      } else {
        return nb_of_removed_elements_; // Sequential
      }
    }

    Index_type create(Container* container)
    {
      size_type& freelist_ = free_list();
      Index_type idx{freelist_};
      Element_type elt(container, idx);
      if(recycle_ && (freelist_ != Index_type::invalid_index)){
        freelist_ = free_list_next_function_(storage_[idx]);
        --local_number_of_removed_elements();
        removed_[idx] = false;
        const auto ec = EraseCounterStrategy<Element_type>::erase_counter(&elt);
        properties_.reset(idx);
        EraseCounterStrategy<Element_type>::restore_erase_counter(&elt, ec);
      } else {
        elt = Element_type(container, idx);
        idx = Index_type{static_cast<size_type>(properties_.push_back())};
      }
      Time_stamper::restore_timestamp(&elt, elt.index().id());
      return idx;
    }

    void remove(Handle ch, Container* container)
    {
      size_type& freelist_ = free_list();
      Index_type idx = ch->index();
      if constexpr(! is_parallel){
        if(idx.idx() == properties_.size()-1){
          properties_.pop_back();
          return;
        }
      }
      removed_[idx] = true; ++local_number_of_removed_elements(); garbage_ = true;
      free_list_next_function_(storage_[idx]) = Index_type{freelist_};
      freelist_ = static_cast<size_type>(idx);
      Element_type elt(container, idx);
      EraseCounterStrategy<Element_type>::increment_erase_counter(elt);
    }

    bool is_valid_index(Index_type idx) const {
      return idx.id() < size() && !removed_[idx];
    }

    bool has_garbage() const { return garbage_; }
  }; // end class Indexed_container

  template <typename Vb = Vertex<>, typename Cb = Cell<>, class ConcurrencyTag = Sequential_tag>
  struct Indexed_storage
  {

    using Self = Indexed_storage<Vb,Cb, ConcurrencyTag>;
    using TDS = Triangulation_data_structure_3<Vb, Cb, ConcurrencyTag, Index_tag>;
    using Concurrency_tag  = ConcurrencyTag;
    using Storage_tag = Index_tag;

    TDS& tds()
    {
      static_assert(std::is_base_of_v<Self, TDS>,
                    "Indexed_storage must be a base class of TDS");
      return *static_cast<TDS*>(this);
    }

    const TDS& tds() const
    {
      static_assert(std::is_base_of_v<Self, TDS>,
                    "Indexed_storage must be a base class of TDS");
      return *static_cast<const TDS*>(this);
    }

    /// The type used to represent an index.
    using size_type = std::uint32_t;
    using difference_type = std::ptrdiff_t;

    // Tools to change the Vertex and Cell types of the TDS.
    template < typename Vb2 >
    struct Rebind_vertex {
      using Other = Indexed_storage<Vb2, Cb, ConcurrencyTag>;
    };

    template < typename Cb2 >
    struct Rebind_cell {
      using Other = Indexed_storage<Vb, Cb2, ConcurrencyTag>;
    };

    // Put this TDS inside the Vertex and Cell types.
    using Vertex = typename Vb::template Rebind_TDS<Self>::Other;
    using Cell = typename Cb::template Rebind_TDS<Self>::Other;

    struct Vertex_index: public Index<Vertex_index>
    {
      using Index<Vertex_index>::Index; // inherit constructors
      auto output_prefix() const { return 'v'; }
    };

    struct Cell_index: public Index<Cell_index>
    {
      using Index<Cell_index>::Index; // inherit constructors
      auto output_prefix() const { return 'c'; }
    };

    using Cell_handle = Index_handle<Cell, Cell_index, Self>;
    using Vertex_handle = Index_handle<Vertex, Vertex_index, Self>;

    using Facet = std::pair<Cell_handle, int>;
    using Edge = Triple<Cell_handle, int, int>;

    class Cell_data {
      unsigned char conflict_state;
    public:
      Cell_data() : conflict_state(0) {}

      void clear()            { conflict_state = 0; }
      void mark_in_conflict() { conflict_state = 1; }
      void mark_on_boundary() { conflict_state = 2; }
      void mark_processed()   { conflict_state = 1; }

      bool is_clear()       const { return conflict_state == 0; }
      bool is_in_conflict() const { return conflict_state == 1; }
      bool is_on_boundary() const { return conflict_state == 2; }
      bool processed() const { return conflict_state == 1; }
    };

    using TDS_data = Cell_data;

    using Vertex_storage = typename Vertex::Storage;
    using Cell_storage = typename Cell::Storage;
    using Vertex_storage_property_map = Property_map<Vertex_index, Vertex_storage, ConcurrencyTag>;
    using Cell_storage_property_map = Property_map<Cell_index, Cell_storage, ConcurrencyTag>;

    static size_type& vertex_free_list_next_function(Vertex_storage& s)
    {
      return s.icell.id();
    }
    static size_type& cell_free_list_next_function(Cell_storage& s)
    {
      return s.ivertices[0].id();
    }

    using Vertex_container = Indexed_container<Vertex_index, Vertex, Vertex_storage, Self, Concurrency_tag,
        &vertex_free_list_next_function, 'v'>;
    using Cell_container = Indexed_container<Cell_index, Cell, Cell_storage, Self, Concurrency_tag,
        &cell_free_list_next_function, 'c'>;

    Cell_storage_property_map& cell_storage() { return cell_container().storage_; }
    const Cell_storage_property_map& cell_storage() const { return cell_container().storage_; }

    Vertex_storage_property_map& vertex_storage() { return vertex_container().storage_; }
    const Vertex_storage_property_map& vertex_storage() const { return vertex_container().storage_; }

    auto cell_tds_data_pmap() { return cell_data_; }
    auto cell_tds_data_pmap() const { return cell_data_; }

    size_type num_vertices() const { return static_cast<size_type>(vertex_container().size()); }
    size_type num_cells() const { return static_cast<size_type>(cell_container().size()); }

    int dimension() const { return dimension_; }

    void set_dimension(int n) { dimension_ = n; }

    template <typename Index_type>
    auto& container()
    {
      if constexpr (std::is_same_v<Index_type, Vertex_index>) {
        return vertex_container();
      } else if constexpr (std::is_same_v<Index_type, Cell_index>) {
        return cell_container();
      } else {
        static_assert(std::is_same_v<Index_type, void>, "Invalid Index type");
      }
    }

    template <typename Index_type>
    const auto& container() const
    {
      if constexpr (std::is_same_v<Index_type, Vertex_index>) {
        return vertex_container();
      } else if constexpr (std::is_same_v<Index_type, Cell_index>) {
        return cell_container();
      } else {
        static_assert(std::is_same_v<Index_type, void>, "Invalid Index type");
      }
    }

    Vertex_handle create_vertex()
    {
      return Vertex_handle{this, vertex_container().create(this)};
    }

    Cell_handle create_cell()
    {
      return Cell_handle{this, cell_container().create(this)};
    }

    Vertex_handle create_vertex(const Vertex& v)
    {
      Vertex_handle new_v = create_vertex();
      new_v->storage() = v.storage();
      return new_v;
    }

    Cell_handle create_cell(const Cell& c)
    {
      Cell_handle new_c = create_cell();
      new_c->storage() = c.storage();
      return new_c;
    }

    // AF:  What about the equivalent to
    // https://doc.cgal.org/latest/TDS_3/classTriangulationDataStructure__3.html#a1432860206073c24ca43dbbdfb13b26e

    Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3)
    {
      Cell_handle c = create_cell();
      c->set_vertices(v0, v1, v2, v3);
      c->set_neighbors(Cell_handle(), Cell_handle(), Cell_handle(), Cell_handle());
      return c;
    }

    Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2, Vertex_handle v3,
                            Cell_handle n0, Cell_handle n1,
                            Cell_handle n2, Cell_handle n3)
    {
      Cell_handle c = create_cell();
      c->set_vertices(v0, v1, v2, v3);
      c->set_neighbors(n0, n1, n2, n3);
      return c;
    }

    Cell_handle create_face()
    {
      CGAL_precondition(dimension()<3);
      return create_cell();
    }

    Cell_handle create_face(Vertex_handle v0, Vertex_handle v1,
                            Vertex_handle v2)
    {
      CGAL_precondition(dimension()<3);
      return create_cell(v0, v1, v2, Vertex_handle());
    }

    void delete_vertex(Vertex_handle vh)
    {
      vertex_container().remove(vh, this);
    }

    void delete_cell(Cell_handle ch)
    {
      cell_container().remove(ch, this);
    }

    void reserve(size_type n_vertices, size_type n_cells)
    {
      vertex_container().reserve(n_vertices);
      cell_container().reserve(n_cells);
    }

    void clear()
    {
      clear_without_removing_property_maps();
      remove_all_property_maps();
    }

    void clear_without_removing_property_maps()
    {
      vertex_container().clear();
      cell_container().clear();
      dimension_ = -2;
    }

    void remove_all_property_maps()
    {
      remove_property_maps<Vertex_index>(2);
      remove_property_maps<Cell_index>(3);
    }

    size_type number_of_removed_vertices() const
    {
      return vertex_container().number_of_removed_elements();
    }

    size_type number_of_removed_cells() const
    {
      return cell_container().number_of_removed_elements();
    }

    size_type number_of_vertices() const
    {
      return num_vertices() - number_of_removed_vertices();
    }

    size_type number_of_cells() const
    {
      return num_cells() - number_of_removed_cells();
    }

    bool is_valid_vertex_index(Vertex_index idx) const {
      return vertex_container().is_valid_index(idx);
    }

    bool is_valid_cell_index(Cell_index idx) const {
      return cell_container().is_valid_index(idx);
    }

    bool is_vertex(Vertex_handle v) const
    {
      return this == v->tds() && is_valid_vertex_index(v->index());
    }

    bool is_valid_cell_handle(Cell_handle c) const
    {
      return this == c->tds()  && is_valid_cell_index(c->index());
    }

    bool is_cell( Cell_handle c ) const
      // returns false when dimension <3
    {
      if (dimension() < 3)
        return false;

      return is_valid_cell_handle(c);
    }

    /// adds a property map named `name` with value type `T` and default `t`
    /// for index type `I`. Returns the property map together with a Boolean
    /// that is `true` if a new map was created. In case it already exists
    /// the existing map together with `false` is returned.


    template<class I, class T>
    std::pair<Property_map<I, T>, bool>
    add_property_map(std::string name=std::string(), const T t=T()) {
      return container<I>().template add_property_map<T>(name, t);
    }

    /// returns an optional property map named `name` with key type `I` and value type `T`.
    template <class I, class T>
    std::optional<Property_map<I, T>> property_map(const std::string& name) const
    {
      return container<I>().properties_.template get<T>(name);
    }


    /// removes property map `p`. The memory allocated for that property map is freed.
    template<class I, class T>
    void remove_property_map(Property_map<I, T>& p)
    {
      container<I>().properties_.template remove<T>(p);
    }

    /// removes all property maps for index type `I` added by a call to `add_property_map<I>()`.
    /// The memory allocated for those property maps is freed.
    template<class I>
    void remove_property_maps(int nb_of_properties_to_keep)
    {
      container<I>().properties_.resize_property_array(nb_of_properties_to_keep);
    }

    /// @cond CGAL_DOCUMENT_INTERNALS
    /// returns the std::type_info of the value type of the
    /// property identified by `name`.  `typeid(void)` if `name`
    /// does not identify any property.
    ///
    /// @tparam I The key type of the property.

    template<class I>
    const std::type_info& property_type(const std::string& name)
    {
      return container<I>().properties_.get_type(name);
    }
    /// @endcond

    /// returns a vector with all strings that describe properties with the key type `I`.
    /// @tparam I The key type of the properties.
    template<class I>
    std::vector<std::string> properties() const
    {
      return container<I>().properties_.properties();
    }

    Indexed_storage() = default;
    Indexed_storage(const Indexed_storage&) = default;
    Indexed_storage(Indexed_storage&&) = default;
    Indexed_storage& operator=(const Indexed_storage&) = default;
    Indexed_storage& operator=(Indexed_storage&& is) = default;

    void swap(Indexed_storage& other)
    {
      std::swap(*this, other);
    }

    ~Indexed_storage()
    {
#if 0
      std::cout << "sizeof(Vertex) " << " " << sizeof(Vertex_storage) << std::endl;
      std::cout << "sizeof(Cell)   " << " " << sizeof(Cell_storage) << std::endl;
      std::cout << "vertex capacity: "  << vprops_.capacity() << std::endl;
      std::cout << "cell capacity: "  << cprops_.capacity() << std::endl;
      std::cout << sizeof(Vertex_storage)  * vprops_.capacity() << " bytes for vertices" << std::endl;
      std::cout << sizeof(Cell_storage)  * cprops_.capacity() << " bytes for cells" << std::endl;
#endif
    }

    void set_adjacency(Cell_handle c0, int i0,
                       Cell_handle c1, int i1) const
    {
      CGAL_assertion(i0 >= 0 && i0 <= dimension());
      CGAL_assertion(i1 >= 0 && i1 <= dimension());
      CGAL_assertion(c0 != c1);
      c0->set_neighbor(i0,c1);
      c1->set_neighbor(i1,c0);
    }

    bool has_garbage() const { return vertex_container().has_garbage() || cell_container().has_garbage(); }

    /// returns whether the index of vertex `v` is valid, that is within the current array bounds.
    bool has_valid_index(Vertex_index v) const
    {
      return (v.id() < num_vertices());
    }

    bool has_valid_index(Cell_index c) const
    {
      return (c.id() < num_cells());
    }

    /// returns whether vertex `v` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Vertex_index v) const
    {
      return vertex_container().removed_[v];
    }

    bool is_removed(Cell_index c) const
    {
      return cell_container().removed_[c];
    }
    //------------------------------------------------------ iterator types
    using Vertex_iterator = Index_iterator<Vertex_handle, Self>;
    using Vertex_range = Iterator_range<Vertex_iterator>;

    Vertex_iterator vertices_begin() const
    {
      return Vertex_iterator(Vertex_index(0), this);
    }

    /// End iterator for vertices.
    Vertex_iterator vertices_end() const
    {
      return Vertex_iterator(Vertex_index(num_vertices()), this);
    }

    Vertex_range vertices() const {
      return make_range(vertices_begin(), vertices_end());
    }


    using Cell_iterator = Index_iterator<Cell_handle, Self>;
    using Cell_range = Iterator_range<Cell_iterator>;

    Cell_iterator cells_begin() const
    {
      return Cell_iterator(Cell_index(0), this);
    }

    /// End iterator for cells.
    Cell_iterator cells_end() const
    {
      return Cell_iterator(Cell_index(num_cells()), this);
    }

    Cell_range cells() const {
      return make_range(cells_begin(), cells_end());
    }


    friend class internal::Triangulation_ds_facet_iterator_3<Self>;
    friend class internal::Triangulation_ds_edge_iterator_3<Self>;

    friend class internal::Triangulation_ds_cell_circulator_3<Self>;
    friend class internal::Triangulation_ds_facet_circulator_3<Self>;



    using Facet_iterator = internal::Triangulation_ds_facet_iterator_3<Self>;
    using Edge_iterator = internal::Triangulation_ds_edge_iterator_3<Self>;
    using Facets = Iterator_range<Facet_iterator>;
    using Edges = Iterator_range<Edge_iterator>;

    using Cell_circulator = internal::Triangulation_ds_cell_circulator_3<Self>;
    using Facet_circulator = internal::Triangulation_ds_facet_circulator_3<Self>;

    Edge_iterator edges_begin() const
    {
      if ( dimension() < 1 )
        return edges_end();
      return Edge_iterator(this);
    }

    Edge_iterator edges_end() const
    {
      return Edge_iterator(this,1);
    }

    Edges edges() const
    {
      return Edges(edges_begin(), edges_end());
    }

    Facet_iterator facets_begin() const
    {
      if ( dimension() < 2 )
        return facets_end();
      return Facet_iterator(this);
    }

    Facet_iterator facets_end() const
    {
      return Facet_iterator(this, 1);
    }

    Facets facets() const
    {
      return Facets(facets_begin(), facets_end());
    }

    Cell_circulator incident_cells(const Edge & e) const
    {
      CGAL_precondition( dimension() == 3 );
      return Cell_circulator(e);
    }

    Facet_circulator incident_facets(const Edge & e) const
    {
      CGAL_precondition( dimension() == 3 );
      return Facet_circulator(e);
    }

    template <class TDS_src>
    Vertex_handle copy_tds(const TDS_src & src, typename TDS_src::Vertex_handle vert)
    {
      internal::TDS_3::Default_index_vertex_converter<typename TDS_src::Vertex,Vertex> setv;
      internal::TDS_3::Default_index_cell_converter<typename TDS_src::Cell,Cell>  setc;
      return tds().copy_tds(src, vert, setv, setc);
    }

  protected:
    Vertex_container& vertex_container()
    {
      return vertex_container_;
    }
    const Vertex_container& vertex_container() const
    {
      return vertex_container_;
    }
    Cell_container& cell_container()
    {
      return cell_container_;
    }
    const Cell_container& cell_container() const
    {
      return cell_container_;
    }

    Vertex_container vertex_container_;
    Cell_container cell_container_;

    Property_map<Cell_index, Cell_data, Concurrency_tag> cell_data_ =
          cell_container().template add_property_map<Cell_data>("c:data").first;

    // in dimension i, number of vertices >= i+2
    // ( the boundary of a simplex in dimension i+1 has i+2 vertices )
    int dimension_ = -2;

  };

} /// namespace CGAL

#include <CGAL/Hidden_point_memory_policy.h>

namespace CGAL {


  template <typename GT,
            typename Cb = Cell<>>
  class Cell4Delaunay
  : public Cb
  {
    public:
    using Cb::Cb; // inherit constructors
    using Point_3 = typename GT::Point_3;

    template < typename TDS2 >
    struct Rebind_TDS {
      using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
      using Other = Cell4Delaunay<GT, Cb2>;
    };


    Point_3 circumcenter(const GT& gt) const
    {
      return gt.construct_circumcenter_3_object()(this->vertex(0)->point(),
                                                  this->vertex(1)->point(),
                                                  this->vertex(2)->point(),
                                                  this->vertex(3)->point());
    }

    Point_3 circumcenter() const
    {
      return circumcenter(GT());
    }

    void set_circumcenter(const Point_3&)
    {
      std::cerr << "ignore set_circumcenter() call";
    }

  };

  template <typename GT,
            typename Cb = Cell4Delaunay<GT>>
  class CellWithCircumcenter
  : public Cb
  {
    public:
    using Cb::Cb; // inherit constructors
    using Point = typename GT::Point_3;
    using TDS = typename Cb::Triangulation_data_structure;
    using Vertex_handle = typename TDS::Vertex_handle;
    using Cell_handle = typename TDS::Cell_handle;

    struct Storage : public Cb::Storage {
      std::optional<Point> circumcenter_;
      /*
      struct C {

        template <class T>
        bool operator==(const T& )const{
          return true;
        }
        template <class T>
        void reset(const T&) const
        {}
      operator bool() const
      {return true;}

      Point operator*() const
      {
        return Point();
      }

      };
      C circumcenter_;
      */
    };

    auto&& storage()
    {
      return this->tds()->cell_storage()[this->index()];
    }

    auto&& storage() const
    { return this->tds()->cell_storage()[this->index()]; }

    template < typename TDS2 >
    struct Rebind_TDS {
      using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
      using Other = CellWithCircumcenter<GT, Cb2>;
    };

    void invalidate_circumcenter()
    {
      if (storage().circumcenter_) {
          storage().circumcenter_.reset();
      }
    }

  void set_vertex(int i, Vertex_handle v)
  {
      invalidate_circumcenter();
      Cb::set_vertex(i, v);
  }

  void set_vertices()
  {
      invalidate_circumcenter();
      Cb::set_vertices();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
      invalidate_circumcenter();
      Cb::set_vertices(v0, v1, v2, v3);
  }

  void set_circumcenter(const Point& p) const
  {
    storage().circumcenter_ = std::make_optional<Point>(p);
  }

  const Point& circumcenter(const GT& gt = GT()) const
  {
      if (! storage().circumcenter_) {
        storage().circumcenter_ = std::make_optional<Point>(Point(this->Cb::circumcenter(gt)));
      } else {
        CGAL_expensive_assertion(
          this->Cb::circumcenter(gt) == storage().circumcenter_.value());
      }

      return storage().circumcenter_.value();
  }
  };


  template <typename GT,
            typename Cb = Cell<>,
            typename Memory_policy = Keep_hidden_points,
            typename C = std::list<typename GT::Weighted_point_3>>
  class Cell4Regular
    : public Cb
  {
  public:
    using Cb::Cb; // inherit constructors
    using Point_3 = typename GT::Point_3;
    using Point   = typename GT::Weighted_point_3;
    using TDS = typename Cb::Triangulation_data_structure;
    using Vertex_handle = typename TDS::Vertex_handle;
    using Cell_handle = typename TDS::Cell_handle;
    using Geom_traits = GT;

    using Point_container = C;
    using Point_iterator = typename Point_container::iterator;
    using  Point_const_iterator = typename Point_container::const_iterator;

    struct Storage : public Cb::Storage {
      Point_container hidden;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
      using Other = Cell4Regular<GT, Cb2, Memory_policy, C>;
    };

    auto&& storage() {
      return this->tds()->cell_storage()[this->index()];
    }

    auto&& storage() const { return this->tds()->cell_storage()[this->index()]; }

    void hide_point(const Point& p)
    {
      storage().hidden.push_back(p);
    }

    void unhide_point(const Point_iterator pit)
    {
      storage().hidden.erase(pit);
    }

    const C& hidden_points() const
    {
      return storage().hidden;
    }

    Point_iterator hidden_points_begin(){ return storage().hidden.begin();}
    Point_const_iterator hidden_points_begin() const { return storage().hidden.begin();}
    Point_iterator hidden_points_end(){ return storage().hidden.end();}
    Point_const_iterator hidden_points_end() const { return storage().hidden.end();}

  template<typename GT_>
  Point_3 weighted_circumcenter(const GT_& gt) const
  {
    return gt.construct_weighted_circumcenter_3_object()(this->vertex(0)->point(),
                                                         this->vertex(1)->point(),
                                                         this->vertex(2)->point(),
                                                         this->vertex(3)->point());
  }

  Point_3 weighted_circumcenter() const
  {
    return weighted_circumcenter(Geom_traits());
  }

  };
} // namespace CGAL

#include <CGAL/disable_warnings.h>

#endif // CGAL_INDEXED_STORAGE_H
