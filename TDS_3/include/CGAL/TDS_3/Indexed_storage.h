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
#include <CGAL/Iterator_range.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/utility.h>

#include <boost/stl_interfaces/iterator_interface.hpp>

#include <array>
#include <iostream>
#include <list>
#include <type_traits>
#include <vector>
#include <optional>

namespace CGAL {

  template <typename TDS_3 = void>
  struct Vertex;

  template <typename TDS_3 = void>
  struct Cell;

  template <typename Vb = Vertex<>, typename Cb = Cell<>, class ConcurrencyTag = Sequential_tag>
  struct Indexed_storage;

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

  template <typename TDS_3>
  struct Cell {

    using Triangulation_data_structure = TDS_3;
    using Vertex_handle = typename TDS_3::Vertex_handle;
    using Cell_handle   = typename TDS_3::Cell_handle;
    using Vertex_index  = typename TDS_3::Vertex_index;
    using Cell_index    = typename TDS_3::Cell_index;
    using Index         = Cell_index;
    using TDS_data      = typename TDS_3::Cell_data;

    struct Storage {
      std::array<Vertex_index,4> ivertices;
      std::array<Cell_index,4>   ineighbors;
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



  template <typename TDS_3>
  struct Vertex{

    using Triangulation_data_structure = TDS_3;
    using Cell_handle   = typename TDS_3::Cell_handle;
    using Cell_index    = typename TDS_3::Cell_index;
    using Vertex_index  = typename TDS_3::Vertex_index;
    using Index         = Vertex_index;
    using Vertex_handle = typename TDS_3::Vertex_handle;

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
  template<typename T>
  class Index
  {
  public:
    using size_type = std::uint32_t;

    /// Constructor. %Default construction creates an invalid index.
    /// We write -1, which is <a href="https://en.cppreference.com/w/cpp/types/numeric_limits">
    /// <tt>(std::numeric_limits<size_type>::max)()</tt></a>
    /// as `size_type` is an unsigned type.
    explicit Index(size_type _idx=(std::numeric_limits<size_type>::max)()) : idx_(_idx) {}

    /// Get the underlying index of this index
    operator size_type() const { return idx_; }

    /// reset index to be invalid (index=(std::numeric_limits<size_type>::max)())
    void reset() { idx_=(std::numeric_limits<size_type>::max)(); }

    /// return whether the index is valid, i.e., the index is not equal to `%std::numeric_limits<size_type>::max()`.
    bool is_valid() const {
      size_type inf = (std::numeric_limits<size_type>::max)();
      return idx_ != inf;
    }

    // Compatibility with OpenMesh handle
    size_type idx() const {
      return idx_;
    }
    // For convenience
    size_type id() const {
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

  protected:
    size_type idx_;
  };

  template <class T>
  std::size_t hash_value(const Index<T>&  i)
  {
    return i.id();
  }

  template <typename T>
  class TDS_handle {
    using Element = T;
    using TDS = typename Element::Triangulation_data_structure;
    using size_type = typename TDS::size_type;
    using Proxy = boost::stl_interfaces::proxy_arrow_result<Element>;
  public:
    using value_type = Element;
    using reference = Element;
    using pointer = Proxy;

    TDS_handle() = default;

    TDS_handle(TDS* tds, size_type idx)
      : tds_(tds), idx_(idx) {}

    using Index = typename Element::Index;

    Element operator*() const {
      return Element(tds_, Index(idx_));
    }

    Proxy operator->() const {
      return Proxy{this->operator*()};
    }

    Index index() const {
      return Index{idx_};
    }

    auto tds() const {
      return tds_;
    }

    bool operator==(const TDS_handle& other) const {
      CGAL_assertion(tds() == nullptr || other.tds() == nullptr ||
                     tds() == other.tds());
      return index() == other.index();
    }

    bool operator!=(const TDS_handle& other) const {
      return !(*this == other);
    }

    bool operator<(const TDS_handle& other) const {
      return (tds() == other.tds()) ? (index() < other.index()) : (tds() < other.tds());
    }

    bool operator==( std::nullptr_t ) const {
      return tds_ == nullptr && idx_ == (std::numeric_limits<size_type>::max)();
    }

    bool operator!=( std::nullptr_t ) const {
      return !(*this == nullptr);
    }

    friend std::ostream& operator<<(std::ostream& os, const TDS_handle& h)
    {
      return os << "#" << h.index();
    }

    friend std::size_t hash_value(const TDS_handle& h) {
      return static_cast<std::size_t>(h.index().id());
    }

  private:
    TDS* tds_ = nullptr;
    size_type idx_ = (std::numeric_limits<size_type>::max)();
  };

} // end namespace CGAL

namespace std {
  template <typename T>
  struct hash<CGAL::TDS_handle<T>> {
    using Handle = CGAL::TDS_handle<T>;
    std::size_t operator()(const Handle& h) const {
      return hash_value(h);
    }
  };
}

namespace CGAL {

  template <typename Vb, typename Cb, typename ConcurrencyTag>
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

    using Cell_handle = TDS_handle<Cell>;
    using Vertex_handle = TDS_handle<Vertex>;

    using Facet = std::pair<Cell_handle, int>;
    using Edge = Triple<Cell_handle, int, int>;

    // AF: factorize as also in Surface_mesh and Point_set_3
    template <class I, class T>
    struct Property_map : Properties::Property_map_base<I, T, Property_map<I, T> >
    {
      using Base = Properties::Property_map_base<I, T, Property_map<I, T> >;
      using reference = typename Base::reference;

      Property_map() = default;
      Property_map(const Base& pm): Base(pm) {}
    };


    // AF: factorize as also in Surface_mesh and Point_set_3
    template <typename Key, typename T>
    struct Get_property_map {
      using type = Property_map<Key, T>;
    };


    class Vertex_index
      : public Index<Vertex_index>
    {
    public:

      Vertex_index() : Index<Vertex_index>((std::numeric_limits<size_type>::max)()) {}

      explicit Vertex_index(size_type _idx) : Index<Vertex_index>(_idx) {}

      template<class T> bool operator==(const T&) const = delete;
      template<class T> bool operator!=(const T&) const = delete;
      template<class T> bool operator<(const T&) const = delete;

      /// are two indices equal?
      bool operator==(const Vertex_index& _rhs) const {
        return this->idx_ == _rhs.idx_;
      }

      /// are two indices different?
      bool operator!=(const Vertex_index& _rhs) const {
        return this->idx_ != _rhs.idx_;
      }

      /// Comparison by index.
      bool operator<(const Vertex_index& _rhs) const {
        return this->idx_ < _rhs.idx_;
      }


      friend std::ostream& operator<<(std::ostream& os, Vertex_index const& v)
      {
        return (os << 'v' << (size_type)v );
      }
    };

    class Cell_index
      : public Index<Cell_index>
    {
    public:

      Cell_index() : Index<Cell_index>((std::numeric_limits<size_type>::max)()) {}

      explicit Cell_index(size_type _idx) : Index<Cell_index>(_idx) {}

      template<class T> bool operator==(const T&) const = delete;
      template<class T> bool operator!=(const T&) const = delete;
      template<class T> bool operator<(const T&) const = delete;

      /// are two indices equal?
      bool operator==(const Cell_index& _rhs) const {
        return this->idx_ == _rhs.idx_;
      }

      /// are two indices different?
      bool operator!=(const Cell_index& _rhs) const {
        return this->idx_ != _rhs.idx_;
      }

      /// Comparison by index.
      bool operator<(const Cell_index& _rhs) const {
        return this->idx_ < _rhs.idx_;
      }


      friend std::ostream& operator<<(std::ostream& os, Cell_index const& v)
      {
        return (os << 'c' << (size_type)v );
      }
    };

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
    using Vertex_storage_property_map = Property_map<Vertex_index, Vertex_storage>;
    using Cell_storage_property_map = Property_map<Cell_index, Cell_storage>;



    Cell_storage_property_map& cell_storage() { return cell_storage_; }
    const Cell_storage_property_map& cell_storage() const { return cell_storage_; }

    Vertex_storage_property_map& vertex_storage() { return vertex_storage_; }
    const Vertex_storage_property_map& vertex_storage() const { return vertex_storage_; }

    auto cell_tds_data_pmap() { return cell_data_; }
    auto cell_tds_data_pmap() const { return cell_data_; }

    size_type num_vertices() const { return (size_type) vprops_.size(); }
    size_type num_cells() const { return (size_type) cprops_.size(); }


    int dimension() const { return dimension_; }

    void set_dimension(int n) { dimension_ = n; }

    Vertex_handle create_vertex()
    {
      size_type inf = (std::numeric_limits<size_type>::max)();
      if(recycle_ && (vertices_freelist_ != inf)){
        size_type idx = vertices_freelist_;
        vertices_freelist_ = (size_type)vertex_storage_[Vertex_index(vertices_freelist_)].icell;
        --removed_vertices_;
        vremoved_[Vertex_index(idx)] = false;
        vprops_.reset(Vertex_index(idx));
        return Vertex_handle(this, Vertex_index(idx));
      } else {
        vprops_.push_back();
        return Vertex_handle(this, Vertex_index(num_vertices()-1));
      }
    }

    Cell_handle create_cell()
    {
      size_type inf = (std::numeric_limits<size_type>::max)();
      if(recycle_ && (cells_freelist_ != inf)){
        size_type idx = cells_freelist_;
        cells_freelist_ = (size_type)cell_storage_[Cell_index(cells_freelist_)].ivertices[0];
        --removed_cells_;
        cremoved_[Cell_index(idx)] = false;
        cprops_.reset(Cell_index(idx));
        return Cell_handle(this, Cell_index(idx));
      } else {
        cprops_.push_back();
        return Cell_handle(this, Cell_index(num_cells()-1));
      }
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
      Vertex_index v = vh->index();
      vremoved_[v] = true; ++removed_vertices_; garbage_ = true;
      vertex_storage_[v].icell = Cell_index{vertices_freelist_};
      vertices_freelist_ = static_cast<size_type>(v);
    }

    void delete_cell(Cell_handle ch)
    {
      Cell_index c = ch->index();
      cremoved_[c] = true; ++removed_cells_; garbage_ = true;
      cell_storage_[c].ivertices[0] = Vertex_index{cells_freelist_};
      cells_freelist_ = static_cast<size_type>(c);
    }

    void reserve(size_type n_vertices, size_type n_cells)
    {
      vprops_.reserve(n_vertices);
      cprops_.reserve(n_cells);
    }

    void clear()
    {
      clear_without_removing_property_maps();
      remove_all_property_maps();
      allocate_tds_properties();
    }

    void clear_without_removing_property_maps()
    {
      vprops_.resize(0);
      cprops_.resize(0);

      vprops_.shrink_to_fit();
      cprops_.shrink_to_fit();

      removed_vertices_ = removed_cells_ = 0;
      vertices_freelist_ = cells_freelist_ = (std::numeric_limits<size_type>::max)();
      garbage_ = false;
      recycle_ = true;
      anonymous_property_nb = 0;
      dimension_ = -2;
    }

    void remove_all_property_maps()
    {
      remove_property_maps<Vertex_index>();
      remove_property_maps<Cell_index>();
    }

    size_type number_of_removed_vertices() const
    {
      return removed_vertices_;
    }

    size_type number_of_removed_cells() const
    {
      return removed_cells_;
    }

    size_type number_of_vertices() const
    {
      return num_vertices() - number_of_removed_vertices();
    }

    size_type number_of_cells() const
    {
      return num_cells() - number_of_removed_cells();
    }


    bool is_vertex(Vertex_handle v) const
    {
      return this == v->tds()  && v->index().id() < num_vertices() && (! vremoved_[v.index()]);
    }

    bool is_valid_cell_handle(Cell_handle c) const
    {
      return this == c->tds()  && c->index().id() < num_cells() && (! cremoved_[c.index()]);
    }

    bool is_cell( Cell_handle c ) const
      // returns false when dimension <3
    {
      if (dimension() < 3)
        return false;

      return is_valid_cell_handle(c);
    }

    //--------------------------------------------------- property handling

    // Property_selector maps an index type to a property_container, the
    // dummy is necessary to make it a partial specialization (full
    // specializations are only allowed at namespace scope).
    template<typename, bool = true>
    struct Property_selector {};

    template<bool dummy>
    struct Property_selector<Vertex_index, dummy> {
      Self * m_;
      Property_selector(Self* m) : m_(m) {}
      Properties::Property_container<Self,
                                     Vertex_index>&
      operator()() { return m_->vprops_; }
      void resize_property_array() { m_->vprops_.resize_property_array(2); }
    };

    template<bool dummy>
    struct Property_selector<Cell_index, dummy> {
      Self * m_;
      Property_selector(Self* m) : m_(m) {}
      Properties::Property_container<Self,
                                     Cell_index>&
      operator()() { return m_->cprops_; }
      void resize_property_array() { m_->cprops_.resize_property_array(3); }
    };

    /// adds a property map named `name` with value type `T` and default `t`
    /// for index type `I`. Returns the property map together with a Boolean
    /// that is `true` if a new map was created. In case it already exists
    /// the existing map together with `false` is returned.


    template<class I, class T>
    std::pair<Property_map<I, T>, bool>
    add_property_map(std::string name=std::string(), const T t=T()) {
      if(name.empty()){
        std::ostringstream oss;
        oss << "anonymous-property-" << anonymous_property_nb++;
        name = std::string(oss.str());
      }
      return Property_selector<I>(this)().template add<T>(name, t);
    }

    /// returns an optional property map named `name` with key type `I` and value type `T`.
    template <class I, class T>
    std::optional<Property_map<I, T>> property_map(const std::string& name) const
    {
      return Property_selector<I>(const_cast<Self*>(this))().template get<T>(name);
    }


    /// removes property map `p`. The memory allocated for that property map is freed.
    template<class I, class T>
    void remove_property_map(Property_map<I, T>& p)
    {
      (Property_selector<I>(this)()).template remove<T>(p);
    }

    /// removes all property maps for index type `I` added by a call to `add_property_map<I>()`.
    /// The memory allocated for those property maps is freed.
    template<class I>
    void remove_property_maps()
    {
      Property_selector<I>(this).resize_property_array();
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
      return Property_selector<I>(this)().get_type(name);
    }
    /// @endcond

    /// returns a vector with all strings that describe properties with the key type `I`.
    /// @tparam I The key type of the properties.
    template<class I>
    std::vector<std::string> properties() const
    {
      return Property_selector<I>(const_cast<Self*>(this))().properties();
    }

    void allocate_tds_properties() {
      vertex_storage_ = add_property_map<Vertex_index, Vertex_storage>("v:storage").first;
      cell_storage_   = add_property_map<Cell_index, Cell_storage>("c:storage").first;
      cell_data_      = add_property_map<Cell_index, Cell_data>("c:data").first;
      vremoved_       = add_property_map<Vertex_index, bool>("v:removed", false).first;
      cremoved_       = add_property_map<Cell_index, bool>("c:removed", false).first;
    }


    Indexed_storage()
      : dimension_(-2)
    {
      allocate_tds_properties();
      removed_vertices_ = removed_cells_ = 0;
      vertices_freelist_ = cells_freelist_  = (std::numeric_limits<size_type>::max)();
      garbage_ = false;
      recycle_ = true;
      anonymous_property_nb = 0;
    }

    Indexed_storage(Indexed_storage&& is) noexcept
      : vprops_(std::move(is.vprops_)),
        cprops_(std::move(is.cprops_)),
        vertex_storage_(std::move(is.vertex_storage_)),
        cell_storage_(std::move(is.cell_storage_)),
        cell_data_(std::move(is.cell_data_)),
        vremoved_(std::move(is.vremoved_)),
        cremoved_(std::move(is.cremoved_)),
        removed_vertices_(std::exchange(is.removed_vertices_, 0)),
        removed_cells_(std::exchange(is.removed_cells_, 0)),
        vertices_freelist_(std::exchange(is.vertices_freelist_, (std::numeric_limits<size_type>::max)())),
        cells_freelist_(std::exchange(is.cells_freelist_, (std::numeric_limits<size_type>::max)())),
        garbage_(std::exchange(is.garbage_, false)),
        recycle_(std::exchange(is.recycle_, true)),
        anonymous_property_nb(std::exchange(is.anonymous_property_nb, 0)),
        dimension_(is.dimension_)
    {
      is.dimension_ = -2;
    }

    Indexed_storage& operator=(const Indexed_storage& rhs)
    {
      if (this != &rhs)
      {
        // clear properties
        vprops_.clear();
        cprops_.clear();

        // allocate standard properties
        allocate_tds_properties();

        // copy properties from other triangulation
        vertex_storage_.array() = rhs.vertex_storage_.array();
        cell_storage_.array()   = rhs.cell_storage_.array();
        cell_data_.array()      = rhs.cell_data_.array();

        vremoved_.array()  = rhs.vremoved_.array();
        cremoved_.array()  = rhs.cremoved_.array();

        // resize (needed by property containers)
        vprops_.resize(rhs.num_vertices());
        cprops_.resize(rhs.num_cells());


        // how many elements are removed?
        removed_vertices_     = rhs.removed_vertices_;
        removed_cells_        = rhs.removed_cells_;
        vertices_freelist_    = rhs.vertices_freelist_;
        cells_freelist_       = rhs.cells_freelist_;
        garbage_              = rhs.garbage_;
        recycle_              = rhs.recycle_;
        anonymous_property_nb = rhs.anonymous_property_nb;
      }
      return *this;
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


    /// move assignment
    Indexed_storage& operator=(Indexed_storage&& is) noexcept
    {
      if (this != &is) {
        dimension_ = is.dimension_;
        vertex_storage_ = std::move(is.vertex_storage_);
        cell_storage_ = std::move(is.cell_storage_);
        cell_data_ = std::move(is.cell_data_);
        vprops_ = std::move(is.vprops_);
        cprops_ = std::move(is.cprops_);
        vremoved_ = std::move(is.vremoved_);
        cremoved_ = std::move(is.cremoved_);
        removed_vertices_ = std::exchange(is.removed_vertices_, 0);
        removed_cells_ = std::exchange(is.removed_cells_, 0);
        vertices_freelist_ = std::exchange(is.vertices_freelist_, (std::numeric_limits<size_type>::max)());
        cells_freelist_ = std::exchange(is.cells_freelist_, (std::numeric_limits<size_type>::max)());
        garbage_ = std::exchange(is.garbage_, false);
        recycle_ = std::exchange(is.recycle_, true);
        anonymous_property_nb = std::exchange(is.anonymous_property_nb, 0);
        is.dimension_ = -2;
      }
      return *this;
    }

    void swap(Indexed_storage& other)
    {
      std::swap(*this, other);
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

    bool has_garbage() const { return garbage_; }

    /// returns whether the index of vertex `v` is valid, that is within the current array bounds.
    bool has_valid_index(Vertex_index v) const
    {
      return ((size_type)v < num_vertices());
    }
    bool has_valid_index(Cell_index c) const
    {
      return ((size_type)c < num_cells());
    }

    /// returns whether vertex `v` is marked removed.
    /// \sa `collect_garbage()`
    bool is_removed(Vertex_index v) const
    {
      return vremoved_[v];
    }

    bool is_removed(Cell_index c) const
    {
      return cremoved_[c];
    }
    //------------------------------------------------------ iterator types
    template<typename Handle_>
    class Index_iterator // make it derive from Handle_
      : public boost::stl_interfaces::v1::proxy_iterator_interface<
                   Index_iterator<Handle_>,
                   std::random_access_iterator_tag,
                   typename Handle_::value_type
                   >
    {
      using Index = typename Handle_::Index;

      using Facade = boost::stl_interfaces::v1::proxy_iterator_interface<
                         Index_iterator<Handle_>,
                         std::random_access_iterator_tag,
                         typename Handle_::value_type
                         >;

    public:
      using value_type = typename Handle_::value_type;
      using reference = value_type;

      Index_iterator() : idx_(), tds_(nullptr) {}
      Index_iterator(const Index& h, const Self* m)
        : idx_(h), tds_(const_cast<Self*>(m)) { // AF: todo make const_cast safe
        if (tds_ && tds_->has_garbage()){
          while (tds_->has_valid_index(idx_) && tds_->is_removed(idx_)) ++idx_;
        }
      }

      auto handle() const
      {
        static_assert(std::is_base_of_v<Facade, Index_iterator<Handle_>>);

        CGAL_assertion(tds_ != nullptr);
        CGAL_assertion(tds_->has_valid_index(idx_));
        return Handle_(tds_, idx_.id());
      }

      operator Handle_() const { return handle(); }

      reference operator*() const { return value_type{tds_, idx_}; }

      typename Handle_::pointer operator->() const
      {
        return typename Handle_::pointer{value_type{tds_, idx_}};
      }

      using Facade::operator++;
      Index_iterator& operator++()
      {
        ++idx_;
        CGAL_assertion(tds_ != nullptr);

        if(tds_->has_garbage())
          while ( tds_->has_valid_index(idx_) && tds_->is_removed(idx_)) ++idx_;
        return *this;
      }

      using Facade::operator--;
      Index_iterator& operator--()
      {
        --idx_;
        CGAL_assertion(tds_ != nullptr);
        if(tds_->has_garbage())
          while ( tds_->has_valid_index(idx_) && tds_->is_removed(idx_)) --idx_;
        return *this;
      }

      Index_iterator& operator+=(std::ptrdiff_t n)
      {
        CGAL_assertion(tds_ != nullptr);

        if (tds_->has_garbage())
          {
            if (n > 0)
              for (std::ptrdiff_t i = 0; i < n; ++ i)
                this->operator++();
            else
              for (std::ptrdiff_t i = 0; i < -n; ++ i)
                this->operator--();
          }
        else
          idx_ += n;
        return *this;
      }

      std::ptrdiff_t operator-(const Index_iterator& other) const
      {
        if (tds_->has_garbage())
          {
            bool forward = (other.idx_ > idx_);

            std::ptrdiff_t out = 0;
            Index_iterator it = *this;
            while (!(it == other))
              {
                if (forward)
                  {
                    ++ it;
                    ++ out;
                  }
                else
                  {
                    -- it;
                    -- out;
                  }
              }
            return out;
          }

        // else
        return std::ptrdiff_t(other.idx_) - std::ptrdiff_t(this->idx_);
      }

      bool operator==(const Index_iterator& other) const
      {
        return this->idx_ == other.idx_;
      }

    private:
      Index idx_;
      Self* tds_;
    };

    using Vertex_iterator = Index_iterator<Vertex_handle>;
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


    using Cell_iterator = Index_iterator<Cell_handle>;
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

    Properties::Property_container<Self, Vertex_index> vprops_;
    Properties::Property_container<Self, Cell_index>   cprops_;

    Vertex_storage_property_map vertex_storage_;
    Cell_storage_property_map cell_storage_;

    Property_map<Cell_index, Cell_data>                cell_data_;

    Property_map<Vertex_index, bool>  vremoved_;
    Property_map<Cell_index, bool>    cremoved_;


    size_type removed_vertices_;
    size_type removed_cells_;

    size_type vertices_freelist_;
    size_type cells_freelist_;
    bool garbage_;
    bool recycle_;

    size_type anonymous_property_nb;

    // in dimension i, number of vertices >= i+2
    // ( the boundary of a simplex in dimension i+1 has i+2 vertices )
    int dimension_;

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
    {}

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
      if (! storage().circumcenter_){
       storage().circumcenter_ = std::make_optional<Point>(p);
      }
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
