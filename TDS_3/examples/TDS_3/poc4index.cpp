#include <CGAL/TDS_3/internal/Dummy_tds_3.h>
#include <CGAL/TDS_3/internal/Triangulation_ds_circulators_3.h>
#include <CGAL/TDS_3/internal/Triangulation_ds_iterators_3.h>
#include <CGAL/Triangulation_utils_3.h>

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


using Point =  int;

namespace CGAL {

  template <typename TDS_3 = void>
  struct Vertex;

  template <typename TDS_3 = void>
  struct Cell;

  template <typename Vb = Vertex<>, typename Cb = Cell<>, class ConcurrencyTag = Sequential_tag>
  struct TDS;

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
      : tds_(nullptr), index_()
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
      assert(false);
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

    TDS_3* tds_;
    Cell_index index_;
  };

  template <typename TDS_3>
  std::ostream& operator<<(std::ostream& os, const Cell<TDS_3>& c)
  {
    os << "Cell " << c.index();
    return os;
  }

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
      Point point;
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

    const Point& point() const
    {
      return storage().point;
    }

    void set_point(const Point& p)
    {
      storage().point = p;
    }

    bool is_valid(bool  = false) const
    {
      assert(false); // This should be implemented
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
  };

  template <typename TDS_3>
  std::ostream& operator<<(std::ostream& os, const Vertex<TDS_3>& v)
  {
    os << "Vertex " << v.index();
    return os;
  }

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
    std::size_t ret = i;
    return ret;
  }



  template <typename Vb, typename Cb, typename ConcurrencyTag>
  struct TDS
  : public Triangulation_utils_3
  {

    using Self = TDS<Vb,Cb, ConcurrencyTag>;
    using Concurrency_tag  = ConcurrencyTag;

    /// The type used to represent an index.
    using size_type = std::uint32_t;

    // Tools to change the Vertex and Cell types of the TDS.
    template < typename Vb2 >
    struct Rebind_vertex {
      using Other = TDS<Vb2, Cb, ConcurrencyTag>;
    };

    template < typename Cb2 >
    struct Rebind_cell {
      using Other = TDS<Vb, Cb2, ConcurrencyTag>;
    };

    // Put this TDS inside the Vertex and Cell types.
    using Vertex = typename Vb::template Rebind_TDS<Self>::Other;
    using Cell = typename Cb::template Rebind_TDS<Self>::Other;

    template <typename T>
    class Handle {
      using Element = T;
      using Proxy = boost::stl_interfaces::proxy_arrow_result<Element>;
    public:
      using value_type = Element;
      using reference = Element;
      using pointer = Proxy;

      Handle() = default;

      Handle(Self* tds, size_type idx)
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

      bool operator==(const Handle& other) const {
        return tds() == other.tds() && index() == other.index();
      }

      bool operator!=(const Handle& other) const {
        return !(*this == other);
      }

      bool operator<(const Handle& other) const {
        return (tds() == other.tds()) ? (index() < other.index()) : (tds() < other.tds());
      }

      friend std::ostream& operator<<(std::ostream& os, const Handle& h)
      {
        return os << "#" << h.index();
      }

    private:
      Self* tds_ = nullptr;
      size_type idx_ = (std::numeric_limits<size_type>::max)();
    };

    using Cell_handle = Handle<Cell>;
    using Vertex_handle = Handle<Vertex>;

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

    void delete_vertex(Vertex_index v)
    {
      vremoved_[v] = true; ++removed_vertices_; garbage_ = true;
      vertex_storage_[v].cell = Cell_index(vertices_freelist_);
      vertices_freelist_ = (size_type)v;
    }

    void delete_cell(Cell_index c)
    {
      cremoved_[c] = true; ++removed_cells_; garbage_ = true;
      cell_storage_[c].ivertices[0] = Vertex_index(cells_freelist_);
      vertices_freelist_ = (size_type)c;
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
      return this == v->tds()  && v->idx() < num_vertices() && (! vremoved_[v.index()]);
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
      void resize_property_array() { m_->vprops_.resize_property_array(3); }
    };

    template<bool dummy>
    struct Property_selector<Cell_index, dummy> {
      Self * m_;
      Property_selector(Self* m) : m_(m) {}
      Properties::Property_container<Self,
                                     Cell_index>&
      operator()() { return m_->cprops_; }
      void resize_property_array() { m_->cprops_.resize_property_array(3); } // AF:  What is the 3 about?
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

    /// removes all property maps for all index types added by a call to `add_property_map()`.
    /// The memory allocated for those property maps is freed.
    void remove_all_property_maps()
    {
      remove_property_maps<Vertex_index>();
      remove_property_maps<Cell_index>();
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


    TDS()
      : dimension_(-2)
    {
      vertex_storage_  = add_property_map<Vertex_index, Vertex_storage>("v:storage").first;
      cell_storage_    = add_property_map<Cell_index, Cell_storage>("c:storage").first;
      cell_data_       = add_property_map<Cell_index, Cell_data>("c:data").first;
      vremoved_        = add_property_map<Vertex_index, bool>("v:removed", false).first;
      cremoved_        = add_property_map<Cell_index, bool>("c:removed", false).first;

      removed_vertices_ = removed_cells_ = 0;
      vertices_freelist_ = cells_freelist_  = (std::numeric_limits<size_type>::max)();
      garbage_ = false;
      recycle_ = true;
      anonymous_property_nb = 0;

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

    void just_incident_cells_3(Vertex_handle v,
                              std::vector<Cell_handle>& cells) const
    {
      CGAL_precondition(dimension() == 3);

      Cell_handle d = v->cell();
      cells.push_back(d);
      d->tds_data().mark_in_conflict();
      int head=0;
      int tail=1;
      do {
        Cell_handle c = cells[head];

        for (int i=0; i<4; ++i) {
          if (c->vertex(i) == v)
            continue;
          Cell_handle next = c->neighbor(i);
          if (! next->tds_data().is_clear())
            continue;
          cells.push_back(next);
          ++tail;
          next->tds_data().mark_in_conflict();
        }
        ++head;
      } while(head != tail);
    }

    Vertex_handle insert_first_finite_cell(Vertex_handle &v0,
                                           Vertex_handle &v1,
                                           Vertex_handle &v2,
                                           Vertex_handle &v3,
                                           Vertex_handle v_infinite = Vertex_handle())
    {
      CGAL_precondition(
                        (v_infinite == Vertex_handle() && dimension() == -2)
                        || (v_infinite != Vertex_handle() && dimension() == -1));

      if (v_infinite == Vertex_handle())
        v_infinite = create_vertex();

      set_dimension(3);

      v0 = create_vertex();
      v1 = create_vertex();
      v2 = create_vertex();
      v3 = create_vertex();

      Cell_handle c0123 = create_cell(v0,         v1,   v2,   v3);
      Cell_handle ci012 = create_cell(v_infinite, v0,   v1,   v2);
      Cell_handle ci103 = create_cell(v_infinite, v1,   v0,   v3);
      Cell_handle ci023 = create_cell(v_infinite, v0,   v2,   v3);
      Cell_handle ci132 = create_cell(v_infinite, v1,   v3,   v2);

      v_infinite->set_cell(ci012);
      v0->set_cell(c0123);
      v1->set_cell(c0123);
      v2->set_cell(c0123);
      v3->set_cell(c0123);

      set_adjacency(c0123, 0, ci132, 0);
      set_adjacency(c0123, 1, ci023, 0);
      set_adjacency(c0123, 2, ci103, 0);
      set_adjacency(c0123, 3, ci012, 0);

      set_adjacency(ci012, 3, ci103, 3);
      set_adjacency(ci012, 2, ci023, 3);
      set_adjacency(ci012, 1, ci132, 2);
      set_adjacency(ci103, 1, ci023, 2);
      set_adjacency(ci023, 1, ci132, 1);
      set_adjacency(ci103, 2, ci132, 3);

      return v_infinite;
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

      using value_type = typename Handle_::value_type;
    public:
      Index_iterator() : hnd_(), tds_(nullptr) {}
      Index_iterator(const Index& h, const Self* m)
        : hnd_(h), tds_(const_cast<Self*>(m)) { // AF: todo make const_cast safe
        if (tds_ && tds_->has_garbage()){
          while (tds_->has_valid_index(hnd_) && tds_->is_removed(hnd_)) ++hnd_;
        }
      }

      auto handle() const
      {
        static_assert(std::is_base_of_v<Facade, Index_iterator<Handle_>>);

        CGAL_assertion(tds_ != nullptr);
        CGAL_assertion(tds_->has_valid_index(hnd_));
        return Handle_(tds_, hnd_.id());
      }

      operator Handle_() const { return handle(); }

      value_type operator*() const { return value_type{tds_, hnd_}; }

      typename Handle_::pointer operator->() const
      {
        return typename Handle_::pointer{value_type{tds_, hnd_}};
      }

      using Facade::operator++;
      Index_iterator& operator++()
      {
        ++hnd_;
        CGAL_assertion(tds_ != nullptr);

        if(tds_->has_garbage())
          while ( tds_->has_valid_index(hnd_) && tds_->is_removed(hnd_)) ++hnd_;
        return *this;
      }

      using Facade::operator--;
      Index_iterator& operator--()
      {
        --hnd_;
        CGAL_assertion(tds_ != nullptr);
        if(tds_->has_garbage())
          while ( tds_->has_valid_index(hnd_) && tds_->is_removed(hnd_)) --hnd_;
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
          hnd_ += n;
        return *this;
      }

      std::ptrdiff_t operator-(const Index_iterator& other) const
      {
        if (tds_->has_garbage())
          {
            bool forward = (other.hnd_ > hnd_);

            std::ptrdiff_t out = 0;
            Index_iterator it = *this;
            while (!it.equal(other))
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
        return std::ptrdiff_t(other.hnd_) - std::ptrdiff_t(this->hnd_);
      }

      bool operator==(const Index_iterator& other) const
      {
        return this->hnd_ == other.hnd_;
      }

    private:
      Index hnd_;
      Self* tds_;
    };

    // AF:  The value type of the iterators should be the same as the value type of the handle. This is not the case now.
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

namespace CGAL {
  template <typename Cb = Cell<>>
  class Cell4Regular
    : public Cb
  {
  public:
    using Cb::Cb; // inherit constructors

    using TDS = typename Cb::Triangulation_data_structure;
    using Vertex_handle = typename TDS::Vertex_handle;
    using Cell_handle = typename TDS::Cell_handle;

    struct Storage : public Cb::Storage {
      std::list<Point> hidden_points;
      Vertex_handle my_other_vertex;
      Cell_handle my_other_cell;
    };

    template < typename TDS2 >
    struct Rebind_TDS {
      using Cb2 = typename Cb::template Rebind_TDS<TDS2>::Other;
      using Other = Cell4Regular<Cb2>;
    };

    auto&& storage() {
      return this->tds()->cell_storage()[this->index()];
    }

    auto&& storage() const { return this->tds()->cell_storage()[this->index()]; }

    void hide_point(const Point& p)
    {
      storage().hidden_points.push_back(p);
    }
  };
} // namespace CGAL

int main() {
  using Vb = CGAL::Vertex<>;
  using Cb = CGAL::Cell4Regular<>;
  using TDS = CGAL::TDS<Vb, Cb>;
  TDS tds;
  using Vertex_handle = TDS::Vertex_handle;
  using Vertex_iterator = TDS::Vertex_iterator;
  using Cell_handle = TDS::Cell_handle;
  using Cell_circulator = TDS::Cell_circulator;
  using Facet_circulator = TDS::Facet_circulator;
  using Edge = TDS::Edge;

  Vertex_handle inf, v0, v1, v2, v3;
  inf = tds.insert_first_finite_cell(v0, v1, v2, v3);

  for(const auto& v : tds.vertices()) {
    std::cout << v << std::endl;
  }

  for(Vertex_iterator vit = tds.vertices_begin(); vit != tds.vertices_end(); ++vit) {
    std::cout << vit->point() << std::endl;
  }

  for(auto c : tds.cells()) {
    std::cout << c.index() << std::endl;
    c.hide_point(Point{});
  }

  inf->cell()->storage().my_other_vertex = v0;
  inf->cell()->storage().my_other_cell = v1->cell();

  std::vector<Cell_handle> incident_cells;
  tds.just_incident_cells_3(inf, incident_cells);
  assert(incident_cells.size() == 4);

  for(auto f : tds.facets()) {
    std::cout << f.first << " " << f.second << std::endl;
  }

  for(auto e : tds.edges()) {
    std::cout << e.first << " " << e.second << " " << e.third << std::endl;
  }

  Cell_handle c = inf->cell();
  Edge e(c, 0, 1);

  std::cout << "Incident cells to edge:\n";
  Cell_circulator cc = tds.incident_cells(e), cdone(cc);
  do {
    Cell_handle ch = cc;
    std::cout << ch << " ";
    ++cc;
  } while (cc != cdone);
  std::cout << std::endl;

  std::cout << "Incident facets to edge:\n";
  Facet_circulator fc = tds.incident_facets(e), fdone(fc);
  do {
    Cell_handle ch = fc->first;
    std::cout << ch << " ";
    ++fc;
  } while (fc != fdone);
  std::cout << std::endl;
  return 0;
}
