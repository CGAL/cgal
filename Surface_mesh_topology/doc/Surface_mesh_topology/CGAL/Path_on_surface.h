
namespace CGAL {
namespace Surface_mesh_topology {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses

    The class `Path_on_surface` represents a walk in a mesh which is either a model of `CombinatorialMap`, a model of `GeneralizedMap` or a model of a `FaceGraph`. Each object of this class is constructed from an external mesh on which the path should lie. A path is represented as a sequence of darts or halfedges, each one representing an oriented edge in the path. The class `Path_on_surface` behaves as a container for this sequence of darts/halfedges. Elements are added in the path one at a time to the path thanks to the `push_back()` method.

    \tparam Mesh a model of `CombinatorialMap`, `GeneralizedMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Path_on_surface
  {
  public:
    /*!
      A handle to `Dart` for combinatorial/generalized maps, or a halfedge descriptor for models of the `FaceGraph` concept.
    */
    typedef unspecified_type halfedge_descriptor;

    /// creates an empty path object which lies on `amesh`.
    Path_on_surface(const Mesh& amesh);

    /// returns `true` iff this path is empty
    bool is_empty() const;

    /// returns `true` iff this path is closed.
    bool is_closed() const;

    /// returns `true` iff this path does not pass twice through a same edge or a same vertex.
    bool is_simple() const;

    /// returns the length of the path, i.e. its number of edges.
    std::size_t length() const;

    /// returns the ith dart of the path.
    /// @pre i<`length()`.
    halfedge_descriptor operator[] (std::size_t i) const;

    /// clears this path.
    void clear();

    /// returns `true` iff `hd` can be added at the end of this path. If `flip` is true, the direction of `hd` is reversed before checking
    bool can_be_pushed(halfedge_descriptor hd, bool flip=false) const;

    /// adds `hd` at the end of this path. If `flip` is true, the opposite of `hd` is considered.
    /// @pre `can_be_pushed(hd)`
    void push_back(halfedge_descriptor hd, bool flip=false);

    /// returns `true` iff the dart/halfedge with index `i` can be added at the end of this path.
    /// If Mesh is a `Polyhedron_3`, takes time proportional to the number of darts/halfedges.
    bool can_be_pushed_by_index(std::size_t i) const;

    /// adds the dart/halfedge with index `i` at the end of this path.
    /// If Mesh is a `Polyhedron_3`, takes time proportional to the number of halfedges.
    /// @pre `can_be_pushed_by_index(i)`
    void push_back_by_index(std::size_t i);

    /// adds successively all dart/halfedges in `l` (a sequence of indices), at the end of the path.
    /// If Mesh is a `Polyhedron_3`, takes time proportional to the number of halfedges.
    /// For each index `i`, `can_be_pushed_by_index(i)` should be true.
    void push_back_by_index(std::initializer_list<std::size_t> l);

    /// adds the dart/halfedge obtained by turning `nb` times around the target vertex of the last dart/halfedge in this path, in the positive circular order. To extend with a positive 1 turn thus amounts to extend with the `next()` pointer. (A zero turn corresponds to the `opposite()` pointer.)
    /// @pre !`is_empty()`
    void extend_positive_turn(std::size_t nb);

    /// adds the dart/halfedge obtained by turning `nb` times around the target vertex of the last dart/halfedge in this path, in the negative circular order. To extend with a negative 1 turn thus amounts to extend with the composite pointer: `opposite(prev(opposite()))`.
    /// @pre !`is_empty()`
    void extend_negative_turn(std::size_t nb);

    /// returns `true` iff the dart/halfedge with label `l` can be added at the end of this path.
    /// @pre Mesh must be a model of `PolygonalSchema` concept.
    bool can_be_pushed_by_label(const std::string& l) const;

    /// adds successively all darts/halfedges in `s` (a sequence of labels separated by spaces) at the end of this path. For each label, l, `can_be_pushed_by_label(l)` should be true.
    /// @pre Mesh must be a model of `PolygonalSchema` concept.
    void push_back_by_label(const std::string& s);

    /// concatenates `other` to this path.
    /// @pre the last vertex of this path should coincide with the first vertex of `other`.
    Self& operator+=(const Self& other);

    /// reverses this path (i.e. negates its orientation).
    void reverse();

    /// creates a random open path with `length` darts/halfedges.
    void generate_random_path(std::size_t length, Random& random=get_default_random());

    /// creates a random closed path with at least `length` darts/halfedges.
    void generate_random_closed_path(std::size_t length, Random& random=get_default_random());
   };
}
}
