
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Path_on_surface` represents a walk in a mesh which is either a 2D `CombinatorialMap` or a model of a `FaceGraph`. Each object of this class is constructed from an external mesh on which the path should lie. A path is represented as a sequence of darts or halfedges, each one representing an oriented edge in the path. The class `Path_on_surface` behaves as a container for this sequence of darts/halfedges. Elements are added in the path one at a time to the path thanks to the `push_back()` method.
    
    \tparam Mesh a model of `CombinatorialMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Path_on_surface
  {
  public:
    /*!
      %halfedge_descriptor type. A handle to Dart for combinatorial maps, or a halfedge descriptor for models of FaceGraph.
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

    /// clears this path.
    void clear();
    
    /// returns `true` iff `hd` can be added at the end of this path.
    bool can_be_pushed(halfedge_descriptor hd) const;
  
    /// adds `hd` at the end of this path.
    /// @pre `can_be_pushed(hd)`
    void push_back(halfedge_descriptor hd);

    /// adds the dart/halfedge with index `i` at the end of this path.
    /// @pre `can_be_pushed_by_index(i)`
    void push_back_by_index(std::size_t i);

    /// returns `true` iff the dart/halfedge with index `i` can be added at the end of this path.
    bool can_be_pushed_by_index(std::size_t i) const;

    /// Add the dart/halfedge obtained by turning `nb` times around the target vertex of the last dart/halfedge in this path, in the positive circular order.
    /// @pre !`is_empty()`
    void extend_positive_turn(std::size_t nb); 

    /// Add the dart/halfedge obtained by turning `nb` times around the target vertex of the last dart/halfedge in this path, in the negative circular order.
    /// @pre !`is_empty()`
    void extend_negative_turn(std::size_t nb); 

    /// concatenates `other` to this path.
    /// @pre the last vertex of this path should coincide with the first vertex of `other`.
    Self& operator+=(const Self& other);

    /// reverses this path (i.e. negate its orientation).
    void reverse();

    /// creates a random open path with `length` darts/halfedges.
    void generate_random_path(std::size_t length, Random& random=get_default_random());

    /// creates a random closed path with at least `length` darts/halfedges.
    void generate_random_closed_path(std::size_t length, Random& random=get_default_random());
    
   };
}
