
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Path_on_surface` represents a walk in a mesh which is either a 2D `CombinatorialMap` or a model of a `FaceGraph`. Each object in this class is constructed from an external mesh on which the path should lie. A path is represented as a sequence of darts or half-edges, each one representing an oriented edge in the path. The class `Path_on_surface` behaves as a container for this sequence of darts/half-edges. Elements are added in the path one at a time to the path thanks to the `push_back()` method.
    
    \tparam Mesh a model of `CombinatorialMap` or of `FaceGraph`
  */
  template<typename Mesh>
  class Path_on_surface
  {
  public:
    /*!
      %Dart const handle type. An handle to Dart for combinatorial maps, or an half-edge descriptor for models of FaceGraph.
    */
    typedef unspecified_type Dart_const_handle;

    /// Constructor. Creates an empty path object which lies on amesh.
    Path_on_surface(const Mesh& amesh);

    /// @return `true` iff the path is empty
    bool is_empty() const;

    /// @return `true` iff the path is closed.
    bool is_closed() const;

    /// @return `true` iff the path does not pass twice through a same edge
    ///              or a same vertex.
    bool is_simple() const;

  /// clear the path.
    void clear();
    
    /// @return `true` iff dh can be added at the end of the path.
    bool can_be_pushed(Dart_const_handle dh) const;
  
    /// Add the given dart/half-edge at the end of this path.
    /// @pre can_be_pushed(dh)
    void push_back(Dart_const_handle dh);

    /// Add the dart/half-edge with given index i at the end of this path.
    /// @pre can_be_pushed_by_index(i)
    void push_back_by_index(std::size_t i);

    /// @return `true` iff the dart/half-edge with index i can be added at the end of the path.
    bool can_be_pushed_by_index(std::size_t i) const;

    /// Add the dart/half-edge obtained by turning nb times around the target vertex of the last dart/half-edge in this path, in the positive circular order.
    /// @pre !is_empty()
    void extend_positive_turn(std::size_t nb); 

    /// Add the dart/half-edge obtained by turning nb times around the target vertex of the last dart/half-edge in this path, in the negative circular order.
    /// @pre !is_empty()
    void extend_negative_turn(std::size_t nb); 

    /// Concatenation operator. Concatenates other to this path.
    /// @pre the last vertex of this path should coincide with the first vertex of other.
    Self& operator+=(const Self& other);

    /// Reverse the path (i.e. negate its orientation).
    void reverse();

    /// Creates a random open path with lenght darts/half-edges.
    void generate_random_path(std::size_t lenght, CGAL::Random& random=CGAL::get_default_random());

    /// Creates a random closed path with at least lenght darts/half-edges.
    void generate_random_closed_path(std::size_t length, CGAL::Random& random=CGAL::get_default_random());
    
   };
}
