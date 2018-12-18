
namespace CGAL {
  /*!
    \ingroup PkgSurfaceMeshTopologyClasses
    
    The class `Path_on_surface` represents XXX
    
    \tparam Map XXX
  */
  template<typename Map_>
  class Path_on_surface
  {
  public:
    /*! XXX
     */
    Path_on_surface(const Map& amap);

    /// @return true iff the path is empty
    bool is_empty() const;

    /// @return true iff the path is closed.
    bool is_closed() const;

    /// @return true iff the path does not pass twice through a same edge
    ///              or a same vertex.
    bool is_simple() const;

  /// clear the path.
    void clear();
    
    /// @return true iff df can be added at the end of the path.
    bool can_be_pushed(Dart_const_handle dh) const;
  
    /// Add the given dart at the end of this path.
    /// @pre can_be_pushed(dh)
    void push_back(Dart_const_handle dh, bool update_isclosed=true);

    /*!
     */
    Self& operator+=(const Self& other);

    /*!
     */
    void generate_random_path(std::size_t length, CGAL::Random& random=CGAL::get_default_random());

     /*!
     */
    void generate_random_closed_path(std::size_t length, CGAL::Random& random=CGAL::get_default_random());
    
    /// Reverse the path (i.e. negate its orientation).
    void reverse();
   };
}
