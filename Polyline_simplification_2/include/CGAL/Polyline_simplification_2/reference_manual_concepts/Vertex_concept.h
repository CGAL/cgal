/// An implementation-defined model of this concept is the vertex of the 2D constrained Delaunay triangulation that is used internally
/// by the polyline simplification data structure.
struct PolylineSimplificationVertex : CGAL::Triangulation_vertex_base_with_id_2
{
    /// @return true if this vertex correspond to an endpoint of an open polyline
    bool is_terminal() const ;

    /// @return Sets whether the vertex is terminal or not.
    void set_is_terminal( bool is ) ;

    /// @return true if this vertex is shared by more than one polyline
    bool is_shared() const ;

    /// @return Sets whether the vertex is shared or not.
    void set_is_shared( bool is ) ;

    /// @return true if this vertex has been manually marked as fixed by the user
    bool is_fixed() const ;

    /// @return Sets whether the vertex is fixed or not.
    void set_is_fixed( bool is );
    
    /// @return The simplification cost for the vertex as computed by the 'PolylineSimplificationCostFunction'
    boost::optional<double> cost() const ;
  
};

