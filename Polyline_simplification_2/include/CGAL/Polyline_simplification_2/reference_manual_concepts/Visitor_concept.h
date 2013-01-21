/// Modles of this concept represent a visitation function object whose methods are called at key steps during
/// the simplification process.
struct PolylineSimplificationVisitor
{
  /// Called at the very beggining, before vertices are being collected for processing
  void OnStarted() const {} 
  
  /// Called at the very end when the simplification finished
  void OnFinished() const {} 
  
  /// Called when the 'PolylineSimplificationStopPredicate' returned true
  void OnStopConditionReached() const {} 
  
  /// Called when a vertex that has been alreay classified as "removable" is put in the processing queue.
  /// @param vertex the collected vertex (as a 'PolylineSimplificationVertex')
  template<class VertexHandle>
  void OnCollected( VertexHandle const& vertex ) const {}                
  
  /// Called when a vertex has been popped off the processing queue and will be removed
  /// @param vertex the processed vertex
  /// @param cost   the cost of removing the current vertex as calculated by the 'PolylineSimplificationCostFunction'
  /// @param initial_count the initial total number of vertices in the polyline set
  /// @param current_count the current total number of vertices in the polyline set
  template<class VertexHandle>
  void OnSelected( VertexHandle const& vertex, boost::optional<double> const& cost, unsigned initial_count, unsigned current_count) const {}                
  
  /// Called just before a selected vertex is removed
  /// @param p the previous vertex along the polyline
  /// @param q the vetex about to be remove
  /// @param r the next vertex along the polyline
  template<class VertexHandle>
  void OnRemoving( VertexHandle const& p, VertexHandle const& q, VertexHandle const& r) const {}          
  
  /// Called right after a selected vertex has been removed
  /// @param p the remaning vertex along the polyline that was right before "q" (now removed)
  /// @param r the remaning vertex along the polyline that was right after  "q" (now removed)
  template<class VertexHandle>
  void OnRemoved( VertexHandle const& p, VertexHandle const& r) const {}        
  
  /// Called when a vertex that has been classified as "non-removable"
  /// @param vertex the non-removable vertex 
  template<class VertexHandle>
  void OnNonRemovable( VertexHandle const& vertex) const {}                

};    
