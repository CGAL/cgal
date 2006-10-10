struct Visitor
{
  Visitor ( ostream& o ) : mOut(&o) {}
  
  void OnStarted( Surface& ) {} 
  
  void OnFinished ( Surface& ) {} 
  
  void OnStopConditionReached( Surface& ) {} 
  
  void OnCollected( Halfedge_handle const& aEdge, bool aIsFixed, Surface& )
  {
    out() << "P" << aEdge->id() << " " << aIsFixed << endl ;
  }                
  
  void OnSelected( Halfedge_handle const& aEdge, Surface&, optional<double> const& aCost, size_t, size_t )
  {
    out() << "S" << aEdge->id() << " " ;
    if ( aCost )
      out() << *aCost ;
    out() << endl ;  
  }                
  
  void OnCollapsing(Halfedge_handle const& aEdge, Surface&, optional<Point> const& aPlacement ) 
  {
    out() << "C" << aEdge->id() << " " ;
    if ( aPlacement )
      out() << aPlacement->x() << " " << aPlacement->y() ;
    out() << endl ;  
  }                
  
  void OnNonCollapsable(Halfedge_handle const& aEdge, Surface& ) 
  {
    out() << "N" << aEdge->id() << endl ;
  }                
  
  ostream& out() { return *mOut; }
  
  ostream* mOut ;
} ;
