struct Visitor
{
  Visitor ( string audit_name ) 
    :
    mOut( new ofstream(audit_name.c_str(), ios::trunc ) ) 
  {
    CHECK_MSG(out(), "Unable to open audit file: " << audit_name);
  }
  
  void OnStarted( Surface& ) {} 
  
  void OnFinished ( Surface& ) {} 
  
  void OnStopConditionReached( Surface& ) {} 
  
  void OnCollected( Halfedge_handle const& aEdge, bool aIsFixed, Surface& )
  {
    out() << "I " << aEdge->id() << " " << aIsFixed << endl ;
  }                
  
  void OnSelected( Halfedge_handle const& aEdge, Surface&, optional<double> const& aCost, size_t, size_t )
  {
    out() << "S " << aEdge->id() << " " ;
    if ( aCost )
      out() << *aCost ;
    out() << endl ;  
  }                
  
  void OnCollapsing(Halfedge_handle const& aEdge, Surface&, optional<Point> const& aPlacement ) 
  {
    out() << "C " << aEdge->id() << " " ;
    if ( aPlacement )
      out() << aPlacement->x() << " " << aPlacement->y() << " " << aPlacement->z() ;
    out() << endl ;  
  }                
  
  void OnNonCollapsable(Halfedge_handle const& aEdge, Surface& ) 
  {
    out() << "N " << aEdge->id() << endl ;
  }                
  
  ostream& out() { return *mOut; }
  
  shared_ptr<ostream> mOut ;
} ;
