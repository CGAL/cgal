/// This classcept represent s wrap per are passed to the polyline simplification algorithm to calculate
/// the "cost" of removing a vertex. Such a cost represents some measure of the deviation error between the
/// polyline sets before and after removal. The smaller the error the lower the cost, and the algoritm processes
/// vertices in increasing cost order to preserve the overall polyline set shape as much as possible
/// concept
class PolylineSimplificationNode
{
  
public:

  /// Returns the cost
  template<class Tr, class VertexDescriptor, class PointIterator>  
  boost::optional<double> operator()( Tr               const& aTr      
                                    , VertexDescriptor const& aPV
                                    , VertexDescriptor const& aQV
                                    , VertexDescriptor const& aRV
                                    , PointIterator           aOriginalSubpolyline_VerticesBegin
                                    , PointIterator           aOriginalSubpolyline_VerticesEnd
                                    ) const ;
  
};

