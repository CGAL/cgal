class Surface_mesh_geodesic_distance
{
   public:
      // METHODS FOR QUERYING DISTANCE =======================
      
      // Distance to a single source
      double getDistance( Vertex y );                         // get distance from the current source set to a vertex y
      void getDistance( VertexSet T, VertexValuePairs dist ); // get distance from the current source set to each vertex in T
      void getDistance( VertexFunction& dist );               // get distance from the current source set to all vertices

      // All-pairs distance query
      void allPairsDistance( VertexSet S, Matrix D ); // compute all distances between pairs of points in S

      // Furthest point queries
      void furthestPoint( Vertex& p, double& furthestDistance ); // find the point furthest from the current source set, and its distance; if there is not a unique furthest point, pick an arbitrary point

      // Distance gradients
      double getDistanceGradient( Vertex y );                         // get distance gradient with respect to the current source set to a vertex y
      void getDistanceGradient( VertexSet T, VertexValuePairs dist ); // get distance gradient with respect to the current source set to each vertex in T
      void getDistanceGradient( VertexFunction& dist );               // get distance gradient with respect to the current source set to all vertices

      // METHODS FOR SPECIFYING SOURCE SET ===================
      bool addSource( Vertex s ); // add s to the source set, returning false if s is already in the set
      bool removeSource( Vertex s ); // remove s to the source set, returning false if s is not in the set
      void getSources( VertexSet& S ); // get the list of current sources
      void clearSources(); // empty the source set
      SourcePointIterator sourcesBegin(); // iterator to beginning of source list
      SourcePointIterator sourcesEnd(); // iterator end end of source list

   protected:
      // precomputed data
      CholeskyFactorization Lheat; // factorization for Step I of heat method
      CholeskyFactorization Lpoisson; // factorization for Step II of heat method
      SurfaceMesh remesh; // intrinsic Delaunay mesh

      VertexSet sources;
};

