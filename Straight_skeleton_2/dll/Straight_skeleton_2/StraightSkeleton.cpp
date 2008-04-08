
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include<CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "StraightSkeleton.h"

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel>    Contour;
typedef boost::shared_ptr<Contour> ContourPtr;
typedef std::vector<ContourPtr>    ContourSequence ;

typedef CGAL::Straight_skeleton_2<Kernel> Ss;

typedef Ss::Face_iterator         Face_iterator;
typedef Ss::Face_handle           Face_handle;
typedef Ss::Halfedge_iterator     Halfedge_iterator;
typedef Ss::Halfedge_handle       Halfedge_handle;
typedef Ss::Halfedge_const_handle Halfedge_const_handle;
typedef Ss::Vertex_handle         Vertex_handle;
typedef Ss::Vertex_const_handle   Vertex_const_handle;

struct Visitor
{
  Visitor ( ProgressCallback progress ) : Progress(progress), mCurr(0), mTotal(0) {}

  void on_contour_edge_entered ( Halfedge_const_handle const& he ) const {}
  
  void on_initialization_started( std::size_t size_of_vertices ) const {}
  
  void on_initial_events_collected( Vertex_const_handle const& v, bool is_reflex, bool is_degenerate )  const
  {
    ++ mCurr ;
    if ( Progress )
    Progress(mCurr,mTotal);
  }

  void on_edge_event_created( Vertex_const_handle const& lnode
  , Vertex_const_handle const& rnode
  )  const {}

  void on_split_event_created( Vertex_const_handle const& node )  const {}

  void on_pseudo_split_event_created( Vertex_const_handle const& lnode
  , Vertex_const_handle const& rnode
  )  const {}

  void on_initialization_finished() const {}
  
  void on_propagation_started() const {}

  void on_anihiliation_event_processed ( Vertex_const_handle const& node0
  , Vertex_const_handle const& node1
  )  const {}


  void on_edge_event_processed( Vertex_const_handle const& lseed
  , Vertex_const_handle const& rseed
  , Vertex_const_handle const&  node
  )  const {} 

  void on_split_event_processed( Vertex_const_handle const& seed
  , Vertex_const_handle const& node0
  , Vertex_const_handle const& node1
  )  const {}

  void on_pseudo_split_event_processed( Vertex_const_handle const& lseed
  , Vertex_const_handle const& rseed
  , Vertex_const_handle const& node0
  , Vertex_const_handle const& node1
  )  const {}

  void on_vertex_processed( Vertex_const_handle const& node ) const 
  {
    if ( node->is_contour() )
    {
      ++ mCurr ;
      if ( Progress )
      Progress(mCurr,mTotal);
    }
  }

  void on_propagation_finished() const {}
  
  void on_cleanup_started() const {}
  
  void on_cleanup_finished() const {}
  
  void on_algorithm_finished ( bool finished_ok ) const {}

  void on_error( char const* msg )  const
  {
    std::cerr << msg << std::endl ;
  }

  public :


  void set_total ( int aTotal ) { mTotal = aTotal ; }

  ProgressCallback Progress ;
  mutable int mCurr ;
  int mTotal ;

} ;


typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>              SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss,Visitor> SsBuilder;

typedef CGAL::Bbox_2 Bbox_2;

extern "C"
{

  void __declspec (dllexport) StraightSkeletonFree(int* numFace_i, double* xf, double* yf)
  {
    delete [] xf;
    delete [] yf;
    delete [] numFace_i;
  }

  int __declspec (dllexport) StraightSkeleton( int np
                                             , int* np_i
                                             , double* xp
                                             , double* yp
                                             , int& numFaces
                                             , int& numVertices
                                             , int*& numFace_i
                                             , double*& xf
                                             , double*& yf
                                             , int dumpEPS
                                             , ProgressCallback progress
                                             )
  {
    int result = 0 ;

    numFace_i = NULL ;
    xf = yf = NULL ;

    try
    {
      double scale = 1.0;

      SsBuilderTraits traits ;
      Visitor         visitor(progress) ; 
      SsBuilder ssb(traits,visitor) ;  

      Bbox_2 bbox;
      
      int currentPoint = 0;
      for(int i = 0; i < np; i++)
      {
        int s = np_i[i];
        std::vector<Point_2> points(s);
        for(int j=0; j < s; j++)
        {
          Point_2 p(xp[currentPoint], yp[currentPoint]);
          if(currentPoint == 0)
          bbox = p.bbox() ;
          else bbox = bbox + p.bbox();

          points[j] = p;
          ++currentPoint;
        }

        if( ! CGAL::is_simple_2(points.begin(),points.end()))
        {
          std::cerr << "Polygon " << i << "  is not simple" << std::endl;
          return 0;
        }

        if(  CGAL::orientation_2(points.begin(),points.end()) != ( i == 0 ? CGAL::COUNTERCLOCKWISE 
                                                                          : CGAL::CLOCKWISE 
                                                                 ) 
            )
             ssb.enter_contour(points.rbegin(),points.rend());
        else ssb.enter_contour(points.begin(),points.end());
      }  
      

      visitor.set_total(currentPoint*2);

      // Construct the skeleton
      boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
      
      // Proceed only if the skeleton was correctly constructed.
      if ( ss )
      {
        // We first count the points
        numFaces= (int)ss->size_of_faces();
        numFace_i = new int[numFaces];

        numVertices = 0;
        int currentFace = 0;
        for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit)
        {
          int count = 0;
          Halfedge_handle h = fit->halfedge();
          Halfedge_handle done;
          done = h;
          do
          {
            count++;
            h = h->next();
          } while(h != done);
          numVertices += count;

          numFace_i[currentFace] = count;
          ++currentFace;
        }
        
        // Allocate the x and y array and traverse the faces again
        xf = new double[numVertices];
        yf = new double[numVertices];

        int currentVertex = 0;

        for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit)
        {
          Halfedge_handle h = fit->halfedge();
          Halfedge_handle done;
          done = h;
          do
          {
            xf[currentVertex] =  h->vertex()->point().x();
            yf[currentVertex] =  h->vertex()->point().y();
            ++currentVertex;
            h = h->next();
          }
          while(h != done);
        }

        if(dumpEPS)
        {
          scale = 1000 / (bbox.xmax() - bbox.xmin()) ;

          std::ofstream dump("dump.eps");
          dump << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" << scale* bbox.xmin()-1 << " " << scale* bbox.ymin()-1 << " " << scale*bbox.xmax()+1 << " "  << scale*bbox.ymax()+1 << std::endl;
          dump << "%%EndComments\n"
                  "gsave\n"
                  "1.0 setlinewidth\n"
                  "/cont { 0 0 0 setrgbcolor } bind def\n"
                  "/cont_w { 0.1 setlinewidth } bind def\n"
                  "/skel { 1 0 0 setrgbcolor } bind def\n"
                  "/skel_w { 1.0 setlinewidth } bind def\n"
                  "% stroke - x1 y1 x2 y2 E\n"
                  "/E {newpath moveto lineto stroke} bind def\n" << std::endl;

          for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit)
          {
            Halfedge_handle h = fit->halfedge();
            Halfedge_handle done;
            done = h;
            do
            {
              if(h->is_bisector())
                   dump << "skel\n";
              else dump << "cont\n";
              
              dump << scale* h->vertex()->point().x() << " " << scale*h->vertex()->point().y() << " "
                   << scale*h->opposite()->vertex()->point().x() << " " 
                   << scale*h->opposite()->vertex()->point().y() << " E\n";
              h = h->next();
            } 
            while(h != done);
          }

          dump << "grestore\nshowpage" << std::endl;

          dump.close();
        }

        result = 1 ;
      }
    }
    catch ( std::exception const& e ) 
    {
      std::cerr << "Exception thrown: " << e.what() << std::endl ;
      StraightSkeletonFree(numFace_i,xf,yf);  
    }
    catch ( ... ) 
    {
      std::cerr << "Unknown exception thrown." << std::endl ;
      StraightSkeletonFree(numFace_i,xf,yf);  
    }

    return result ;
  }

} // extern "C"


