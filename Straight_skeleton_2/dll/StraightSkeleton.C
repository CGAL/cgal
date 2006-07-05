
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Straight_skeleton_2.h>
#include<CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::Point_2 Point_2;
typedef CGAL::Polygon_2<Kernel>    Contour;
typedef boost::shared_ptr<Contour> ContourPtr;
typedef std::vector<ContourPtr>    ContourSequence ;

typedef CGAL::Straight_skeleton_2<Kernel> Ss;

typedef Ss::Face_iterator Face_iterator;
typedef Ss::Halfedge_iterator Halfedge_iterator;
typedef Ss::Halfedge_handle   Halfedge_handle;
typedef Ss::Vertex_handle     Vertex_handle;

typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>      SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

typedef CGAL::Bbox_2 Bbox_2;

extern "C"
__declspec (dllexport)
bool 
StraightSkeleton(int np, int* np_i, double* xp, double* yp,
		 int& numFaces, int& numVertices, int*& numFace_i,
                double*& xf, double*& yf, bool dumpEPS)
{
  double scale = 1.0;
  SsBuilder ssb ;  

  Bbox_2 bbox;

  int currentPoint = 0;
  for(int i = 0; i < np; i++){
    std::vector<Point_2> points(np_i[i]);
    for(int j=0; j < np_i[i]; j++){
      Point_2 p(xp[currentPoint], yp[currentPoint]);
      if(currentPoint == 0){
	bbox = p.bbox() ;
      } else {
	bbox = bbox + p.bbox();
      }
      points[j] = p;
      ++currentPoint;
    }

    scale = 1000 / (bbox.xmax() - bbox.xmin()) ;
    if(! CGAL::is_simple_2(points.begin(),points.end())){
      std::cerr << "Polygon " << i << "  is not simple" << std::endl;
      return false;
    }
    if( ( (i == 0) && (CGAL::orientation_2(points.begin(),points.end()) != CGAL::COUNTERCLOCKWISE) )
	|| ( (i != 0) && (CGAL::orientation_2(points.begin(),points.end()) != CGAL::CLOCKWISE)) ){
      std::cerr << "Warning: polygon " << i << " has wrong orientation" << std::endl;
    ssb.enter_contour(points.rbegin(),points.rend());
    } else {
      ssb.enter_contour(points.begin(),points.end());
    }
  }  

  // Construct the skeleton
  boost::shared_ptr<Ss> ss = ssb.construct_skeleton();
    
  // Proceed only if the skeleton was correctly constructed.
  if ( ss )
    {
      //std::cerr <<  "We first count the points" << std::endl;
      numFaces= ss->size_of_faces();
      //std::cerr << "numFaces = " << numFaces << std::endl;
      numFace_i = new int[numFaces];
      if(numFace_i == NULL){
	std::cerr << "Allocation failed" << std::endl;
	return false;
      }
      numVertices = 0;
      int currentFace = 0;
      for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit){
	int count = 0;
	Halfedge_handle h = fit->halfedge();
	Halfedge_handle done;
	done = h;
	do {
	  count++;
	  h = h->next();
	} while(h != done);
	numVertices += count;
	numFace_i[currentFace] = count;
	//std::cerr << "numFace_i[" << currentFace << "] = " << numFace_i[currentFace] << std::endl;
	++currentFace;
      }
      //std::cerr <<  "numVertices = " << numVertices << std::endl;
        
      //std::cerr <<  " Allocate the x and y array and traverse the faces again" << std::endl;
      xf = new double[numVertices];
      if(xf == NULL){
	std::cerr << "Allocation failed" << std::endl;
	delete [] numFace_i;
	return false;
      }
      yf = new double[numVertices];
      
      if(yf == NULL){
	std::cerr << "Allocation failed" << std::endl;
	delete [] numFace_i;
	delete [] xf;
	return false;
      }
      int currentVertex = 0;
      for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit){
	Halfedge_handle h = fit->halfedge();
	Halfedge_handle done;
	done = h;
	do {
	  //	  std::cerr << "xf[" << currentVertex << "] = " << h->vertex()->point().x()  << std::endl;
	  //std::cerr << "yf[" << currentVertex << "] = " << h->vertex()->point().y()  << std::endl;
	  xf[currentVertex] =  h->vertex()->point().x();
	  yf[currentVertex] =  h->vertex()->point().y();
	  ++currentVertex;
	  h = h->next();
	} while(h != done);
      }

      if(dumpEPS){
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

	for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit){
	  Halfedge_handle h = fit->halfedge();
	  Halfedge_handle done;
	  done = h;
	  do {
	    if(h->is_bisector()){
	      dump << "skel\n";
	    } else {
	      dump << "cont\n";
	    }
	    dump << scale* h->vertex()->point().x() << " " << scale*h->vertex()->point().y() << " "
		 << scale*h->opposite()->vertex()->point().x() << " " << scale*h->opposite()->vertex()->point().y() << " E\n";
	    h = h->next();
	  } while(h != done);
	}

	dump << "grestore\nshowpage" << std::endl;

	dump.close();
      }

      std::cerr << "leave dll" << std::endl;

      return true;
    } else{
      std::cerr << "something went wrong" << std::endl;
      return false;
    }
}



extern "C"
__declspec (dllexport)

void
StraightSkeletonFree(int*& numFace_i,
		     double*& xf, double*& yf)
{
  if(xf != NULL) delete [] xf;
  if(yf != NULL) delete [] yf;
  if(numFace_i != NULL) delete [] numFace_i;
}

