
#include <CGAL/Timer.h>

#include <nanoflann.hpp>

#include <boost/lexical_cast.hpp>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Memory_sizer.h>

using namespace std;
using namespace nanoflann;

template <typename T>
struct Point
{
  Point(T x, T y, T z)
    : x(x), y(y), z(z)
  {}

  T   x,y,z;
};

template <typename T>
ostream& operator <<(ostream& os, const Point<T>& p)
{
  os << p.x << " " << p.y << " " << p.z ;
  return os;
}


// This is an exampleof a custom data set class
template <typename T>
struct PointCloud
{
  std::vector<Point<T> >  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t size) const
	{
		const T d0=p1[0]-pts[idx_p2].x;
		const T d1=p1[1]-pts[idx_p2].y;
		const T d2=p1[2]-pts[idx_p2].z;
		return d0*d0+d1*d1+d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx].x;
		else if (dim==1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

template <typename T>
void generateRandomPointCloud(PointCloud<T> &point, istream& is, int n)
{
  T x, y, z;
  for(int i=0; i < n; i++){
    CGAL::read(is,x);
    CGAL::read(is,y);
    CGAL::read(is,z);
    
    point.pts.push_back(Point<T>(x,y,z));
  }
  std::cout << "Read "<< point.pts.size() << " points\n";
}

template <typename num_t>
void kdtree_demo(int argc, char** argv)
{
	PointCloud<num_t> cloud;
        int n;

	// Generate points:
        std::ifstream input(argv[1], std::ios::in | std::ios::binary);
        CGAL::set_binary_mode(input);
        //        input >> n >> n; // dimension and # of points
        CGAL::read(input,n);
        CGAL::read(input,n);
	generateRandomPointCloud(cloud, input, n);

        std::vector<Point<double> > queries;
        std::ifstream queries_stream(argv[2], std::ios::in | std::ios::binary);
        CGAL::set_binary_mode(queries_stream);
        CGAL::read(queries_stream,n);
        CGAL::read(queries_stream,n);
        // queries_stream >> n >> n;
        double x,y,z;
        for(int i=0; i < n; i++){
          CGAL::read(queries_stream,x);
          CGAL::read(queries_stream,y);
          CGAL::read(queries_stream,z);
          queries.push_back(Point<double>(x,y,z));
        }

        int runs = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 1;
        std::cerr << "runs = "  << runs <<std::endl;
 
        int bucketsize = (argc>4) ? boost::lexical_cast<int>(argv[4]) : 10;
        std::cerr << "bucketsize = "  << bucketsize <<std::endl;

	num_t query_pt[3] = { 0, 0, 0};

        CGAL::Timer timer;
        timer.start();
	// construct a kd-tree index:
	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<num_t, PointCloud<num_t> > ,
		PointCloud<num_t>,
		3 /* dim */
		> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(bucketsize /* max leaf */) );
	index.buildIndex();

        timer.stop();
        std::cout << "construction time: " << timer.time() << " sec" << std::endl;
        std::cerr << "Tree statistics:" << std::endl;
        std::cerr << "Number of items stored: "
          << index.items << std::endl;
        std::cerr << "Number of nodes: "
          << index.internals << std::endl;
        std::cerr << " Tree depth: " << index.depth() << std::endl;
	// do a knn search
	const size_t num_results = 10;
	size_t ret_index[num_results];
	num_t out_dist_sqr[num_results];
        int size = queries.size();
        
        std::cout << "start search" << std::endl;
        bool dump = true;
        double sum = 0;

        for(int i=0;i<runs;++i){
          
          nanoflann::KNNResultSet<num_t> resultSet(num_results);
         
          for(int i = 0 ; i < size; i++){
            query_pt[0] = queries[i].x;
            query_pt[1] = queries[i].y;
            query_pt[2] = queries[i].z;
	    resultSet.init(ret_index, out_dist_sqr );
            timer.reset();
          timer.start();
            index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10,0));
            timer.stop();

            for (int k=0; k<num_results; ++k){
	      if(dump)
	        std::cerr <<cloud.pts[ret_index[k]] << std::endl;
	    }
		    dump=false;
                    sum += timer.time();
          }
          
          
        }
        std::cerr << index.count_items <<" items\n";
        std::cerr << index.count_leafs <<" leaf\n";
        std::cerr << index.count_internals <<" internals visited\n";

        std::cerr<<std::endl << "total: " << sum << " sec\n";
        if(runs>1){
          std::cerr << "average: " << sum/runs << " sec\n";
        }
          std::cerr << "done\n";
}

int main(int argc, char** argv)
{
  kdtree_demo<double>(argc, argv);
  
  return 0;
}

