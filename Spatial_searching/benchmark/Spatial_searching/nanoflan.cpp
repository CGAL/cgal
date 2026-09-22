
#include <CGAL/Timer.h>

#include <nanoflann.hpp>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

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
        inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t /* size */) const
        {
                const T d0=p1[0]-pts[idx_p2].x;
                const T d1=p1[1]-pts[idx_p2].y;
                const T d2=p1[2]-pts[idx_p2].z;
                return d0*d0+d1*d1+d2*d2;
        }

        // Returns the dim-th component of the idx-th point in the class:
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
        bool kdtree_get_bbox(BBOX & /* bb */) const { return false; }

};

template <typename T>
void generateRandomPointCloud(PointCloud<T> &point, istream& is)
{
  T x, y, z;
  while(is >> x >> y  >> z){
    point.pts.push_back(Point<T>(x,y,z));
  }
  std::cout << "Read "<< point.pts.size() << " points\n";
}

template <typename num_t>
void kdtree_demo(const size_t /* N */)
{
        PointCloud<num_t> cloud;


        // Generate points:
        std::ifstream input("points.xyz");
        generateRandomPointCloud(cloud, input);

        std::vector<Point<double> > queries;
        std::ifstream queries_stream("queries.xyz");

        double x,y,z;
        while( queries_stream >> x >> y >> z){
          queries.push_back(Point<double>(x,y,z));
        }

        num_t query_pt[3] = { 0, 0, 0};

        CGAL::Timer timer;
        timer.start();
        // construct a kd-tree index:
        typedef KDTreeSingleIndexAdaptor<
                L2_Simple_Adaptor<num_t, PointCloud<num_t> > ,
                PointCloud<num_t>,
                3 /* dim */
                > my_kd_tree_t;

        my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index.buildIndex();

        timer.stop();
        std::cout << "construction time: " << timer.time() << " sec" << std::endl;

        timer.reset();
        timer.start();
        // do a knn search
        const size_t num_results = 50;
        size_t ret_index[50];
        num_t out_dist_sqr[50];
        int size = cloud.pts.size();

        std::cout << "start search" << std::endl;
        double sum = 0;
        for(int i = 0 ; i < size; i++){
          query_pt[0] = queries[i].x;
          query_pt[1] = queries[i].y;
          query_pt[2] = queries[i].z;
          nanoflann::KNNResultSet<num_t> resultSet(num_results);
          resultSet.init(ret_index, out_dist_sqr );
          index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10,0));

          //std::cout << "knnSearch(nn="<<num_results<<"): \n";
          //std::cout << "ret_index=" << ret_index << " out_dist_sqr=" << out_dist_sqr << endl;
          //std::cout << cloud.pts[ret_index] << std::endl;
          for (size_t k=0; k<num_results; ++k)
            sum += cloud.pts[ret_index[k]].x;
        }
        timer.stop();
        std::cout << sum << " done in " << timer.time() << " sec."<< std::endl;
}

int main(int /* argc */, char** /* argv */)
{
  //kdtree_demo<float>(100000);
        kdtree_demo<double>(100000);

        return 0;
}

