#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point          Point;

int main()
{
  {
    std::vector< std::pair<Point,unsigned> > points;
    points.push_back( std::make_pair(Point(0,0),0)  );
    points.push_back( std::make_pair(Point(1,0),1)  );
    points.push_back( std::make_pair(Point(0,1),2)  );
    points.push_back( std::make_pair(Point(14,4),3) );
    points.push_back( std::make_pair(Point(2,2),4)  );
    points.push_back( std::make_pair(Point(-4,0),5) );


    CDT T;
    T.insert( points.begin(),points.end() );

    assert( T.number_of_vertices() == 6 );

    // check that the info was correctly set.
    CDT::Finite_vertices_iterator vit;
    for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
      if( points[ vit->info() ].first != vit->point() ){
        std::cerr << "Error different info" << std::endl;
        exit(EXIT_FAILURE);
      }
    std::cout << "OK" << std::endl;
  }

  //test zip iterator
  {
    std::vector<unsigned> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    indices.push_back(3);
    indices.push_back(4);
    indices.push_back(5);

    std::vector<Point> points;
    points.push_back(Point(0,0));
    points.push_back(Point(1,0));
    points.push_back(Point(0,1));
    points.push_back(Point(1,47));
    points.push_back(Point(2,2));
    points.push_back(Point(-1,0));



    CDT T;
    T.insert( boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )  );

    assert( T.number_of_vertices() == 6 );


    // check that the info was correctly set.
    CDT::Finite_vertices_iterator vit;
    for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
      if( points[ vit->info() ] != vit->point() ){
        std::cerr << "Error different info" << std::endl;
        exit(EXIT_FAILURE);
      }
  }

  return 0;
}
