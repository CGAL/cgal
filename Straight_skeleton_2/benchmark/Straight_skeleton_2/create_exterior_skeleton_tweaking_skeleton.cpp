#include<vector>

#include<boost/shared_ptr.hpp>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_offset_polygons_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;
typedef Ss::Halfedge_const_iterator Halfedge_const_iterator ;

typedef boost::shared_ptr<Polygon_2> PolygonPtr ;
typedef boost::shared_ptr<Ss> SsPtr ;

typedef std::vector<PolygonPtr> PolygonPtrVector ;

int main()
{
  Polygon_2 poly ;
  poly.push_back( Point(192, 768) ) ;
  poly.push_back( Point(192, 704 ) ) ;
  poly.push_back( Point(224, 688 ) ) ;
  poly.push_back( Point(192, 672 ) ) ;
  poly.push_back( Point(192, 624 ) ) ;
  poly.push_back( Point(224, 624 ) ) ;
  poly.push_back( Point(240, 656 ) ) ;
  poly.push_back( Point(256, 624 ) ) ;
  poly.push_back( Point(288, 624 ) ) ;
  poly.push_back( Point(288, 768 ) ) ;
  poly.push_back( Point(256, 768 ) ) ;
  poly.push_back( Point(240, 736 ) ) ;
  poly.push_back( Point(224, 768 ) ) ;

  double offset = 144 ;
  SsPtr ss =
    CGAL::CGAL_SS_i::create_partial_exterior_straight_skeleton_2(
      offset,
      poly.vertices_begin(),
      poly.vertices_end(),
      K() );

  typedef CGAL::Polygon_offset_builder_traits_2<K>   OffsetBuilderTraits;
  typedef CGAL::Polygon_offset_builder_2<Ss,OffsetBuilderTraits,Polygon_2> OffsetBuilder;

  OffsetBuilder builder(*ss);

  for ( Halfedge_const_iterator i =  ss->halfedges_begin();
                                i != ss->halfedges_end();
                                ++i )
    if (  (i->is_bisector() && (i->id()%2)==0) ){
      Point p,q;

      if ( i->opposite()->vertex()->has_infinite_time() )
      {
        boost::optional<Point> op=
          builder.Construct_offset_point(offset , i->opposite());
        if(!op) continue;
        p=*op;
      }
      else
        p=i->opposite()->vertex()->point();

      if( i->vertex()->has_infinite_time() )
      {
        boost::optional<Point> op=
          builder.Construct_offset_point(offset , i);
        if(!op) continue;
        q=*op;
      }
      else
        q=i->vertex()->point();

      std::cout << "2 " << p << " 0 " <<  q << " 0" << std::endl;
    }
  return 0;
}


