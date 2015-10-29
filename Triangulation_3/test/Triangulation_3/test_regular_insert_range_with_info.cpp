#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>              Traits;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, Traits>  Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,CGAL::Regular_triangulation_cell_base_3<Traits> > Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds>                     Regular;
typedef Traits::Weighted_point                                         Weighted_point;
typedef Traits::Bare_point                                             Point;

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Spatial_lock_grid_3<CGAL::Tag_priority_blocking>         Lock_ds;
typedef CGAL::Triangulation_data_structure_3<Vb, 
  CGAL::Regular_triangulation_cell_base_3<Traits>, CGAL::Parallel_tag> Tds_parallel;
typedef CGAL::Regular_triangulation_3<Traits, Tds_parallel, Lock_ds>   RT_parallel;
#endif

template <bool is_const>
void test_iterator_on_pair(){
  typedef std::vector< std::pair<Weighted_point,unsigned> > Container;
  typedef typename boost::mpl::if_< boost::mpl::bool_<is_const>,boost::add_const<Container>::type,Container >::type Cast_type;
  Container points;
  

  points.push_back( std::make_pair( Weighted_point(Point(0,0,0),1),0  ) );
  points.push_back( std::make_pair( Weighted_point(Point(1,0,0),2),1  ) );
  points.push_back( std::make_pair( Weighted_point(Point(0,1,0),3),2  ) );
  points.push_back( std::make_pair( Weighted_point(Point(0,0,1),4),3  ) );
  points.push_back( std::make_pair( Weighted_point(Point(2,2,2),5),4  ) );
  points.push_back( std::make_pair( Weighted_point(Point(-1,0,1),6),5 ) );

  
  Regular R( static_cast<Cast_type&>(points).begin(),static_cast<Cast_type&>(points).end() );

  assert( R.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Regular::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert(points[vit->info()].first == vit->point());

#ifdef CGAL_LINKED_WITH_TBB
  {
  // Construct the locking data-structure, using the bounding-box of the points
  typename RT_parallel::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);
  // Contruct the triangulation in parallel
  RT_parallel R(
    static_cast<Cast_type&>(points).begin(), 
    static_cast<Cast_type&>(points).end(),
    &locking_ds);

  assert(R.number_of_vertices() == 6);

  // check that the info was correctly set.
  RT_parallel::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert(points[vit->info()].first == vit->point());
  }
#endif
}

void toto(int){}

template <bool is_const>
void test_zip_iterator(){
  typedef std::vector< Weighted_point > Container;
  Container points;
  typedef typename boost::mpl::if_< boost::mpl::bool_<is_const>,boost::add_const<Container>::type,Container >::type Cast_type;
  
  points.push_back( Weighted_point(Point(0,0,0),1) );
  points.push_back( Weighted_point(Point(1,0,0),2) );
  points.push_back( Weighted_point(Point(0,1,0),3) );
  points.push_back( Weighted_point(Point(0,0,1),4) );
  points.push_back( Weighted_point(Point(2,2,2),5) );
  points.push_back( Weighted_point(Point(-1,0,1),6) );
  
  std::vector<unsigned> indices;
  indices.push_back(0);
  indices.push_back(1);
  indices.push_back(2);
  indices.push_back(3);
  indices.push_back(4);
  indices.push_back(5); 

  Regular R( boost::make_zip_iterator(boost::make_tuple( static_cast<Cast_type&>(points).begin(),indices.begin() )),
             boost::make_zip_iterator(boost::make_tuple( static_cast<Cast_type&>(points).end(),indices.end() ) )  );
  assert( R.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Regular::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert( points[ vit->info() ] == vit->point() );
  
#ifdef CGAL_LINKED_WITH_TBB
  {
  // Construct the locking data-structure, using the bounding-box of the points
  typename RT_parallel::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);
  // Contruct the triangulation in parallel
  RT_parallel R(
    boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).begin(), indices.begin())),
    boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).end(), indices.end())),
    &locking_ds);

  assert(R.number_of_vertices() == 6);

  // check that the info was correctly set.
  RT_parallel::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert(points[vit->info()] == vit->point());
  }
#endif
}

struct Auto_count : public std::unary_function<const Weighted_point&,std::pair<Weighted_point,unsigned> >{
  mutable unsigned i;
  Auto_count() : i(0){}
  std::pair<Weighted_point,unsigned> operator()(const Weighted_point& p) const {
    return std::make_pair(p,i++);
  }
};

template <bool is_const>
void test_transform_iterator(){
  typedef std::vector< Weighted_point > Container;
  Container points;
  typedef typename boost::mpl::if_< boost::mpl::bool_<is_const>,boost::add_const<Container>::type,Container >::type Cast_type;
  
  points.push_back( Weighted_point(Point(0,0,0),1) );
  points.push_back( Weighted_point(Point(1,0,0),2) );
  points.push_back( Weighted_point(Point(0,1,0),3) );
  points.push_back( Weighted_point(Point(0,0,1),4) );
  points.push_back( Weighted_point(Point(2,2,2),5) );
  points.push_back( Weighted_point(Point(-1,0,1),6) );

  Regular R( boost::make_transform_iterator(static_cast<Cast_type&>(points).begin(),Auto_count()),
             boost::make_transform_iterator(static_cast<Cast_type&>(points).end(),  Auto_count() )  );

  assert( R.number_of_vertices() == 6 );

  // check that the info was correctly set.
  Regular::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert( points[ vit->info() ] == vit->point() );  
  
#ifdef CGAL_LINKED_WITH_TBB
  {
  // Construct the locking data-structure, using the bounding-box of the points
  typename RT_parallel::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);
  // Contruct the triangulation in parallel
  RT_parallel R(
    boost::make_transform_iterator(static_cast<Cast_type&>(points).begin(), Auto_count()),
    boost::make_transform_iterator(static_cast<Cast_type&>(points).end(), Auto_count()),
    &locking_ds);

  assert(R.number_of_vertices() == 6);

  // check that the info was correctly set.
  RT_parallel::Finite_vertices_iterator vit;
  for (vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
    assert(points[vit->info()] == vit->point());
  }
#endif
}

int main()
{
  test_iterator_on_pair<false>();
  test_iterator_on_pair<true>();
  test_zip_iterator<false>();
  test_zip_iterator<true>();
  test_transform_iterator<false>();
  test_transform_iterator<true>();
  return 0;
}

