#include "test_dependencies.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/functional.h>

#include <iostream>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel            K;

template <typename K>
struct Tester
{
  typedef K                                                                 Traits;
  typedef CGAL::Regular_triangulation_vertex_base_3<Traits>                 Vbb;
  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, Traits,Vbb> Vb;
  typedef CGAL::Regular_triangulation_cell_base_3<Traits>                   Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                      Tds;

  typedef CGAL::Regular_triangulation_3<Traits, Tds>                        RT;
  typedef typename RT::Bare_point                                           Bare_point;
  typedef typename RT::Weighted_point                                       Weighted_point;

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Spatial_lock_grid_3<CGAL::Tag_priority_blocking>           Lock_ds;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag> Tds_parallel;
  typedef CGAL::Regular_triangulation_3<Traits, Tds_parallel, Lock_ds>     RT_parallel;
#endif

  template <bool is_const>
  void test_iterator_on_pair() const
  {
    typedef std::vector<std::pair<Weighted_point, unsigned> >              Container;
    typedef std::conditional_t<is_const,
                               std::add_const_t<Container>,
                               Container>                                  Cast_type;

    Container points;
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.160385, 0.599679, 0.374932), -0.118572), 0));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.17093, 0.82228, 0.51697), -0.0226131), 1));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.15773, 0.66293, 0.458541), -0.0753074), 2));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.388417, 0.685989, 0.401349), -0.0616195), 3));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.380061, 0.852124, 0.538984), -0.0145638), 4));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.0402467, 0.519724, 0.417205), -0.0698374), 5));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.0270472, 0.360373, 0.358776), -0.0109982), 6));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.257734, 0.383432, 0.301584), -0.0601458), 7));
    points.push_back(std::make_pair(Weighted_point(Bare_point(0.142091, 0.643406, 0.61943), -0.251061), 8));

    RT R(static_cast<Cast_type&>(points).begin(), static_cast<Cast_type&>(points).end());
    assert(R.number_of_vertices() == 9);

    R.clear();
    R.insert(static_cast<Cast_type&>(points).begin(), static_cast<Cast_type&>(points).end());
    assert(R.number_of_vertices() == 9);

    // check that the info was correctly set.
    typename RT::Finite_vertices_iterator vit;
    for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
      assert(points[vit->info()].first == vit->point());

  #ifdef CGAL_LINKED_WITH_TBB
    {
      // Construct the locking data-structure, using the bounding-box of the points
      typename RT_parallel::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);

      // Construct the triangulation in parallel
      RT_parallel R(static_cast<Cast_type&>(points).begin(), static_cast<Cast_type&>(points).end(), &locking_ds);
      assert(R.number_of_vertices() == 9);

      R.clear();
      R.insert(static_cast<Cast_type&>(points).begin(), static_cast<Cast_type&>(points).end(), &locking_ds);
      assert(R.number_of_vertices() == 9);

      // check that the info was correctly set.
      typename RT_parallel::Finite_vertices_iterator vit;
      for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit) {
        assert(points[vit->info()].first == vit->point());
      }
    }
#endif
  }

  template <bool is_const>
  void test_zip_iterator() const
  {
    typedef std::vector<Weighted_point>                                   Container;
    typedef std::conditional_t<is_const,
                               std::add_const_t<Container>,
                               Container >                                Cast_type;

    Container points;
    points.push_back(Weighted_point(Bare_point(0,0,0),1));
    points.push_back(Weighted_point(Bare_point(1,0,0),2));
    points.push_back(Weighted_point(Bare_point(0,1,0),3));
    points.push_back(Weighted_point(Bare_point(0,0,1),4));
    points.push_back(Weighted_point(Bare_point(2,2,2),5));
    points.push_back(Weighted_point(Bare_point(-1,0,1),6));

    std::vector<unsigned> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    indices.push_back(3);
    indices.push_back(4);
    indices.push_back(5);

    RT R(boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).begin(), indices.begin())),
         boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).end(), indices.end())));
    assert(R.number_of_vertices() == 6);

    // check that the info was correctly set.
    typename RT::Finite_vertices_iterator vit;
    for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
      assert(points[ vit->info() ] == vit->point());

#ifdef CGAL_LINKED_WITH_TBB
    {
      // Construct the locking data-structure, using the bounding-box of the points
      typename RT_parallel::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);

      // Construct the triangulation in parallel
      RT_parallel R(boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).begin(), indices.begin())),
                    boost::make_zip_iterator(boost::make_tuple(static_cast<Cast_type&>(points).end(), indices.end())),
                    &locking_ds);
      assert(R.number_of_vertices() == 6);

      // check that the info was correctly set.
      typename RT_parallel::Finite_vertices_iterator vit;
      for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit) {
        assert(points[vit->info()] == vit->point());
      }
    }
#endif
  }

  struct Auto_count
    : public CGAL::cpp98::unary_function<const Weighted_point&,
                                         std::pair<Weighted_point,
                                         unsigned> >
  {
    mutable unsigned i;
    Auto_count() : i(0){}
    std::pair<Weighted_point,unsigned> operator()(const Weighted_point& p) const
    {
      return std::make_pair(p,i++);
    }
  };

  template <bool is_const>
  void test_transform_iterator() const
  {
    typedef std::vector< Weighted_point >                                   Container;
    typedef std::conditional_t<is_const,
                               std::add_const_t<Container>,
                               Container >                                  Cast_type;

    Container points;
    points.push_back(Weighted_point(Bare_point(0,0,0),1));
    points.push_back(Weighted_point(Bare_point(1,0,0),2));
    points.push_back(Weighted_point(Bare_point(0,1,0),3));
    points.push_back(Weighted_point(Bare_point(0,0,1),4));
    points.push_back(Weighted_point(Bare_point(2,2,2),5));
    points.push_back(Weighted_point(Bare_point(-1,0,1),6));

    RT R(boost::make_transform_iterator(static_cast<Cast_type&>(points).begin(), Auto_count()),
         boost::make_transform_iterator(static_cast<Cast_type&>(points).end(), Auto_count()));

    assert(R.number_of_vertices() == 6);

    // check that the info was correctly set.
    typename RT::Finite_vertices_iterator vit;
    for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit)
      assert(points[ vit->info() ] == vit->point());

#ifdef CGAL_LINKED_WITH_TBB
    {
      // Construct the locking data-structure, using the bounding-box of the points
      typename RT_parallel::Lock_data_structure locking_ds(CGAL::Bbox_3(-1., 0., 0., 2, 2, 2), 50);

      // Construct the triangulation in parallel
      RT_parallel R(boost::make_transform_iterator(static_cast<Cast_type&>(points).begin(), Auto_count()),
                    boost::make_transform_iterator(static_cast<Cast_type&>(points).end(), Auto_count()),
                    &locking_ds);
      assert(R.number_of_vertices() == 6);

      // check that the info was correctly set.
      typename RT_parallel::Finite_vertices_iterator vit;
      for(vit = R.finite_vertices_begin(); vit != R.finite_vertices_end(); ++vit) {
        assert(points[vit->info()] == vit->point());
      }
    }
#endif
  }

  void operator()() const
  {
    test_iterator_on_pair<false>();
    test_iterator_on_pair<true>();
    test_zip_iterator<false>();
    test_zip_iterator<true>();
    test_transform_iterator<false>();
    test_transform_iterator<true>();
  }
};

int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
  typedef CGAL::Exact_predicates_exact_constructions_kernel   Epeck;

  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<Epick> test_epic;
  test_epic();

  std::cerr << "TESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  Tester<Epeck> test_epec;
  test_epec();

  return EXIT_SUCCESS;
}
