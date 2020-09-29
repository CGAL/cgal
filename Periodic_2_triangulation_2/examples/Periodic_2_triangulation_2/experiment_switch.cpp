#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#define CGAL_GENERIC_P2T2 // @todo still needed but to remove eventually

#include <CGAL/internal/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>

#include <CGAL/Random.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel                              K;
typedef K::Vector_2                                                                    Vector;

typedef typename CGAL::Periodic_2_offset_2                                             Offset;
typedef CGAL::Lattice_2<K>                                                             Lattice;
typedef CGAL::Periodic_2_triangulations_2::internal::Lattice_construct_point_2<K>      CP2;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_base_2<K, Offset, Lattice, CP2> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2_generic<GT>                       Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2_generic<GT>                         Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                                   Tds;

typedef CGAL::Periodic_2_Delaunay_triangulation_2_generic<GT, Tds>                     PDT;

typedef PDT::Vertex_handle                                                             Vertex_handle;
typedef PDT::Point                                                                     Point;

std::pair<double,double> get_mean_and_std(const std::vector<int>& v)
{
  double mean = std::accumulate(v.begin(), v.end(), 0) / double(v.size());
  int sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0);
  double std = std::sqrt(sq_sum / double(v.size()) - mean * mean);
  return std::make_pair(mean, std);
}

double get_median(std::vector<int> v)
{
  std::sort(v.begin(), v.end());
  int n = int(v.size());
  if(n % 2) {
    return v[n/2];
  } else {
    return (v[n/2 - 1] + v[n/2]) / 2.0f;
  }
}

int main(int, char**)
{
  std::ofstream log = std::ofstream("log.txt");
  std::ofstream res = std::ofstream("results.txt");
  std::ofstream agg = std::ofstream("aggregate.txt");

  auto rng = CGAL::Random(2995176);
  // auto rng = CGAL::Random();

  std::vector<double> lengths = {1.0, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 10};
  std::vector<double> skews = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  int ntrials = 200;

  for(double bl : lengths)
  {
    for(double skew : skews)
    {
      double bx = (skew / 2) * bl;
      CGAL::cpp11::array<Vector, 2> basis;
      basis = CGAL::make_array(Vector(bl, 0), Vector(bx, 1));
      Lattice lat = Lattice(basis[0], basis[1]);
      std::cout << bl << " " << skew << std::endl;
      agg << bl << " " << skew << " " << lat.systole_sq_length() << std::endl;
      res << "Length: " << bl << ", Skew: " << skew << ", Systole: " << lat.systole_sq_length() << std::endl;
      log << "Length: " << bl << ", Skew: " << skew << ", Systole: " << lat.systole_sq_length() << std::endl;

      std::vector<int> v_until_criterion;
      std::vector<int> v_until_simplicial;
      std::vector<int> v_until_last_nonsimplicial;

      for(int trial=0; trial<ntrials; ++trial)
      {
        log << "Trial " << trial << std::endl;
        bool criterion_fulfilled = false;
        int n_until_criterion = 1;
        int n_until_simplicial = -1;
        int n_until_last_nonsimplicial = 1;

        std::vector<Point> pts = { Point(rng.get_double(0, bl), rng.get_double()) };
        PDT T(pts.begin(), pts.end(), basis);

        log << T.number_of_vertices() << " " << T.is_simplicial_complex() << " " << T.too_big_faces() << std::endl;

        if(T.too_big_faces() == 0)
          criterion_fulfilled = true;

        while(!criterion_fulfilled)
        {
          T.insert(Point(rng.get_double(0, bl), rng.get_double()));
          ++n_until_criterion;

          bool is_simplicial = T.is_simplicial_complex();
          if(!is_simplicial)
            n_until_last_nonsimplicial = n_until_criterion;
          else if(n_until_simplicial == -1)
            n_until_simplicial = n_until_criterion;

          log << T.number_of_vertices() << " " << T.is_simplicial_complex() << " " << T.too_big_faces() << "\n";
          if(T.too_big_faces() == 0)
            criterion_fulfilled = true;
        }

        v_until_criterion.push_back(n_until_criterion);
        v_until_simplicial.push_back(n_until_simplicial);
        v_until_last_nonsimplicial.push_back(n_until_last_nonsimplicial);

        // Write results: until simplicial, until last non-simplicial, until criterion fulfilled.
        res << n_until_simplicial << " " << n_until_last_nonsimplicial << " " << n_until_criterion << std::endl;
      }

      // Write aggregate data: medians, means, stds
      agg << get_median(v_until_simplicial) << " "
          << get_median(v_until_last_nonsimplicial) << " "
          << get_median(v_until_criterion) << " ";

      auto ms1 = get_mean_and_std(v_until_simplicial);
      auto ms2 = get_mean_and_std(v_until_last_nonsimplicial);
      auto ms3 = get_mean_and_std(v_until_criterion);
      agg << ms1.first << " " << ms2.first << " " << ms3.first << " ";
      agg << ms1.second << " " << ms2.second << " " << ms3.second << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
