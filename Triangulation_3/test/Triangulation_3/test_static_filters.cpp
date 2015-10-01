#define CGAL_PROFILE
#undef CGAL_NO_STATIC_FILTERS


#include <CGAL/FPU.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Random.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

typedef CGAL::Simple_cartesian<double>                                    Double_kernel;
typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> >           Exact_kernel;
typedef CGAL::Filtered_kernel<Double_kernel, true>                        FK_with_SF;
typedef CGAL::Filtered_kernel<Double_kernel, false>                       FK_without_SF;

typedef CGAL::Regular_triangulation_euclidean_traits_3<Exact_kernel>      Exact_traits;
typedef CGAL::Regular_triangulation_euclidean_traits_3<FK_without_SF>     FTr_without_SF;
typedef CGAL::Regular_triangulation_euclidean_traits_3<FK_with_SF>        FTr_with_SF;
typedef std::pair<double,double>                                          NT_pair;

template < class K1, class K2, class Cmp = CGAL::dont_check_equal >
class Regular_traits_checker : public CGAL::Kernel_checker<K1,K2,Cmp>
{
public:

    typedef CGAL::Comparison_result   Comparison_result;
    typedef CGAL::Oriented_side       Oriented_side;

    typedef K1                        Kernel1;
    typedef K2                        Kernel2;
    typedef Cmp                       Comparator;

    // Kernel objects are defined as pairs, with primitives run in parallel.
#define CGAL_kc_pair(X) typedef std::pair<typename K1::X, typename K2::X> X;

    CGAL_kc_pair(Weighted_point_3)

#undef CGAL_kc_pair

#define CGAL_Kernel_pred(X, Y) \
    typedef CGAL::Primitive_checker<typename K1::X, typename K2::X, Cmp> X; \
    X Y() const { return X(this->k1.Y(), this->k2.Y(), this->cmp); }

#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

public:
  
  CGAL_Kernel_pred(Compare_weighted_squared_radius_3,compare_weighted_squared_radius_3_object)
  CGAL_Kernel_pred(Power_test_3,power_test_3_object)
};


typedef Regular_traits_checker< FTr_with_SF,
                                FTr_without_SF>                             K3;

typedef K3::Weighted_point_3                                                Weighted_point_3;


CGAL::Random *r;

double rand_base()
{
  return r->get_double(0, 1);
}

// Random double almost in [0;1].
double my_rand()
{
  // Ensure 53 random bits, not 48.
  return CGAL_IA_FORCE_TO_DOUBLE(rand_base() + rand_base()/1024);
}

// Random point in unit cube.
Weighted_point_3 my_rand_wp3()
{
  double x = my_rand(), y = my_rand(), z = my_rand(), r=my_rand();
  return Weighted_point_3( FTr_with_SF::Weighted_point_3(FTr_with_SF::Bare_point(x, y, z),r) , FTr_without_SF::Weighted_point_3(FTr_without_SF::Bare_point(x, y, z),r) );
}

// Perturbation with given maximum relative epsilon.
void perturb(Weighted_point_3 &p, double rel_eps)
{
  double x = CGAL_IA_FORCE_TO_DOUBLE(p.first.x()*(1+rand_base()*rel_eps));
  double y = CGAL_IA_FORCE_TO_DOUBLE(p.first.y()*(1+rand_base()*rel_eps));
  double z = CGAL_IA_FORCE_TO_DOUBLE(p.first.z()*(1+rand_base()*rel_eps));
  double r = CGAL_IA_FORCE_TO_DOUBLE(p.first.weight()*(1+rand_base()*rel_eps));
  p=Weighted_point_3( FTr_with_SF::Weighted_point_3(FTr_with_SF::Bare_point(x, y, z),r) , FTr_without_SF::Weighted_point_3(FTr_without_SF::Bare_point(x, y, z),r) );
}

void toto (int){}

void test_compare_weighted_squared_radius_3(){
  Weighted_point_3 p=my_rand_wp3();
  Weighted_point_3 q=my_rand_wp3();
  Weighted_point_3 r=my_rand_wp3();
  Weighted_point_3 s=my_rand_wp3();
  double alpha=my_rand();
  
  
  //test with random points + random alpha
  K3().compare_weighted_squared_radius_3_object()(p,q,r,s,NT_pair(alpha,alpha));
  K3().compare_weighted_squared_radius_3_object()(p,q,r,  NT_pair(alpha,alpha));
  K3().compare_weighted_squared_radius_3_object()(p,q,    NT_pair(alpha,alpha));
  K3().compare_weighted_squared_radius_3_object()(p,      NT_pair(alpha,alpha));

  CGAL::Weighted_converter_3<CGAL::Cartesian_converter<FK_with_SF,Exact_kernel>,FTr_with_SF,Exact_traits >    convert_to_exact;
  Exact_traits::Weighted_point_3 p_e=convert_to_exact(p.first);
  Exact_traits::Weighted_point_3 q_e=convert_to_exact(q.first);
  Exact_traits::Weighted_point_3 r_e=convert_to_exact(r.first);
  Exact_traits::Weighted_point_3 s_e=convert_to_exact(s.first);  
  
  Exact_traits::Compute_squared_radius_smallest_orthogonal_sphere_3 radius=Exact_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object();  
  double alpha_pqrs=CGAL::to_double( radius(p_e,q_e,r_e,s_e) );
  double alpha_pqr =CGAL::to_double( radius(p_e,q_e,r_e) );
  double alpha_pq  =CGAL::to_double( radius(p_e,q_e) );
  double alpha_p   = - p.first.weight();
  
  //test with random points + alpha limit
  K3().compare_weighted_squared_radius_3_object()(p,q,r,s,NT_pair(alpha_pqrs,alpha_pqrs));
  K3().compare_weighted_squared_radius_3_object()(p,q,r,  NT_pair(alpha_pqr,alpha_pqr));
  K3().compare_weighted_squared_radius_3_object()(p,q,    NT_pair(alpha_pq,alpha_pq));
  K3().compare_weighted_squared_radius_3_object()(p,      NT_pair(alpha_p,alpha_p));
  //test correct result
  assert( K3().compare_weighted_squared_radius_3_object()(p,q,r,s,NT_pair(alpha_pqrs-0.1,alpha_pqrs-0.1)) ==CGAL::POSITIVE );
  assert( K3().compare_weighted_squared_radius_3_object()(p,q,r,  NT_pair(alpha_pqr-0.1,alpha_pqr-0.1)) ==CGAL::POSITIVE );
  assert( K3().compare_weighted_squared_radius_3_object()(p,q,    NT_pair(alpha_pq-0.1,alpha_pq-0.1)) ==CGAL::POSITIVE );
  assert( K3().compare_weighted_squared_radius_3_object()(p,      NT_pair(alpha_p-0.1,alpha_p-0.1)) ==CGAL::POSITIVE );

  // Then with some perturbation on coordinates and weight.
  perturb(p, 1.0/(1<<25)/(1<<20)); // 2^-45
  
  K3().compare_weighted_squared_radius_3_object()(p,q,r,s,NT_pair(alpha_pqrs,alpha_pqrs));
  K3().compare_weighted_squared_radius_3_object()(p,q,r,  NT_pair(alpha_pqr,alpha_pqr));
  K3().compare_weighted_squared_radius_3_object()(p,q,    NT_pair(alpha_pq,alpha_pq));
  K3().compare_weighted_squared_radius_3_object()(p,      NT_pair(alpha_p,alpha_p));
  
}

Weighted_point_3 convert_to_pair(const Exact_traits::Weighted_point_3& wp)
{
  CGAL::Weighted_converter_3<CGAL::Cartesian_converter<Exact_kernel,FK_with_SF>,Exact_traits,FTr_with_SF >    convert_to_double;
  CGAL::Weighted_converter_3<CGAL::Cartesian_converter<FK_with_SF,FK_without_SF>,FTr_with_SF,FTr_without_SF > convert_to_double_noSF;
  
  FTr_with_SF::Weighted_point_3 wp_with_sf = convert_to_double(wp);
  FTr_without_SF::Weighted_point_3 wp_without_sf = convert_to_double_noSF(wp_with_sf);
  return Weighted_point_3(wp_with_sf, wp_without_sf);
}

void test_power_test_3(){
  CGAL::Weighted_converter_3<CGAL::Cartesian_converter<FK_without_SF,Exact_kernel>,FTr_without_SF,Exact_traits > convert_to_exact;
  
  Weighted_point_3 p=my_rand_wp3();
  Weighted_point_3 q=my_rand_wp3();
  Weighted_point_3 r=my_rand_wp3();
  Weighted_point_3 s=my_rand_wp3();  
  Weighted_point_3 query_pt=my_rand_wp3();  

  //test with random points
  K3().power_test_3_object()(p,q,r,s,query_pt);  
  K3().power_test_3_object()(p,q,r,query_pt);  
  K3().power_test_3_object()(p,q,query_pt);  
  K3().power_test_3_object()(p,query_pt);
  
  //test in degenerate case
  Exact_traits::Weighted_point_3::Point origin(0,0,0);
  Exact_traits::Weighted_point_3 p_e=convert_to_exact(p.second);
  Exact_traits::Weighted_point_3 q_e=convert_to_exact(q.second);
  Exact_traits::Weighted_point_3 r_e=convert_to_exact(r.second);
  Exact_traits::Weighted_point_3 s_e=convert_to_exact(s.second);
  
  Exact_traits::Weighted_point_3 tmp=Exact_traits::Weighted_point_3(
    origin,
      CGAL::squared_distance(origin,Exact_traits().construct_weighted_circumcenter_3_object()(p_e,q_e,r_e,s_e))-
      Exact_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(p_e,q_e,r_e,s_e)
  );
  assert(Exact_traits().power_test_3_object()(p_e,q_e,r_e,s_e,tmp)==CGAL::ON_ORIENTED_BOUNDARY);
  
  Weighted_point_3 ortho_pqrs = convert_to_pair(tmp);
  tmp=Exact_traits::Weighted_point_3(
      Exact_traits().construct_weighted_circumcenter_3_object()(p_e,q_e,r_e),
      -Exact_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(p_e,q_e,r_e)
    );
  assert(Exact_traits().power_test_3_object()(p_e,q_e,r_e,tmp)==CGAL::ON_ORIENTED_BOUNDARY);    
  Weighted_point_3 ortho_pqr = convert_to_pair(tmp);
  tmp=Exact_traits::Weighted_point_3(
      Exact_traits().construct_weighted_circumcenter_3_object()(p_e,q_e),
      -Exact_traits().compute_squared_radius_smallest_orthogonal_sphere_3_object()(p_e,q_e)
    );
  assert(Exact_traits().power_test_3_object()(p_e,q_e,tmp)==CGAL::ON_ORIENTED_BOUNDARY);        
  Weighted_point_3 ortho_pq = convert_to_pair(tmp);
  
  
  K3().power_test_3_object()(p,q,r,s,ortho_pqrs);
  K3().power_test_3_object()(p,q,r  ,ortho_pqr);
  K3().power_test_3_object()(p,q    ,ortho_pq);
  // Then with some perturbation on coordinates and weight.
  perturb(p, 1.0/(1<<25)/(1<<20)); // 2^-45
  
  K3().power_test_3_object()(p,q,r,s,ortho_pqrs);
  K3().power_test_3_object()(p,q,r  ,ortho_pqr);
  K3().power_test_3_object()(p,q    ,ortho_pq);
}


int main(int argc, char **argv)
{
  assert(!   Exact_traits::Has_filtered_predicates);
  assert(    FTr_with_SF::Has_filtered_predicates);
  assert(    FTr_without_SF::Has_filtered_predicates);
  

  assert(!  FTr_without_SF::Has_static_filters);
  assert(   FTr_with_SF::Has_static_filters);

  
  
  int loops = (argc < 2) ? 2000 : std::atoi(argv[1]);
  int seed  = (argc < 3) ? CGAL::get_default_random().get_int(0, 1<<30)
                         : std::atoi(argv[2]);

  std::cout << "Initializing random generator with seed = " << seed
            << std::endl
            << "#loops = " << loops << " (can be changed on the command line)"
            << std::endl;

  CGAL::Random rnd(seed);
  r = &rnd;

  std::cout.precision(20);
  std::cerr.precision(20);

  std::cout << "ulp(1) = " << CGAL::internal::Static_filter_error::ulp() << std::endl;

  std::cout << "Testing Compare_weighted_squared_radius_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_compare_weighted_squared_radius_3();  

  std::cout << "Testing Power_test_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_power_test_3();  

  return 0;
}












