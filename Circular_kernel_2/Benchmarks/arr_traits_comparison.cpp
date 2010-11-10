#include <CGAL/Cartesian.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/IO/Dxf_variant_reader.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Lazy_circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <fstream>
#include <CGAL/Timer.h>

//bbox filtering
#include <CGAL/Filtered_bbox_circular_kernel_2.h>

typedef CGAL::Gmpq                                              Number_type;
typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>                         Lazy_nt;
typedef CGAL::Cartesian<Number_type>                            Kernel;        
typedef CGAL::Cartesian<Lazy_nt>                                Kernel_with_LK;
typedef CGAL::Lazy_kernel< Kernel >                             Lazy_Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Lazy_Kernel>          Lazy_Traits;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>               Traits;
typedef CGAL::Arr_circle_segment_traits_2<Kernel_with_LK>       Traits_with_LK;

typedef CGAL::Algebraic_kernel_for_circles_2_2<Number_type>     AK_exact;
typedef CGAL::Circular_kernel_2<Kernel,AK_exact>                CK_exact;

typedef CGAL::Interval_nt_advanced                              ANT;
typedef CGAL::Cartesian<ANT>                                    ALK;
typedef CGAL::Algebraic_kernel_for_circles_2_2<ANT>             AAK;
typedef CGAL::Circular_kernel_2 <ALK,AAK>                       ACK;
typedef CGAL::Lazy_circular_kernel_2<CK_exact,ACK>              Lazy_CK;

typedef CGAL::Algebraic_kernel_for_circles_2_2<Lazy_nt>         AK_lazy;
typedef CGAL::Circular_kernel_2<Kernel_with_LK,AK_lazy>         CK_with_LK;


typedef CGAL::Arr_circular_line_arc_traits_2<Lazy_CK>           Lazy_Traits_with_CK;
typedef CGAL::Arr_circular_line_arc_traits_2<CK_exact>          Traits_with_CK;
typedef CGAL::Arr_circular_line_arc_traits_2<CK_with_LK>        Traits_with_CK_with_LK;

typedef CGAL::Filtered_bbox_circular_kernel_2<CK_exact>         Fbb_CK_exact;
typedef CGAL::Arr_circular_line_arc_traits_2<Fbb_CK_exact>      Traits_Fbb_CK_exact;

typedef CGAL::Filtered_bbox_circular_kernel_2<Lazy_CK>          Fbb_Lazy_CK;
typedef CGAL::Arr_circular_line_arc_traits_2<Fbb_Lazy_CK>       Traits_Fbb_Lazy_CK;

template <class CK>
struct Get_exact{
  typename CK::Root_of_2::Exact_type operator()(const typename CK::Root_of_2& t) const {
    return CGAL::exact(t);
  }
};

template <>
struct Get_exact<CK_exact>{
  const CK_exact::Root_of_2& operator()(const CK_exact::Root_of_2& t) const {
    return t;
  }
};



template <class Traits,class Traits_with_CK>
struct Benchmark
{
  
  typedef CGAL::Arrangement_2<Traits>                             Arrangement_2;
  typedef CGAL::Arrangement_2<Traits_with_CK>                     Arrangement_2_with_CK;
  typedef typename Traits_with_CK::Kernel::Circular_arc_2         C2;
  typedef typename Traits_with_CK::Kernel::Line_arc_2             L2;
  typedef typename Traits_with_CK::Curve_2                        Curve_2_with_CK;
  typedef std::vector<Curve_2_with_CK>                            Arccontainer;
  typedef typename Traits_with_CK::Circular_arc_point_2           Circular_arc_point_2;
  typedef typename Traits::CoordNT                                         ORN;
  
  ORN convert (const typename Traits_with_CK::Kernel::Root_of_2& n)
  {
    Get_exact<typename Traits_with_CK::Kernel> exact;
    return CGAL::sign(exact(n).gamma())==CGAL::ZERO?ORN(exact(n).alpha()):ORN(exact(n).alpha(),exact(n).beta(),exact(n).gamma());
  }
  
  typename Traits::Point_2 convert(const Circular_arc_point_2& p)
  {
    ORN x=convert(p.x());
    ORN y=convert(p.y());
    return typename Traits::Point_2(x,y);
  }


  template <class Output_iterator>
  void read_from_dxf(const char* fname,Arccontainer& ac, Output_iterator out)
  {
    std::ifstream fin(fname);
    CGAL::variant_load<typename Traits_with_CK::Kernel, C2, L2>(fin,std::back_inserter(ac));
    fin.close();
    
    for (typename Arccontainer::iterator it= ac.begin();it!=ac.end();++it)
    {
      L2* seg=boost::get<L2>(&(*it));
      if (seg!=NULL){
        typename Traits::Point_2 s=convert(seg->source());
        typename Traits::Point_2 t=convert(seg->target());
        typename Traits::Kernel::Line_2 line(seg->supporting_line().a(),seg->supporting_line().b(),seg->supporting_line().c());
        *out++ = typename Traits::Curve_2( line,s,t );
      }
      else{
        C2* arc=boost::get<C2>(&(*it));
        assert(arc!=NULL);
        Circular_arc_point_2 ck_s=arc->source();
        Circular_arc_point_2 ck_t=arc->target();
        typename Traits::Kernel::Point_2 center(typename Traits::Kernel::FT(arc->supporting_circle().center().x()),
                               typename Traits::Kernel::FT(arc->supporting_circle().center().y()));
        typename Traits::Kernel::FT sqrad(arc->supporting_circle().squared_radius());
        *out++ = typename Traits::Curve_2(typename Traits::Kernel::Circle_2(center,sqrad,CGAL::COUNTERCLOCKWISE),convert(ck_s),convert(ck_t));
      }
    }
  }
  
  void run(const char* fname)
  {
    std::list<typename Traits::Curve_2>  curves;
    Arccontainer ac;

    read_from_dxf(fname,ac,std::back_inserter(curves));

  //using  Arr_circle_segment_traits_2
    {
      CGAL::Timer time; time.start();
      // Construct the arrangement of the curves.
      Arrangement_2    arr;
      insert (arr, curves.begin(), curves.end());
      // Print the size of the arrangement.
      std::cout << "(" << arr.number_of_vertices()
                << "," << arr.number_of_edges() 
                << "," << arr.number_of_faces() << ")" << std::endl;  
      time.stop();
      std::cout << "Time for Arr_circle_segment_traits_2 " <<  time.time() << "\n";
    }
  //using Arr_circular_line_arc_traits_2
    {
      CGAL::Timer time; time.start();
      Arrangement_2_with_CK      ck_arr;
      insert (ck_arr, ac.begin(), ac.end(),boost::false_type());
      std::cout << "(" << ck_arr.number_of_vertices()
                << "," << ck_arr.number_of_edges() 
                << "," << ck_arr.number_of_faces() << ")" << std::endl;  
      time.stop();
      std::cout << "Time for Arr_circular_line_arc_traits_2 " <<  time.time() << "\n";
    }
  }
  
};

int main( int argc, char* argv[] )
{
  std::cout << "Running with Lazy_kernel\n";
  Benchmark<Lazy_Traits,Lazy_Traits_with_CK> bench_lazy;
  bench_lazy.run(argv[1]);
  std::cout << "Running with Exact\n";
  Benchmark<Traits,Traits_with_CK> bench;
  bench.run(argv[1]);
  std::cout << "Running with Lazy_exact_nt\n";
  Benchmark<Traits_with_LK,Traits_with_CK_with_LK> bench_with_LK;
  bench_with_LK.run(argv[1]);
  return 0;
}

