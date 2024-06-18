#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double> SCd;
typedef CGAL::Simple_cartesian<CGAL::Exact_rational> SCrat;

typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false>> SCinterval;

template <typename K>
struct Convert {

    static auto get_type() {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            return typename K::Exact_kernel{};
        }
        else {
            return CGAL::Simple_cartesian<CGAL::Exact_rational>{};
        }
    }

        using Has = decltype(get_type());


   using Exact = std::conditional_t < std::is_floating_point_v<typename K::FT>, Has, CGAL::Simple_cartesian<CGAL::Exact_rational>>;

   typename Exact::Point_2 exact(const typename K::Point_2& p)
   {
       if constexpr (typename K::Has_filtered_predicates_tag::value) {
           std::cout << "Has_filtered_predicates_tag" << std::endl;
           return typename K::C2E()(p);
       }

       if constexpr (std::is_floating_point<typename K::FT>::type::value) {
           std::cout << "is_floating_point but not filtered" << std::endl;
           return Exact::Point_2(p.x(), p.y());
       }
       // We could check if K::FT is an interval and then check that is tight.

       // Otherwise we assume that K::FT converts to Rational
       return typename Exact::Point_2();

   }

SCinterval::Point_2 approx(const typename K::Point_2& p)
{
    std::cout << "___________________\n" << typeid(K).name() << std::endl;

    if constexpr (typename K::Has_filtered_predicates_tag::value) {
        std::cout << "Has_filtered_predicates_tag" << std::endl;
        return typename K::C2F()(p);
    }

    if constexpr (std::is_floating_point<typename K::FT>::type::value) {
        std::cout << "is_floating_point but not filtered" << std::endl;
        return SCinterval::Point_2(p.x(), p.y());
    }
    // Neither floating point, nor filtered
    return SCinterval::Point_2(CGAL::to_interval(p.x()), CGAL::to_interval(p.y()));

/*
    if constexpr (typename K::Has_filtered_predicates_tag::value) {
        std::cout << "Has_filtered_predicates_tag" << std::endl;
        std::cout << typeid(typename K::Approximate_kernel).name() << std::endl;
        std::cout << typeid(typename K::Exact_kernel::FT).name() << std::endl;
        std::cout << typeid(typename K::C2F).name() << std::endl;
        std::cout << typeid(typename K::C2E).name() << std::endl;
    }
    std::cout << std::endl;
    */
}



};

int main()
{
  SCinterval::Point_2 p;
  p = Convert<SCd>().approx(SCd::Point_2(78,12));
  p = Convert<SCrat>().approx(SCrat::Point_2(78, 12));
  p = Convert<CGAL::Epick>().approx(CGAL::Epick::Point_2(78,12));
  p = Convert<CGAL::Epeck>().approx(CGAL::Epeck::Point_2(78,12));


  Convert<SCd>::Exact::Point_2 pd = Convert<SCd>().exact(SCd::Point_2(0.3, 12));
  Convert<CGAL::Epick>::Exact::Point_2 pepick = Convert<CGAL::Epick>().exact(CGAL::Epick::Point_2(78, 12));
  return 0;
}
