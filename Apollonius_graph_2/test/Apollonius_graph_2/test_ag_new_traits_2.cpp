#include <CGAL/basic.h>

#include <cassert>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/new_traits/Apollonius_graph_new_traits_2.h>

typedef CGAL::Sign Sign;

template <class GT, class GT_test>
class Check_traits : public GT
{
    GT_test test;
public:

    typedef typename GT::Site_2 Site_2;

    Check_traits (const GT &gt = GT(), const GT_test &gt_test = GT_test())
        : GT(), test (gt_test)
    {}

    struct Vertex_conflict_2 : public GT::Vertex_conflict_2
    {
        typename GT_test::Vertex_conflict_2 other;

        Vertex_conflict_2 (const typename GT::Vertex_conflict_2 &base,
                           const typename GT_test::Vertex_conflict_2 &test)
            : GT::Vertex_conflict_2(base), other(test)
        {}

        Sign operator() (const Site_2 &s1, const Site_2 &s2, const Site_2 &s3,
                const Site_2 &q) const
        {
            Sign r1 = GT::Vertex_conflict_2::operator() (s1, s2, s3, q);
            Sign r2 = other (s1, s2, s3, q);
            if (r1 != r2) {
                std::cerr
                    << "Vertex_conflict_2 "
                    << "(" << s1 << ")"
                    << "(" << s2 << ")"
                    << "(" << s3 << ")"
                    << "(" << q << ")"
                    << "should be " << r1 << ", not " << r2
                    << std::endl;
                abort();
            }
            return r1;
        }
        Sign operator() (const Site_2 &s1, const Site_2 &s2, const Site_2 &q) const
        {
            Sign r1 = GT::Vertex_conflict_2::operator() (s1, s2, q);
            Sign r2 = other (s1, s2, q);
            if (r1 != r2) {
                std::cerr
                    << "Vertex_conflict_2 "
                    << "(" << s1 << ")"
                    << "(" << s2 << ")"
                    << "(" << q << ")"
                    << "should be " << r1 << ", not " << r2
                    << std::endl;
                abort();
            }
            return r1;
        }
    };

    Vertex_conflict_2
        vertex_conflict_2_object() const
        {
            return Vertex_conflict_2 (GT::vertex_conflict_2_object(),
                    test.vertex_conflict_2_object());
        }
};

typedef CGAL::MP_Float NT;
typedef CGAL::Simple_cartesian<NT> K;
typedef CGAL::Apollonius_graph_traits_2<K> GT_old;
typedef CGAL::Apollonius_graph_new_traits_2<K> GT_new;
typedef Check_traits<GT_old,GT_new> GT;
typedef CGAL::Apollonius_graph_2<GT> AG;

typedef AG::Site_2 Site;

int main (void)
{

    std::ifstream is ("data/algo.dat");

    assert (is);

    std::cout << "Apollonius diagram construction" << std::endl;
    AG ag;
    Site s; while (is >> s) ag.insert (s);

    std::cout << "Apollonius diagram validation" << std::endl;
    assert (ag.is_valid());

    std::cout << "OK" << std::endl;
    return 0;
}
