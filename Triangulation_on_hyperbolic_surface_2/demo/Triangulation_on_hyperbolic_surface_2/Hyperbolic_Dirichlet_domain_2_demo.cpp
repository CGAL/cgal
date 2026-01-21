#include "window.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>

#include <CGAL/Hyperbolic_Dirichlet_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

typedef CGAL::Gmpq                                                                                              NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                         Traits;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;

int main(int argc, char **argv)
{
    int seed = time(NULL);
    std::cout << "Using seed " << seed << std::endl;
    Factory factory;
    Domain domain = factory.make_hyperbolic_fundamental_domain_g2(seed);

    QApplication app(argc, argv);
    app.setApplicationName("Hyperbolic Dirichlet domain 2 Demo");
    DemoWindow window;
    window.item().draw_Dirichlet(domain);
    window.show();
    QStringList args = app.arguments();
    args.removeAt(0);
    return app.exec();
}