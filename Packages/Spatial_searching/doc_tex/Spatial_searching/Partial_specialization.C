// define a small kernel for your point type
template <int dim>
class My_kernel
{
public:
  typedef double               FT;
  typedef double               RT;
  typedef Point_float_d<dim> Point_d;
};


namespace CGAL {
  // Specialize CGAL::Kernel_traits<> for my point type(s).
  template <int dim>
  class Kernel_traits < Point_float_d<dim> > {
    public:
    typedef My_kernel<dim> Kernel;
  };
}
