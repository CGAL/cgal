namespace CGAL {
#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
  // Specialize CGAL::Kernel_traits<> for my point type(s).
  template <int dim>
  struct Kernel_traits < Point_float_d<dim> > {
    public:
    typedef My_kernel<dim> Kernel;
  };
#else
  //for compilers that doesn't support partial specialization
  //of class templates we provide the full specialization for
  //the specific types used
  struct Kernel_traits < Point_float_d<4> > {
    public:
    typedef My_kernel<4> Kernel;
  };
#endif
}


