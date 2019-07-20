///////////////// Definitions of several famous surfaces /////////////////
double sphere_function (double, double, double);  // (c=(0,0,0), r=1)
double ellipsoid_function (double, double, double);  // (c=(0,0,0), r=1)
double torus_function (double, double, double);  // (c=(0,0,0), r=2)
double chair_function (double, double, double);  // (c=(0,0,0), r=6)
double tanglecube_function (double, double, double);  // (c=(0,0,0), r=4)
double octic_function (double, double, double);  // (c=(0,0,0), r=2)
double heart_function (double, double, double);  // (c=(0,0,0), r=2)
double klein_function (double, double, double);  // (c=(0,0,0), r=4)
double ring_function (double, double, double);  // (c=(0,0,0), r=?)
double false_knot_function (double, double, double);  // (c=(0,0,0), r=1)
double knot1_function (double, double, double);  // (c=(0,0,0), r=4)
double knot2_function (double, double, double);  // (c=(0,0,0), r=4)
double knot3_function (double, double, double);  // (c=(0,0,0), r=4)
double cube_function (double, double, double);  // (c=(0,0,0), r=2)


template <int Sq_radius>
double sphere_function (double x, double y, double z) // (c=(0,0,0), r=Sq_radius)
{
  double x2=x*x, y2=y*y, z2=z*z;
  return (x2+y2+z2)/Sq_radius - 1;
}



template <typename FT, typename P>
class FT_to_point_function_wrapper : public CGAL::cpp98::unary_function<P, FT>
{
  typedef FT (*Implicit_function)(FT, FT, FT);
  Implicit_function function;
public:
  typedef P Point;
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
};
