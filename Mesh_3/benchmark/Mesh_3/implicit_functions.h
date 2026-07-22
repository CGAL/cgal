struct Klein_function
{
  typedef ::FT           FT;
  typedef ::Point        Point;

  FT operator()(const Point& query) const
  {
    const FT x = query.x();
    const FT y = query.y();
    const FT z = query.z();

    return   (x*x+y*y+z*z+2*y-1)
           * ( (x*x+y*y+z*z-2*y-1) *(x*x+y*y+z*z-2*y-1)-8*z*z)
           + 16*x*z* (x*x+y*y+z*z-2*y-1);
  }
};

struct Tanglecube_function
{
  typedef ::FT           FT;
  typedef ::Point        Point;

  FT operator()(const Point& query) const
  {
    const FT x = query.x();
    const FT y = query.y();
    const FT z = query.z();

    double x2 = x*x, y2 = y*y, z2 = z*z;
    double x4 = x2*x2, y4 = y2*y2, z4 = z2*z2;
    return x4 - 5*x2 + y4 - 5*y2 + z4 - 5*z2 + 11.8;
  }
};

struct Sphere_function
{
  typedef ::FT           FT;
  typedef ::Point        Point;

  Sphere_function(double radius = 1.)
    : m_squared_radius(radius*radius)
  {}

  FT operator()(const Point& query) const
  {
    const FT x = query.x();
    const FT y = query.y();
    const FT z = query.z();

    return (x*x + y*y + z*z - m_squared_radius);
  }

protected:
  FT m_squared_radius;
};

struct Cylinder_function
{
  typedef ::FT           FT;
  typedef ::Point        Point;

  Cylinder_function(double radius = 0.5, double height = 2.)
    : m_radius(radius), m_height(height)
  {}

  FT operator()(const Point& query) const
  {
    const FT x = query.x();
    const FT y = query.y();
    const FT z = query.z();

    if(z > 0.5*m_height)
      return z - 0.5*m_height;
    else if(z < -0.5*m_height)
      return -z + 0.5*m_height;
    else
      return (x*x + y*y - m_radius*m_radius);
  }

protected:
  FT m_radius;
  FT m_height;
};