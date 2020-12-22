#ifndef POLYHEDRON_DEMO_CREATE_SPHERE_H
#define POLYHEDRON_DEMO_CREATE_SPHERE_H

#include <vector>
#include <cmath>
#include <CGAL/number_type_config.h>

template <class FLOAT>
void create_flat_sphere(FLOAT R,
                        std::vector<FLOAT>& positions_spheres,
                        std::vector<FLOAT>& normals_spheres,
                        int prec = 36)
{
  //The more small they are, the more precise the Sphere will be.
  // Must be divisors of 180 and 360.
  const float rings=float(prec)/2;
  const float sectors=float(prec);
  const float to_rad = static_cast<float>(CGAL_PI / 180.0);

  float T, P;
  float x[4],y[4],z[4];

  using std::cos;
  using std::sin;

  //Top of the sphere
  for(float t=0; t<360.f; t+=sectors)
  {

    positions_spheres.push_back(0);
    positions_spheres.push_back(0);
    positions_spheres.push_back(R);

    normals_spheres.push_back(0);
    normals_spheres.push_back(0);
    normals_spheres.push_back(1);

    P = rings*to_rad;
    T = t*to_rad;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);
    positions_spheres.push_back(R * x[1]);
    positions_spheres.push_back(R * y[1]);
    positions_spheres.push_back(R * z[1]);

    normals_spheres.push_back(x[1]);
    normals_spheres.push_back(y[1]);
    normals_spheres.push_back(z[1]);

    //
    P = rings*to_rad;
    T = (t+sectors)*to_rad;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_spheres.push_back(R * x[2]);
    positions_spheres.push_back(R * y[2]);
    positions_spheres.push_back(R * z[2]);

    normals_spheres.push_back(x[2]);
    normals_spheres.push_back(y[2]);
    normals_spheres.push_back(z[2]);
  }

  //Body of the sphere
  for (float p=rings; p<180.f-rings; p+=rings)
    for(float t=0; t<360.f; t+=sectors)
    {
      //A
      P = p*to_rad;
      T = t*to_rad;
      x[0] = sin(P) * cos(T) ;
      y[0] = sin(P) * sin(T) ;
      z[0] = cos(P);

      positions_spheres.push_back(R * x[0]);
      positions_spheres.push_back(R * y[0]);
      positions_spheres.push_back(R * z[0]);

      normals_spheres.push_back(x[0]);
      normals_spheres.push_back(y[0]);
      normals_spheres.push_back(z[0]);

      //B
      P = (p+rings)*to_rad;
      T = t*to_rad;
      x[1] = sin(P) * cos(T) ;
      y[1] = sin(P) * sin(T) ;
      z[1] = cos(P);
      positions_spheres.push_back(R * x[1]);
      positions_spheres.push_back(R * y[1]);
      positions_spheres.push_back(R * z[1]);

      normals_spheres.push_back(x[1]);
      normals_spheres.push_back(y[1]);
      normals_spheres.push_back(z[1]);

      //C
      P = p*to_rad;
      T = (t+sectors)*to_rad;
      x[2] = sin(P) * cos(T) ;
      y[2] = sin(P) * sin(T) ;
      z[2] = cos(P);
      positions_spheres.push_back(R * x[2]);
      positions_spheres.push_back(R * y[2]);
      positions_spheres.push_back(R * z[2]);


      normals_spheres.push_back(x[2]);
      normals_spheres.push_back(y[2]);
      normals_spheres.push_back(z[2]);
      //D
      P = (p+rings)*to_rad;
      T = (t+sectors)*to_rad;
      x[3] = sin(P) * cos(T) ;
      y[3] = sin(P) * sin(T) ;
      z[3] = cos(P);
      positions_spheres.push_back(R * x[3]);
      positions_spheres.push_back(R * y[3]);
      positions_spheres.push_back(R * z[3]);


      normals_spheres.push_back(x[3]);
      normals_spheres.push_back(y[3]);
      normals_spheres.push_back(z[3]);

      positions_spheres.push_back(R * x[1]);
      positions_spheres.push_back(R * y[1]);
      positions_spheres.push_back(R * z[1]);

      normals_spheres.push_back(x[1]);
      normals_spheres.push_back(y[1]);
      normals_spheres.push_back(z[1]);

      positions_spheres.push_back(R * x[2]);
      positions_spheres.push_back(R * y[2]);
      positions_spheres.push_back(R * z[2]);

      normals_spheres.push_back(x[2]);
      normals_spheres.push_back(y[2]);
      normals_spheres.push_back(z[2]);

    }
  //Bottom of the sphere
  for(float t=0; t<360.f; t+=sectors)
  {
    positions_spheres.push_back(0);
    positions_spheres.push_back(0);
    positions_spheres.push_back(-R);

    normals_spheres.push_back(0);
    normals_spheres.push_back(0);
    normals_spheres.push_back(-1);


    P = (180-rings)*to_rad;
    T = t*to_rad;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);
    positions_spheres.push_back(R * x[1]);
    positions_spheres.push_back(R * y[1]);
    positions_spheres.push_back(R * z[1]);

    normals_spheres.push_back(x[1]);
    normals_spheres.push_back(y[1]);
    normals_spheres.push_back(z[1]);


    P = (180-rings)*to_rad;
    T = (t+sectors)*to_rad;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_spheres.push_back(R * x[2]);
    positions_spheres.push_back(R * y[2]);
    positions_spheres.push_back(R * z[2]);

    normals_spheres.push_back(x[2]);
    normals_spheres.push_back(y[2]);
    normals_spheres.push_back(z[2]);

  }
}

template <class FLOAT>
void create_flat_and_wire_sphere(FLOAT R,
                   std::vector<FLOAT>& positions_spheres,
                   std::vector<FLOAT>& normals_spheres,
                   std::vector<FLOAT>& positions_wire_spheres,
                   int prec = 36)
{
  //The smaller they are, the more precise the Sphere will be.
  // Must be a divisor of 360 and 180.
  const float rings=float(prec)/2;
  const float sectors=float(prec);
  const float to_rad = static_cast<float>(CGAL_PI / 180.0);

  create_flat_sphere(R, positions_spheres, normals_spheres);

  float T, P;
  float x[4],y[4],z[4];

  using std::cos;
  using std::sin;

  //Top of the sphere
  for(float t=0; t<360; t+=sectors)
  {
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(R);

    P = rings*to_rad;
    T = t*to_rad;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);

    //
    P = rings*to_rad;
    T = (t+sectors)*to_rad;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);

    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(R);
  }

  //Body of the sphere
  for (float p=rings; p<180.f-rings; p+=rings)
    for(float t=0; t<360.f; t+=sectors)
    {
      //A
      P = p*to_rad;
      T = t*to_rad;
      x[0] = sin(P) * cos(T) ;
      y[0] = sin(P) * sin(T) ;
      z[0] = cos(P);

      //B
      P = (p+rings)*to_rad;
      T = t*to_rad;
      x[1] = sin(P) * cos(T) ;
      y[1] = sin(P) * sin(T) ;
      z[1] = cos(P);

      //C
      P = p*to_rad;
      T = (t+sectors)*to_rad;
      x[2] = sin(P) * cos(T) ;
      y[2] = sin(P) * sin(T) ;
      z[2] = cos(P);

      //D
      P = (p+rings)*to_rad;
      T = (t+sectors)*to_rad;
      x[3] = sin(P) * cos(T) ;
      y[3] = sin(P) * sin(T) ;
      z[3] = cos(P);

      positions_wire_spheres.push_back(R * x[0]);
      positions_wire_spheres.push_back(R * y[0]);
      positions_wire_spheres.push_back(R * z[0]);


      positions_wire_spheres.push_back(R * x[1]);
      positions_wire_spheres.push_back(R * y[1]);
      positions_wire_spheres.push_back(R * z[1]);


      positions_wire_spheres.push_back(R * x[1]);
      positions_wire_spheres.push_back(R * y[1]);
      positions_wire_spheres.push_back(R * z[1]);



      positions_wire_spheres.push_back(R * x[3]);
      positions_wire_spheres.push_back(R * y[3]);
      positions_wire_spheres.push_back(R * z[3]);



      positions_wire_spheres.push_back(R * x[3]);
      positions_wire_spheres.push_back(R * y[3]);
      positions_wire_spheres.push_back(R * z[3]);



      positions_wire_spheres.push_back(R * x[2]);
      positions_wire_spheres.push_back(R * y[2]);
      positions_wire_spheres.push_back(R * z[2]);


      positions_wire_spheres.push_back(R * x[2]);
      positions_wire_spheres.push_back(R * y[2]);
      positions_wire_spheres.push_back(R * z[2]);


      positions_wire_spheres.push_back(R * x[0]);
      positions_wire_spheres.push_back(R * y[0]);
      positions_wire_spheres.push_back(R * z[0]);

    }
  //Bottom of the sphere
  for(float t=0; t<360.f; t+=sectors)
  {
    P = (180-rings)*to_rad;
    T = t*to_rad;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);

    P = (180-rings)*to_rad;
    T = (t+sectors)*to_rad;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(-R);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);


    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);

    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);


    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(-R);
  }
}

#endif //POLYHEDRON_DEMO_DRAW_SPHERE_H
