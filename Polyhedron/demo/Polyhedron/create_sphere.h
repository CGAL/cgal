#ifndef POLYHEDRON_DEMO_CREATE_SPHERE_H
#define POLYHEDRON_DEMO_CREATE_SPHERE_H

#include <vector>
#include <cmath>

template <class FLOAT>
void create_flat_sphere(double R,
                        std::vector<FLOAT>& positions_spheres,
                        std::vector<FLOAT>& normals_spheres,
                        bool four_pos=false)
{
  //The more small they are, the more precise the Sphere will be.
  // Must be a multiple of 360 and 180.
  const int rings=18;
  const int sectors=38;

  float T, P;
  float x[4],y[4],z[4];

  //Top of the sphere
  for(int t=0; t<360; t+=sectors)
  {

    positions_spheres.push_back(0);
    positions_spheres.push_back(0);
    positions_spheres.push_back(R);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(0);
    normals_spheres.push_back(0);
    normals_spheres.push_back(1);

    P = rings*M_PI/180.0;
    T = t*M_PI/180.0;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);
    positions_spheres.push_back(R * x[1]);
    positions_spheres.push_back(R * y[1]);
    positions_spheres.push_back(R * z[1]);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(x[1]);
    normals_spheres.push_back(y[1]);
    normals_spheres.push_back(z[1]);

    //
    P = rings*M_PI/180.0;
    T = (t+sectors)*M_PI/180.0;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_spheres.push_back(R * x[2]);
    positions_spheres.push_back(R * y[2]);
    positions_spheres.push_back(R * z[2]);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(x[2]);
    normals_spheres.push_back(y[2]);
    normals_spheres.push_back(z[2]);
  }

  //Body of the sphere
  for (int p=rings; p<180-rings; p+=rings)
    for(int t=0; t<360; t+=sectors)
    {
      //A
      P = p*M_PI/180.0;
      T = t*M_PI/180.0;
      x[0] = sin(P) * cos(T) ;
      y[0] = sin(P) * sin(T) ;
      z[0] = cos(P);

      positions_spheres.push_back(R * x[0]);
      positions_spheres.push_back(R * y[0]);
      positions_spheres.push_back(R * z[0]);
      if (four_pos) positions_spheres.push_back(1.0);

      normals_spheres.push_back(x[0]);
      normals_spheres.push_back(y[0]);
      normals_spheres.push_back(z[0]);

      //B
      P = (p+rings)*M_PI/180.0;
      T = t*M_PI/180.0;
      x[1] = sin(P) * cos(T) ;
      y[1] = sin(P) * sin(T) ;
      z[1] = cos(P);
      positions_spheres.push_back(R * x[1]);
      positions_spheres.push_back(R * y[1]);
      positions_spheres.push_back(R * z[1]);
      if (four_pos) positions_spheres.push_back(1.0);

      normals_spheres.push_back(x[1]);
      normals_spheres.push_back(y[1]);
      normals_spheres.push_back(z[1]);

      //C
      P = p*M_PI/180.0;
      T = (t+sectors)*M_PI/180.0;
      x[2] = sin(P) * cos(T) ;
      y[2] = sin(P) * sin(T) ;
      z[2] = cos(P);
      positions_spheres.push_back(R * x[2]);
      positions_spheres.push_back(R * y[2]);
      positions_spheres.push_back(R * z[2]);
      if (four_pos) positions_spheres.push_back(1.0);


      normals_spheres.push_back(x[2]);
      normals_spheres.push_back(y[2]);
      normals_spheres.push_back(z[2]);
      //D
      P = (p+rings)*M_PI/180.0;
      T = (t+sectors)*M_PI/180.0;
      x[3] = sin(P) * cos(T) ;
      y[3] = sin(P) * sin(T) ;
      z[3] = cos(P);
      positions_spheres.push_back(R * x[3]);
      positions_spheres.push_back(R * y[3]);
      positions_spheres.push_back(R * z[3]);
      if (four_pos) positions_spheres.push_back(1.0);


      normals_spheres.push_back(x[3]);
      normals_spheres.push_back(y[3]);
      normals_spheres.push_back(z[3]);

      positions_spheres.push_back(R * x[1]);
      positions_spheres.push_back(R * y[1]);
      positions_spheres.push_back(R * z[1]);
      if (four_pos) positions_spheres.push_back(1.0);

      normals_spheres.push_back(x[1]);
      normals_spheres.push_back(y[1]);
      normals_spheres.push_back(z[1]);

      positions_spheres.push_back(R * x[2]);
      positions_spheres.push_back(R * y[2]);
      positions_spheres.push_back(R * z[2]);
      if (four_pos) positions_spheres.push_back(1.0);

      normals_spheres.push_back(x[2]);
      normals_spheres.push_back(y[2]);
      normals_spheres.push_back(z[2]);

    }
  //Bottom of the sphere
  for(int t=0; t<360; t+=sectors)
  {
    positions_spheres.push_back(0);
    positions_spheres.push_back(0);
    positions_spheres.push_back(-R);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(0);
    normals_spheres.push_back(0);
    normals_spheres.push_back(-1);


    P = (180-rings)*M_PI/180.0;
    T = t*M_PI/180.0;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);
    positions_spheres.push_back(R * x[1]);
    positions_spheres.push_back(R * y[1]);
    positions_spheres.push_back(R * z[1]);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(x[1]);
    normals_spheres.push_back(y[1]);
    normals_spheres.push_back(z[1]);


    P = (180-rings)*M_PI/180.0;
    T = (t+sectors)*M_PI/180.0;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_spheres.push_back(R * x[2]);
    positions_spheres.push_back(R * y[2]);
    positions_spheres.push_back(R * z[2]);
    if (four_pos) positions_spheres.push_back(1.0);

    normals_spheres.push_back(x[2]);
    normals_spheres.push_back(y[2]);
    normals_spheres.push_back(z[2]);

  }
}

template <class FLOAT>
void create_flat_and_wire_sphere(double R,
                   std::vector<FLOAT>& positions_spheres,
                   std::vector<FLOAT>& normals_spheres,
                   std::vector<FLOAT>& positions_wire_spheres,
                   bool four_pos=false)
{
  //The more small they are, the more precise the Sphere will be.
  // Must be a multiple of 360 and 180.
  const int rings=18;
  const int sectors=38;

  create_flat_sphere(R, positions_spheres, normals_spheres, four_pos);
  
  float T, P;
  float x[4],y[4],z[4];

  //Top of the sphere
  for(int t=0; t<360; t+=sectors)
  {
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(R);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    P = rings*M_PI/180.0;
    T = t*M_PI/180.0;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    //
    P = rings*M_PI/180.0;
    T = (t+sectors)*M_PI/180.0;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);
    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(R);
    if (four_pos) positions_wire_spheres.push_back(1.0);
  }

  //Body of the sphere
  for (int p=rings; p<180-rings; p+=rings)
    for(int t=0; t<360; t+=sectors)
    {
      //A
      P = p*M_PI/180.0;
      T = t*M_PI/180.0;
      x[0] = sin(P) * cos(T) ;
      y[0] = sin(P) * sin(T) ;
      z[0] = cos(P);

      //B
      P = (p+rings)*M_PI/180.0;
      T = t*M_PI/180.0;
      x[1] = sin(P) * cos(T) ;
      y[1] = sin(P) * sin(T) ;
      z[1] = cos(P);

      //C
      P = p*M_PI/180.0;
      T = (t+sectors)*M_PI/180.0;
      x[2] = sin(P) * cos(T) ;
      y[2] = sin(P) * sin(T) ;
      z[2] = cos(P);

      //D
      P = (p+rings)*M_PI/180.0;
      T = (t+sectors)*M_PI/180.0;
      x[3] = sin(P) * cos(T) ;
      y[3] = sin(P) * sin(T) ;
      z[3] = cos(P);

      positions_wire_spheres.push_back(R * x[0]);
      positions_wire_spheres.push_back(R * y[0]);
      positions_wire_spheres.push_back(R * z[0]);
      if (four_pos) positions_wire_spheres.push_back(1.0);

      positions_wire_spheres.push_back(R * x[1]);
      positions_wire_spheres.push_back(R * y[1]);
      positions_wire_spheres.push_back(R * z[1]);
      if (four_pos) positions_wire_spheres.push_back(1.0);

      positions_wire_spheres.push_back(R * x[1]);
      positions_wire_spheres.push_back(R * y[1]);
      positions_wire_spheres.push_back(R * z[1]);
      if (four_pos) positions_wire_spheres.push_back(1.0);


      positions_wire_spheres.push_back(R * x[3]);
      positions_wire_spheres.push_back(R * y[3]);
      positions_wire_spheres.push_back(R * z[3]);
      if (four_pos) positions_wire_spheres.push_back(1.0);


      positions_wire_spheres.push_back(R * x[3]);
      positions_wire_spheres.push_back(R * y[3]);
      positions_wire_spheres.push_back(R * z[3]);
      if (four_pos) positions_wire_spheres.push_back(1.0);


      positions_wire_spheres.push_back(R * x[2]);
      positions_wire_spheres.push_back(R * y[2]);
      positions_wire_spheres.push_back(R * z[2]);
      if (four_pos) positions_wire_spheres.push_back(1.0);

      positions_wire_spheres.push_back(R * x[2]);
      positions_wire_spheres.push_back(R * y[2]);
      positions_wire_spheres.push_back(R * z[2]);
      if (four_pos) positions_wire_spheres.push_back(1.0);

      positions_wire_spheres.push_back(R * x[0]);
      positions_wire_spheres.push_back(R * y[0]);
      positions_wire_spheres.push_back(R * z[0]);
      if (four_pos) positions_wire_spheres.push_back(1.0);
    }
  //Bottom of the sphere
  for(int t=0; t<360; t+=sectors)
  {
    P = (180-rings)*M_PI/180.0;
    T = t*M_PI/180.0;
    x[1] = sin(P) * cos(T) ;
    y[1] = sin(P) * sin(T) ;
    z[1] = cos(P);

    P = (180-rings)*M_PI/180.0;
    T = (t+sectors)*M_PI/180.0;
    x[2] = sin(P) * cos(T) ;
    y[2] = sin(P) * sin(T) ;
    z[2] = cos(P);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(-R);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);
    if (four_pos) positions_wire_spheres.push_back(1.0);


    positions_wire_spheres.push_back(R * x[1]);
    positions_wire_spheres.push_back(R * y[1]);
    positions_wire_spheres.push_back(R * z[1]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);
    if (four_pos) positions_wire_spheres.push_back(1.0);


    positions_wire_spheres.push_back(R * x[2]);
    positions_wire_spheres.push_back(R * y[2]);
    positions_wire_spheres.push_back(R * z[2]);
    if (four_pos) positions_wire_spheres.push_back(1.0);

    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(0);
    positions_wire_spheres.push_back(-R);
    if (four_pos) positions_wire_spheres.push_back(1.0);
  }
}

#endif //POLYHEDRON_DEMO_DRAW_SPHERE_H
