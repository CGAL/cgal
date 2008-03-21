// Define kernel, points, vectors, ... for MFC classes
//
#ifndef _GyrovizKernel_
#define _GyrovizKernel_
#pragma once


// CGAL
#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>



// kernel
struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef K::FT FT;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Sphere_3 Sphere;
typedef K::Vector_2 Vector_2;
typedef K::Vector_3 Vector_3;


#endif // _GyrovizKernel_