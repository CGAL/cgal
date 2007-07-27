#ifndef CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#define CGAL_ARRANGEMENT_OF_SPHERES_3_BASIC_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Simple_cartesian.h>

//#define CGAL_ARRANGEMENT_OF_SPHERES_3_USE_TEMPLATES





#define CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE CGAL_BEGIN_NAMESPACE namespace CGAL_AOS3_internal {
#define CGAL_AOS3_END_INTERNAL_NAMESPACE CGAL_END_NAMESPACE }

#define CGAL_AOS3_INTERNAL_NS CGAL::CGAL_AOS3_internal

#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>



#ifndef CGAL_AOS3_USE_TEMPLATES
#define CGAL_AOS3_TYPENAME

#define CGAL_AOS3_TEMPLATE
#define CGAL_AOS3_TARG
#define CGAL_AOS3_TRAITS public: typedef Arrangement_of_spheres_traits_3 Traits; private: friend class Semicolon_eater
 
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
typedef struct Geom_traits: Simple_cartesian<double>{} Arrangement_of_spheres_3_geom_traits;
CGAL_AOS3_END_INTERNAL_NAMESPACE

#include <CGAL/Arrangement_of_spheres_traits_3.h>





#else
#define CGAL_AOS3_TYPENAME typename
#define CGAL_AOS3_TEMPLATE template <class Traits_t>
#define CGAL_AOS3_TARG <Traits_t>
#define CGAL_AOS3_TRAITS public: typedef Traits_t Traits; private: friend class Semicolon_eater


#endif



#endif
