// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
//
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_DUMMY_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_DUMMY_H

#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Aff_transformation_2.h>

namespace CGAL {
    
    // this class might to be moved to hyperbolic traits
    template<class GT>
    class Construct_reflection
    {
        typedef typename GT::Point_2  Point_2;
        typedef typename GT::Circle_2 Circle_2;
        typedef typename GT::Vector_2 Vector_2;
        typedef typename GT::FT       FT;
        
        typedef typename GT::Compute_squared_distance_2     Compute_squared_distance_2;
        typedef typename GT::Construct_vector_2             Construct_vector_2;
        typedef typename GT::Construct_scaled_vector_2      Construct_scaled_vector_2;
        typedef typename GT::Construct_translated_point_2   Construct_translated_point_2;
        
        typedef CGAL::Regular_triangulation_filtered_traits_2<GT> Rtr_2;
        typedef typename Rtr_2::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
        typedef typename Rtr_2::Weighted_point_2                  Weighted_point_2;
        typedef typename Rtr_2::Bare_point                        Bare_point;
        
    public:
        Construct_reflection(const Circle_2& c = Circle_2(Point_2(0.0, 0.0), 1.)): _unit_circle(c)
        {
        }
        
        // exact arithmetic! very nice!
        Point_2 operator ()(const Point_2& p, const Point_2& q, const Point_2& r) const
        {
            typedef typename GT::Collinear_2 Collinear_2;
            if(Collinear_2()(p, q, _unit_circle.center())){
                assert(false);
            }
            Weighted_point_2 wp(p);
            Weighted_point_2 wq(q);
            Weighted_point_2 wo(_unit_circle.center(), _unit_circle.squared_radius());
            
            Bare_point  center = Construct_weighted_circumcenter_2()(wp, wo, wq);
            FT          radius = Compute_squared_distance_2()(p, center);
            FT          dist   = Compute_squared_distance_2()(r,center);
            Vector_2    vect   = Construct_vector_2()(center,r);
                        vect   = Construct_scaled_vector_2()(vect,radius/dist);
            Point_2     image  = Construct_translated_point_2()(center,vect);
            return image;
        }
        
    private:
        const Circle_2& _unit_circle;
    };
    
    template<class GT>
    typename GT::Point_2 apply_rotation(const typename GT::Point_2& p)
    {
        // Note Iordan: This is rotation by \pi/4 -- the second and third argument are the values of cos(\theta) and sin(\theta), respectively
        CGAL::Aff_transformation_2<GT> rotate(CGAL::ROTATION, std::sqrt(0.5), std::sqrt(0.5));
        
        return rotate(p);
    }
    
    
    

    
    
    template<class GT>
    void compute_dummy_points(  std::vector<typename GT::Point_2>& inner_points,
                                std::vector<typename GT::Point_2>& points_on_boundary,
                                std::vector<typename GT::Point_2>& points_on_vertex)
    {
        assert(inner_points.size()        == 0);
        assert(points_on_boundary.size()  == 0);
        assert(points_on_vertex.size()    == 0);
        
        typedef typename GT::Kernel  K;
        typedef typename GT::Point_2 Point_2;
        
        double phi = CGAL_PI / 8.;
        double psi = CGAL_PI / 3.;
        double rho = std::sqrt(cos(psi)*cos(psi) - sin(phi)*sin(phi));
        
        const Point_2 o(0.0, 0.0);
        const Point_2 a(cos(phi)*cos(phi + psi)/rho, sin(phi)*cos(phi + psi)/rho);
        const Point_2 b(a.x(), -a.y());
        
        inner_points.push_back(b);
        Point_2 c = Construct_reflection<K>()(a, b, o);
        Point_2 d = Construct_reflection<K>()(a, c, b);
        //inner_points.push_back(d);
        Point_2 e = Construct_reflection<K>()(d, c, a);
        Point_2 f = Construct_reflection<K>()(d, e, c);
        
        int size = inner_points.size();
        for(int i = 0; i < 7; i++) {
            for(int j = 0; j < size; j++) {
                inner_points.push_back(apply_rotation<K>(inner_points[i*size + j]));
            }
        }
        inner_points.push_back(o);
        
        Point_2 cr = c;
        Point_2 fr = f;
        for (int i = 0; i < 2; i++) {
            cr = apply_rotation<K>(cr);
            fr = apply_rotation<K>(fr);
        }

        points_on_boundary.push_back(cr);
        for(int i = 1; i < 4; i++) {
           points_on_boundary.push_back(apply_rotation<K>(points_on_boundary[i-1]));
        }
        

        points_on_vertex.push_back(apply_rotation<K>(fr));
    }
    
    /*

    template<class GT>
    void recursive_translate(Diametric_translations<GT> g,
                             std::vector<typename GT::Point_2>& points,
                             int depth,
                             int start,
                             int end) {
        
        if (depth > 0) {
            int my_start = points.size();
            int my_end   = start;
            
            for (int i = start; i <= end; i++){
                typename GT::Point_2 subject = points[i];
                points.push_back( g.a().DoAction(subject) );
                points.push_back( g.b().inverse().DoAction(subject) );
                points.push_back( g.c().DoAction(subject) );
                points.push_back( g.d().inverse().DoAction(subject) );
                points.push_back( g.a().inverse().DoAction(subject) );
                points.push_back( g.b().DoAction(subject) );
                points.push_back( g.c().inverse().DoAction(subject) );
                points.push_back( g.d().DoAction(subject) );
                my_end += 8;
            }
            
            recursive_translate(g, points, depth - 1, my_start, my_end);
        }
        
    }
    
    
    
    template<class GT>
    void recursive_translate(Diametric_translations<GT> g,
                             std::vector<typename GT::Point_2>& points,
                             typename GT::Point_2 subject,
                             int depth) {
        
        if (depth > 0) {
            int start = points.size();
            int end   = start + 8;
            
            // Add points in the order indicated by the group -- not necessary, but seems logical
            points.push_back( g.a().DoAction(subject) );
            points.push_back( g.b().inverse().DoAction(subject) );
            points.push_back( g.c().DoAction(subject) );
            points.push_back( g.d().inverse().DoAction(subject) );
            points.push_back( g.a().inverse().DoAction(subject) );
            points.push_back( g.b().DoAction(subject) );
            points.push_back( g.c().inverse().DoAction(subject) );
            points.push_back( g.d().DoAction(subject) );
            
            recursive_translate(g, points, depth - 1, start, end);
        }
        
    }
*/
    
    template < class Gt, class Tds >
    void Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
    insert_dummy_points(std::vector<typename Gt::Point_2>& all_points) {
        
        std::vector<typename Gt::Point_2> inner_points, points_on_boundary, points_on_vertex;
        compute_dummy_points<Gt>(inner_points, points_on_boundary, points_on_vertex);
        all_points = inner_points;
        
        
        for (int i = 0; i < points_on_boundary.size(); i++) {
            all_points.push_back(points_on_boundary[i]);
        }
        
        for (int i = 0; i < points_on_vertex.size(); i++) {
            all_points.push_back(points_on_vertex[i]);
        }
        

        
        std::cout << "All points length: " << all_points.size() << std::endl;
        
        }
        
    } // namespace CGAL
    
#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_DUMMY_H
    
