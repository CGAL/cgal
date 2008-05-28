#ifndef _APSS_ORACLE_H_
#define _APSS_ORACLE_H_

#include <CGAL/Surface_mesher/Null_oracle_visitor.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Surface_mesher/Sphere_oracle_3.h>

#include <queue>

#ifdef CGAL_SURFACE_MESHER_DEBUG_CLIPPED_SEGMENT
#  include <string>
#  include <sstream>
#endif

template <class GT, class Surface, class Point_creator = CGAL::Creator_uniform_3<typename GT::FT,
          typename GT::Point_3>,
          class Visitor = CGAL::Surface_mesher::Null_oracle_visitor> 
class ApssOracle
{
    // private types
    typedef ApssOracle<GT,Surface,Point_creator,Visitor> Self;

    typedef CGAL::Surface_mesher::Sphere_oracle_3<GT, Point_creator> Sphere_oracle;
    typedef typename GT::Point_3 Point;
    typedef typename GT::FT FT;
    typedef typename GT::Sphere_3 Sphere_3;
public:

    // Public types
    typedef GT Geom_traits;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Segment_3 Segment_3;
    typedef typename GT::Ray_3 Ray_3;
    typedef typename GT::Line_3 Line_3;
    typedef Surface Surface_3;


private:
    // Private members
    Visitor visitor; // a visitor that can modify a point, before returning it.

public:

    // Constructors
    ApssOracle (Visitor visitor_ = Visitor() )
        : visitor(visitor_)
    {
    }

    class Intersect_3 
    {
        Visitor visitor;
    public:
        Intersect_3(Visitor visitor) : visitor(visitor)
        {
        }
    
        CGAL::Object operator()(const Surface_3& surface, Segment_3 s) const
        // s is passed by value, because it is clipped below
        {
            typename GT::Construct_point_on_3 point_on = GT().construct_point_on_3_object();
    
            typename Sphere_oracle::Intersect_3 clip = Sphere_oracle().intersect_3_object();
    
            const Sphere_3& sphere = surface.bounding_sphere();
    
            Point_3 a = point_on(s, 0);
            Point_3 b = point_on(s, 1);
    
            // if both extremities are on the same side of the surface, return
            // no intersection
            if(surf_equation(surface, a) * surf_equation(surface, b) > 0)
                return CGAL::Object();
    
            // Code for surfaces with boundaries
            // First rescale segment s = [a, b]
            if( clip.clip_segment(sphere, a, b) )
                return intersect_clipped_segment(surface, a, b, surface.squared_error_bound());
            // else
            return CGAL::Object();
        } // end operator()(Surface_3, Segment_3)
    
        CGAL::Object operator()(const Surface_3& surface, const Ray_3& r) const
        {
            typename Sphere_oracle::Intersect_3 clip =
            Sphere_oracle().intersect_3_object();
    
            const Sphere_3& sphere = surface.bounding_sphere();
    
            Point_3 a, b;
            if(clip.clip_ray(sphere, r, a, b))
            {
                return intersect_clipped_segment(surface, a, b, surface.squared_error_bound());
            }
            // else
            return CGAL::Object();
        } // end operator()(Surface_3, Ray_3)
    
        CGAL::Object operator()(const Surface_3& surface, const Line_3& l) const
        {
            typename Sphere_oracle::Intersect_3 clip = Sphere_oracle().intersect_3_object();
    
            const Sphere_3& sphere = surface.bounding_sphere();
    
            Point_3 a, b;
            if(clip.clip_line(sphere, l, a, b))
            {
                return intersect_clipped_segment(surface, a, b, surface.squared_error_bound());
            }
            else
                return CGAL::Object();
        }; // end operator()(Surface_3, Line_3)
    
        // debug function
        static std::string debug_point(const Surface_3& surface,
                                        const Point& p) 
        {
            std::stringstream s;
            s << p << " (distance=" 
                << CGAL::sqrt(CGAL::squared_distance(p,
                                        surface.bounding_sphere().center()))
                << ", sign=" << surf_equation(surface, p)
                << ")";
            return s.str();
        }
    
        static CGAL::Sign surf_equation (Surface_3 surface,
                                        const Point& p) 
        {
            return CGAL::sign(surface(p));
        } // @TODO, @WARNING: we use x(), y() and z()
    
    private:
        // Private functions
        CGAL::Object intersect_clipped_segment(const Surface_3& surface,
                                        Point p1,
                                        Point p2,
                                        const FT& squared_distance_bound) const
        {
            typename GT::Compute_squared_distance_3 squared_distance = 
                GT().compute_squared_distance_3_object();
            typename GT::Construct_midpoint_3 midpoint =
                GT().construct_midpoint_3_object();
    
            CGAL::Sign sign_at_p1 = surf_equation(surface, p1);
            CGAL::Sign sign_at_p2 = surf_equation(surface, p2);
    
            if( sign_at_p1 == CGAL::ZERO )
            {
                visitor.new_point(p1);
                return make_object(p1);
            }
            if( sign_at_p2 == CGAL::ZERO )
            {
                visitor.new_point(p2);
                return make_object(p2);
            }
    
            // if both extremities are on the same side of the surface, return
            // no intersection
            if(sign_at_p1 * sign_at_p2 > 0)
            return CGAL::Object();
    
            while(true)
            {
                Point mid = midpoint(p1, p2);
                const CGAL::Sign sign_at_mid = surf_equation(surface, mid);
        
                if ( sign_at_mid == CGAL::ZERO || 
                     squared_distance(p1, p2) < squared_distance_bound )
                // If the two points are close, then we must decide
                {
                    visitor.new_point(mid);
                    return make_object(mid);
                }
        
                // Else we must go on
                if ( sign_at_p1 * sign_at_mid < 0 )
                {
                    p2 = mid;
                    sign_at_p2 = sign_at_mid;
                }
                else
                {
                    p1 = mid;
                    sign_at_p1 = sign_at_mid;
                }
            }
        } // end intersect_clipped_segment

    }; // end nested class Intersect_3

    class Construct_initial_points
    {
        const Self& oracle;
    public:
        Construct_initial_points(const Self& oracle) : oracle(oracle)
        {
        }
        
        /** Picks n points from the input point cloud, and project them onto the MLS surface.
        */
        template <typename OutputIteratorPoints>
        OutputIteratorPoints operator() (const Surface_3& surface, OutputIteratorPoints out, int n = 10000)
        {
            std::vector<unsigned int> ids(surface.m->points.size());
            for (unsigned int i=0 ; i<ids.size() ; ++i)
                ids[i] = i;
            
            if ((unsigned int)n >= surface.m->points.size())
            {
                n = surface.m->points.size()-1;
            }
    
            while (n>0)
            {
                double random_number = double(rand())/double(RAND_MAX); // random number in [0,1]
                unsigned int i = (unsigned int)(random_number*ids.size());
                unsigned int id = ids[i];
                ids[i] = ids.back();
                ids.pop_back();
                
                Point p = surface.m->points[id];
                surface.project(p);
                
                *out++= p;
                --n;
            }
            return out;
        }
    }; // end nested class Construct_initial_points

    Construct_initial_points construct_initial_points_object() const
    {
        return Construct_initial_points(*this);
    }

    Intersect_3 intersect_3_object() const
    {
        return Intersect_3(visitor);
    }

    bool is_in_volume(const Surface_3& surface, const Point& p) const
    {
        return Intersect_3::surf_equation(surface, p)<0.;
    }

};

template <typename FT>
FT approximate_sqrt(const FT x) {
  return FT (CGAL_NTS sqrt(CGAL_NTS to_double(x)));
}



#endif
