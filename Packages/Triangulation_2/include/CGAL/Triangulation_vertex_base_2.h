

#ifndef CGAL_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/Triangulation_short_names_2.h>

template < class GT >
class CGAL_Triangulation_vertex_base_2 {

public:
    typedef typename GT::Point Point;


  CGAL_Triangulation_vertex_base_2 ()
        : _f(NULL)
    {}
    
    CGAL_Triangulation_vertex_base_2(const Point & p)
        :  _p(p), _f(NULL)
    {}
    
    CGAL_Triangulation_vertex_base_2(const Point & p, void* f)
        :  _p(p), _f(f)
    {}

    inline void set_point(const Point & p)
    {
        _p = p;
    }
    
    inline void set_face(void* f)
    {
        _f = f;
    }

    inline Point point() const
    {
        return _p;
    }
    
    
    inline void* face() const
    {
        return _f;
    }
    


private:
        Point _p;
        void * _f;

};

#endif CGAL_TRIANGULATION_VERTEX_BASE_2_H
