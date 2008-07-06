// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*!\file CGAL/IO/Qt_widget_Curve_renderer_2.h
 * \brief
  * provides \c CGAL::Qt_widget interface for the curve renderer
 */

#ifndef CGAL_QT_WIDGET_CURVE_RENDERER_2_H
#define CGAL_QT_WIDGET_CURVE_RENDERER_2_H

#ifndef CGAL_CKVA_NO_QT_WIDGET_INTERFACE
#include <qpainter.h>
#include <CGAL/IO/Qt_widget.h>
#endif

// #include <CGAL/Interval_traits.h>
// #include <CGAL/Bigfloat_interval_traits.h>
// #include <CGAL/convert_to_bfi.h>
// #include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief 
 * represents a single curve renderer instance to speed up rendering of several
 * objects using the same parameter set
 * 
 * @warning not recommended to use for multi-threaded applications
 */
template <class CurvedKernelViaAnalysis_2>
class Curve_renderer_singleton
{
    Curve_renderer_singleton() { // private constructor
    }

public:
    typedef CGAL::Interval_nt<true> Interval_double;

    typedef Curve_renderer_facade<CurvedKernelViaAnalysis_2, Interval_double>
         Curve_renderer_inst;

    typedef typename Curve_renderer_inst::Coord_vec_2 Coord_vec_2;

    typedef typename Curve_renderer_inst::Coord_2 Coord_2;
    
    inline static void draw(const CGAL::Bbox_2& bbox, int res_w, int res_h,
        const typename Curve_renderer_inst::Arc_2& arc,
        std::list<Coord_vec_2>& points, 
        std::pair<Coord_2, Coord_2>& endpoints) {
        
        _instance()._setup(bbox, res_w, res_h);
        _instance().gfx_instance.draw(arc, points, endpoints);
    }

    //! returns \c false if the point does not fall within drawing box
    inline static bool draw(const CGAL::Bbox_2& bbox, int res_w, int res_h,
        const typename Curve_renderer_inst::Point_2& pt, Coord_2& coord) {
        
        _instance()._setup(bbox, res_w, res_h);
        return _instance().gfx_instance.draw(pt, coord);
    }

#ifndef CGAL_CKVA_NO_QT_WIDGET_INTERFACE
    static void draw_qt(Qt_widget& ws, 
       const typename Curve_renderer_inst::Arc_2& arc)
    {
        std::list<Coord_vec_2> points;
        std::pair<Coord_2, Coord_2> end_points;

        draw(CGAL::Bbox_2(ws.x_min(), ws.y_min(), ws.x_max(), ws.y_max()),
            ws.width(), ws.height(), arc, points, end_points);
        if(points.empty()) 
            return;
        
        QPainter *ppnt = &ws.get_painter();
        int height = ws.height();

        std::list<Coord_vec_2>::const_iterator lit = points.begin();
        while(lit != points.end()) {

            const Coord_vec_2& vec = *lit;
            Coord_vec_2::const_iterator vit = vec.begin();
            
            if(vec.size() == 2) {
                ppnt->moveTo(vit->x, height - vit->y);
                vit++;
                ppnt->lineTo(vit->x, height - vit->y);
                
            } else {
                ppnt->moveTo(vit->x, height - vit->y);
                //std::cout << "(" << (*it).x << "; " << height - (*it).y <<
                    //")\n";
                while(vit != vec.end()) {
                    ppnt->lineTo(vit->x, height - vit->y);
                //std::cout << "(" << (*it).x << "; " << height - (*it).y <<
                    //")\n";
                    vit++;
                }
            }
            lit++;
        }
        
        QPen old_pen = ppnt->pen();
        ppnt->setPen(QPen(Qt::NoPen)); // avoid drawing outlines
        // draw with the current brush attributes
        ppnt->drawEllipse(end_points.first.x-3,height-end_points.first.y-3,
              6, 6);
        ppnt->drawEllipse(end_points.second.x-3,height-end_points.second.y-3,
              6, 6);
        ppnt->setPen(old_pen);
    }
    
    static void draw_qt(Qt_widget& ws, 
        const typename Curve_renderer_inst::Point_2& pt)
    {
        Coord_2 coord;
        if(!draw(CGAL::Bbox_2(ws.x_min(), ws.y_min(), ws.x_max(), ws.y_max()),
                ws.width(), ws.height(), pt, coord))
            return;
       
        QPainter *ppnt = &ws.get_painter();
        int height = ws.height();
       
        QPen old_pen = ppnt->pen();
        ppnt->setPen(QPen(Qt::NoPen));
        ppnt->drawEllipse(coord.x - 3, height - coord.y - 3, 6, 6);
        ppnt->setPen(old_pen);
    }
#endif // !CGAL_CKVA_NO_QT_WIDGET_INTERFACE

private:

    static Curve_renderer_singleton& _instance() {
        static Curve_renderer_singleton _this;
        return _this;
    }

    void _setup(const CGAL::Bbox_2& bbox, int res_w, int res_h) {
        int _w, _h;
        CGAL::Bbox_2 tmp;
        gfx_instance.renderer().get_resolution(_w, _h);
        gfx_instance.renderer().get_window(tmp);

        if(bbox != tmp || res_w != _w || res_h != _h) {
            gfx_instance.setup(bbox, res_w, res_h);
        }
    }

    Curve_renderer_inst gfx_instance;    
}; 

} // namespace CGALi

#ifndef CGAL_CKVA_NO_QT_WIDGET_INTERFACE
/*! \brief
 * outputs a curve arc to \c Qt_widget
 */
template <class CKvA>
Qt_widget& operator << (Qt_widget& ws, const CGALi::Arc_2<CKvA>& arc) {
    
    CGALi::Curve_renderer_singleton<CKvA>::draw_qt(ws, arc);
    return ws;
}

/*! \brief
 *  outputs a curve point to \c Qt_widget
 */
template <class CKvA>
Qt_widget& operator << (Qt_widget& ws, const CGALi::Point_2<CKvA>& pt) {
    
    CGALi::Curve_renderer_singleton<CKvA>::draw_qt(ws, pt);
    return ws;
}
#endif // !CGAL_CKVA_NO_QT_WIDGET_INTERFACE

CGAL_END_NAMESPACE

#endif // CGAL_QT_WIDGET_CURVE_RENDERER_2_H

