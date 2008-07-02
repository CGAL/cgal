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

#include <vector>

#include <qpainter.h>
#include <CGAL/IO/Qt_widget.h>

#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>
#include <CGAL/convert_to_bfi.h>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

/*!\brief 
 * singleton providing drawing routines for \c CurvedKernelViaAnalysis_2
 * 
 * @warning not recommended to use for multi-threaded applications
 */
template <class CurvedKernelViaAnalysis_2>
class Renderer_instance
{
    Renderer_instance() { // private constructor
    }

public:
    typedef typename CurvedKernelViaAnalysis_2::Curve_kernel_2::Coefficient
        Coefficient;

    typedef CGAL::Interval_nt<true> Interval_double;

    typedef Curve_renderer_facade<CurvedKernelViaAnalysis_2, Interval_double>
         Curve_renderer_inst;

    static void resize_window(const Qt_widget& ws)
    {
        CGAL::Bbox_2 new_box(ws.x_min(), ws.y_min(), ws.x_max(), ws.y_max());
        if(bbox != new_box) {
            bbox = new_box;
            gfx_instance().setup(bbox, ws.width(), ws.height());
        }
    }
    
    static void draw(Qt_widget& ws, 
        const typename Curve_renderer_inst::Arc_2& arc)
    {
        resize_window(ws);
        QPainter *ppnt = &ws.get_painter();
        int height = ws.height();
        std::list<Coord_vec_2> points;
        std::pair<Coord_2, Coord_2> end_points;
        
        gfx_instance().draw(arc, points, end_points);
        if(points.empty()) 
            return;

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
    
    static void draw(Qt_widget& ws, 
        const typename Curve_renderer_inst::Point_2& pt)
    {
        resize_window(ws);
        QPainter *ppnt = &ws.get_painter();
        int height = ws.height();
        
        Coord_2 coord;
        if(!gfx_instance().draw(pt, coord))
            return;

        QPen old_pen = ppnt->pen();
        ppnt->setPen(QPen(Qt::NoPen));
        ppnt->drawEllipse(coord.x - 3, height - coord.y - 3, 6, 6);
        ppnt->setPen(old_pen);
    }

private:
    static CGAL::Bbox_2 bbox;

    inline static Curve_renderer_inst& gfx_instance() {
        static Curve_renderer_inst _instance;
        return _instance;
    }
};

template <class CKvA>
CGAL::Bbox_2 Renderer_instance<CKvA>::bbox(0.0, 0.0, 0.0, 0.0);

} // namespace CGALi

/*! \brief
 * outputs a curve arc to \c Qt_widget
 */
template <class CurvedKernelViaAnalysis_2>
Qt_widget& operator << (Qt_widget& ws, 
    const CGALi::Arc_2<CurvedKernelViaAnalysis_2>& arc) {
    
    CGALi::Renderer_instance<CurvedKernelViaAnalysis_2>::draw(ws, arc);
    return ws;
}

/*! \brief
 *  outputs a curve point to \c Qt_widget
 */
template <class CurvedKernelViaAnalysis_2>
Qt_widget& operator << (Qt_widget& ws, 
    const CGALi::Point_2<CurvedKernelViaAnalysis_2>& pt) {
    
    CGALi::Renderer_instance<CurvedKernelViaAnalysis_2>::draw(ws, pt);
    return ws;
}

CGAL_END_NAMESPACE

#endif // CGAL_QT_WIDGET_CURVE_RENDERER_2_H

