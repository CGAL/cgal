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

#include <qpainter.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

CGAL_BEGIN_NAMESPACE

/*! \brief
 * outputs a curve arc to \c Qt_widget
 */
template <class CKvA_2>
Qt_widget& operator << (Qt_widget& ws, const CGALi::Arc_2< CKvA_2 >& arc) {
    
    typedef Curve_renderer_facade<CKvA_2> Facade;

    typedef CGAL::Twotuple< int > Coord_2;
    typedef std::vector< Coord_2 > Coord_vec_2;

    CGAL::Twotuple< Coord_2 > end_points;
    std::list<Coord_vec_2> points;
   
    Facade::setup(CGAL::Bbox_2(ws.x_min(), ws.y_min(), ws.x_max(), ws.y_max()),
            ws.width(), ws.height());
    
    Facade::instance().draw(arc, points, end_points);
    if(points.empty()) 
        return ws;
        
    QPainter *ppnt = &ws.get_painter();
    int height = ws.height();

   // std::cerr << ws.width() << " and " <<  ws.height() << "\n";

    typename std::list<Coord_vec_2>::const_iterator lit = points.begin();
    while(lit != points.end()) {

        const Coord_vec_2& vec = *lit;
        typename Coord_vec_2::const_iterator vit = vec.begin();
            
        if(vec.size() == 2) {
            ppnt->moveTo(vit->e0, height - vit->e1);
            vit++;
            ppnt->lineTo(vit->e0, height - vit->e1);
                
        } else {
            ppnt->moveTo(vit->e0, height - vit->e1);
            //std::cerr << "(" << vit->e0 << "; " << vit->e1 << "\n";
            while(vit != vec.end()) {
                ppnt->lineTo(vit->e0, height - vit->e1);
                vit++;
                //if(vit != vec.end())
                //std::cerr << "(" << vit->e0 << "; " << vit->e1 << "\n";
            }
        }
        lit++;
    }
        
    QPen old_pen = ppnt->pen();
    ppnt->setPen(QPen(Qt::NoPen)); // avoid drawing outlines
    // draw with the current brush attributes
    ppnt->drawEllipse(end_points.e0.e0-3,height-end_points.e0.e1-3, 6, 6);
    ppnt->drawEllipse(end_points.e1.e0-3,height-end_points.e1.e1-3, 6, 6);
    ppnt->setPen(old_pen);

    return ws;
}

/*! \brief
 *  outputs a curve point to \c Qt_widget
 */
template <class CKvA_2>
Qt_widget& operator << (Qt_widget& ws, const CGALi::Point_2< CKvA_2 >& pt) {
    
    typedef Curve_renderer_facade<CKvA_2> Facade;
   
    CGAL::Twotuple< int > coord;
    Facade::setup(CGAL::Bbox_2(ws.x_min(), ws.y_min(), ws.x_max(), ws.y_max()),
            ws.width(), ws.height());

    if(!Facade::instance().draw(pt, coord))
        return ws;
       
    QPainter *ppnt = &ws.get_painter();
    QPen old_pen = ppnt->pen();
    ppnt->setPen(QPen(Qt::NoPen));
    ppnt->drawEllipse(coord.x - 3, ws.height() - coord.y - 3, 6, 6);
    ppnt->setPen(old_pen);
    return ws;
}

CGAL_END_NAMESPACE

#endif // CGAL_QT_WIDGET_CURVE_RENDERER_2_H

