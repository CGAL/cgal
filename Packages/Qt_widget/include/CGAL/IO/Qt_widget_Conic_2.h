// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_Conic_2.h
// package       : Qt_widget
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================


#ifndef CGAL_QT_WIDGET_CONIC_2_H
#define CGAL_QT_WIDGET_CONIC_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Conic_2.h>

namespace CGAL{

template< class R >
Qt_widget&
operator << ( Qt_widget& ws, const CGAL::Conic_2<R>& c)
{
    // length of a pixel in window-coordinates
    double pixel_x = 1/ws.x_scal();
    double pixel_y = 1/ws.y_scal();

    // pixel dimensions of window
    int width  = (int)((ws.x_max() - ws.x_min()) * ws.x_scal()) + 1,
        height = (int)((ws.y_max() - ws.y_min()) * ws.y_scal()) + 1,
        dim    = std::max( width, height);

    // pixel coordinates, stored for faster output
    double *X = new double [2*dim];
    double *Y = new double [2*dim];

    // actual number of pixels to be drawn
    int pixels, ind;

    // conic coordinates
    double r = CGAL::to_double (c.r()),
           s = CGAL::to_double (c.s()),
           t = CGAL::to_double (c.t()),
           u = CGAL::to_double (c.u()),
           v = CGAL::to_double (c.v()),
           w = CGAL::to_double (c.w());

    // Phase I (drawing in x-direction)
    pixels = 0;
    // solve conic equation for y
    if (s != 0.0)
        for (double x = ws.x_min(); x <= ws.x_max(); x+=pixel_x) {
            double discr = (t*t-4.0*r*s)*(x*x) + (2.0*t*v-4.0*s*u)*x +
                             v*v - 4.0*s*w;
            if (discr >= 0.0) {
                double y1 = (-t*x - v - CGAL::sqrt(discr))/(2.0*s);
                double y2 = (-t*x - v + CGAL::sqrt(discr))/(2.0*s);
                X[pixels] = x; Y[pixels++] = y1;
                X[pixels] = x; Y[pixels++] = y2; } }
    else
        for (double x = ws.x_min(); x <= ws.x_max(); x+=pixel_x) {
            double denom = t*x + v;
            if (denom != 0.0) {
                double y = -(r*x*x + u*x + w)/denom;
                X[pixels] = x; Y[pixels++] = y; } }
    //ws.draw_pixels (pixels, X, Y);
    for (ind = 0; ind < pixels; ind++)
    {
      typedef Cartesian<double> Rep;
      typedef Point_2<Rep>	Point;
      ws << Point(X[ind], Y[ind]);
    }

    // Phase II (drawing in y-direction)
    pixels = 0;
    // solve conic equation for x
    if (r != 0.0)
        for (double y = ws.y_min(); y <= ws.y_max(); y+=pixel_y) {
            double discr = (t*t-4.0*r*s)*(y*y) + (2.0*t*u-4.0*r*v)*y +
                             u*u - 4.0*r*w;
            if (discr >= 0.0) {
                double x1 = (-t*y - u - CGAL::sqrt(discr))/(2.0*r);
                double x2 = (-t*y - u + CGAL::sqrt(discr))/(2.0*r);
                X[pixels] = x1; Y[pixels++] = y;
                X[pixels] = x2; Y[pixels++] = y; } }
    else
        for (double y = ws.y_min(); y <= ws.y_max(); y+=pixel_y) {
            double denom = t*y + u;
            if (denom != 0.0) {
                double x = -(s*y*y + v*y + w)/denom;
                X[pixels] = x; Y[pixels++] = y; } }
    //ws.draw_pixels (pixels, X, Y);
    for (ind = 0; ind < pixels; ind++)
    {
      typedef Cartesian<double> Rep;
      typedef Point_2<Rep>	Point;
      ws << Point(X[ind], Y[ind]);
    }


    // free memory
    delete[] Y;
    delete[] X;

    return( ws);
}

}//end namespace CGAL

#endif
