// Copyright (c) 2000,2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Bernd Gaertner, Sven Schoenherr <sven@inf.ethz.ch>

// Each of the following operators is individually
// protected against multiple inclusion.

// Window_stream I/O operators
// ===========================

// Conic_2
// -------
#ifdef CGAL_CONIC_2_H
#ifndef CGAL_IO_WINDOW_STREAM_CONIC_2
#define CGAL_IO_WINDOW_STREAM_CONIC_2

template< class R >
CGAL::Window_stream&
operator << ( CGAL::Window_stream& ws, const CGAL::Conic_2<R>& c)
{
    // length of a pixel in window-coordinates
    double pixel = 1/ws.scale();

    // pixel dimensions of window
    int width  = (int)((ws.xmax() - ws.xmin()) * ws.scale()) + 1,
        height = (int)((ws.ymax() - ws.ymin()) * ws.scale()) + 1,
        dim    = std::max( width, height);

    // pixel coordinates, stored for faster output
    double *X = new double [2*dim];
    double *Y = new double [2*dim];

    // actual number of pixels to be drawn
    int pixels;

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
        for (double x = ws.xmin(); x <= ws.xmax(); x+=pixel) {
            double discr = (t*t-4.0*r*s)*(x*x) + (2.0*t*v-4.0*s*u)*x +
                             v*v - 4.0*s*w;
            if (discr >= 0.0) {
                double y1 = (-t*x - v - CGAL::sqrt(discr))/(2.0*s);
                double y2 = (-t*x - v + CGAL::sqrt(discr))/(2.0*s);
                X[pixels] = x; Y[pixels++] = y1;
                X[pixels] = x; Y[pixels++] = y2; } }
    else
        for (double x = ws.xmin(); x <= ws.xmax(); x+=pixel) {
            double denom = t*x + v;
            if (denom != 0.0) {
                double y = -(r*x*x + u*x + w)/denom;
                X[pixels] = x; Y[pixels++] = y; } }
    ws.draw_pixels (pixels, X, Y);

    // Phase II (drawing in y-direction)
    pixels = 0;
    // solve conic equation for x
    if (r != 0.0)
        for (double y = ws.ymin(); y <= ws.ymax(); y+=pixel) {
            double discr = (t*t-4.0*r*s)*(y*y) + (2.0*t*u-4.0*r*v)*y +
                             u*u - 4.0*r*w;
            if (discr >= 0.0) {
                double x1 = (-t*y - u - CGAL::sqrt(discr))/(2.0*r);
                double x2 = (-t*y - u + CGAL::sqrt(discr))/(2.0*r);
                X[pixels] = x1; Y[pixels++] = y;
                X[pixels] = x2; Y[pixels++] = y; } }
    else
        for (double y = ws.ymin(); y <= ws.ymax(); y+=pixel) {
            double denom = t*y + u;
            if (denom != 0.0) {
                double x = -(s*y*y + v*y + w)/denom;
                X[pixels] = x; Y[pixels++] = y; } }
    ws.draw_pixels (pixels, X, Y);

    // free memory
    delete[] Y;
    delete[] X;

    return( ws);
}

#endif // CGAL_IO_WINDOW_STREAM_CONIC_2
#endif // CGAL_CONIC_2_H

// ===== EOF ==================================================================
