// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/IO/Conic_2_Window_stream.h
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Conic_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: graphical output to `leda_window' for Conic_2 algo.
// ============================================================================

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
