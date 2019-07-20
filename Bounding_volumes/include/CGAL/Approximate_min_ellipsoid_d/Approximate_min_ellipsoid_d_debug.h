// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#ifndef CGAL_APPROX_MIN_ELL_D_DEBUG_H
#define CGAL_APPROX_MIN_ELL_D_DEBUG_H

#include <CGAL/license/Bounding_volumes.h>


#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/random/linear_congruential.hpp>

namespace CGAL {

  // We hide the following debugging utilities in a "private" namespace:
  namespace Approximate_min_ellipsoid_d_impl {

    // For debugging only:
    template<typename Iterator>
    void print_matrix(int d, const char *name, Iterator A)
    {
      std::cout << name << ":= Matrix([\n";
      for (int i=0; i<d; ++i) {
	std::cout << "  [ ";
	for (int j=0; j<d; ++j) {
	  std::cout << std::setprecision(30) << A[i+j*d];
	  if (j<d-1)
	    std::cout << ", ";
	}
	std::cout << " ]";
	if (i<d-1)
	  std::cout << ",";
	std::cout << "\n";
      }
      std::cout << "]);\n";	
    }

    // For debugging only:
    template<typename Iterator>
    void print_vector(int d, const char *name, Iterator v)
    {
      std::cout << name << ":= Vector([\n";
      for (int j=0; j<d; ++j) {
	std::cout << std::setprecision(30) << v[j];
	if (j<d-1)
	  std::cout << ", ";
      }
      std::cout << "]);\n";	
    }

    #ifdef CGAL_APPEL_LOG_MODE
    class Logger
    // A singleton class which sends debugging comments to log files.
    // (A class is called a "singleton" if at most one instance of it
    // may ever exist.)  By calling void log(channel,msg) you can send
    // the string msg to the file with name channel.
    {
    private: // private members:
      typedef std::map<std::string,std::ofstream*> Streams;
      Streams channels;           // a collection of pairs (k,v) where
                                  // k is the file-name and v is the
                                  // (open) stream associated with k
  
    private: // (Construction and destruction are private to prevent
  	   // more than one instantiation.)
  
      Logger() {}
      Logger(const Logger&);
      Logger& operator=(const Logger&);
  
      ~Logger()
      {
        // we walk through the list of all opened files and close them:
        for (Streams::iterator it = channels.begin();
  	   it != channels.end(); ++it) {
  	(*it).second->close();
  	delete (*it).second;
        }
      }
      
    public: // access and routines:
  
      static Logger& instance()
      // Returns a reference to the only existing instance of this class:
      {
        // Here's where we maintain the only instance: (Notice that it
        // gets constructed automatically the first time instance() is
        // called, and that it gets disposed of (if ever contructed) at
        // program termination.)
        static Logger instance;
        return instance;
      }
  
      void log(const char *channel,const std::string& msg)
      // If this is the first call to log with string channel as the
      // first parameter, then the file with name channel.log is
      // opened (at the beginning) for writing and msg is written to
      // it.  Otherwise, the string msg is opened to the already open
      // file channel.log.
      {
        const std::string name(channel);
        Streams::iterator it = channels.find(name);
        
        // have we already opened this file?
        if (it != channels.end()) {
  	// If so, then just append the message:
  	*(*it).second << msg;
  	(*it).second->flush();
  	
        } else {
  	// If we haven't seen 'name' before, we create a new file:
  	using std::ofstream;
  	ofstream *o = new ofstream((name+".log").c_str(),
  				   ofstream::out|ofstream::trunc);
  	channels[name] = o;
  	*o << msg;
        }
      }
    };
    #endif // CGAL_APPEL_LOG_MODE

  } // end of namespace Approximate_min_ellipsoid_d_impl
} // end of namespace CGAL

#ifdef CGAL_APPEL_TIMER_MODE
#include <sys/time.h>
#include <sys/resource.h>

namespace CGAL {
  namespace Approximate_min_ellipsoid_d_impl {

    // The following routine is taken from file mptimeval.h from
    // "Matpack Library Release 1.7.1" which is copyright (C) 1991-2002
    // by Berndt M. Gammel.  It works on the timeval struct defined in
    // sys/time.h:
    //
    //       struct timeval {
    //         long tv_sec;        /* seconds */
    //         long tv_usec;       /* microseconds */
    //       };
    //
    inline timeval& operator-=(timeval &t1,const timeval &t2)
    {
      t1.tv_sec -= t2.tv_sec;
      if ( (t1.tv_usec -= t2.tv_usec) < 0 ) {
        --t1.tv_sec;
        t1.tv_usec += 1000000;
      }
      return t1;
    }

    // (Adapted from the above.)
    inline timeval& operator+=(timeval &t1,const timeval &t2)
    {
      t1.tv_sec += t2.tv_sec;
      if ( (t1.tv_usec += t2.tv_usec) > 1000000) {
        ++t1.tv_sec;
        t1.tv_usec -= 1000000;
      }
      return t1;
    }
    
    class Timer
      // A singleton class which maintains a collection of named timers.
      // (A class is called a "singleton" if at most one instance of it
      // may ever exist.)  The following routines are provided:
      //
      // - start(name): If this is the first time start() has been
      //   called with name as the first parameter, then a new timer is
      //   created and started.  Otherwise, the timer with name name is
      //   restarted.
      //
      // - lapse(name): Retuns the number of seconds which have elapsed
      //   since start(name) was called last.
      //   Precondition: start(name) has been called once.
    {
    private: // (Construction and destruction are private to prevent
  	   // more than one instantiation.)
  
      Timer() {}
      Timer(const Timer&);
      Timer& operator=(const Timer&);
  
    public: // access and routines:
  
      static Timer& instance()
      // Returns a reference to the only existing instance of this class:
      {
        // Here's where we maintain the only instance: (Notice that it
        // gets constructed automatically the first time instance() is
        // called, and that it gets disposed of (if ever contructed) at
        // program termination.)
        static Timer instance;
        return instance;
      }
  
      void start(const char *timer_name)
      {
        // fetch current usage:
        rusage now;
        int status = getrusage(RUSAGE_SELF,&now);
        CGAL_assertion(status == 0);
        
        // save it:
        timers[std::string(timer_name)] = now.ru_utime;
      }
  
      float lapse(const char *name)
      {
        // assert that start(name) has been called before:
        CGAL_assertion(timers.find(std::string(name)) != timers.end());
        
        // get current usage:
        rusage now;
        int status = getrusage(RUSAGE_SELF,&now);
        CGAL_assertion(status == 0);
        
        // compute elapsed usage:
        now.ru_utime -= (*timers.find(std::string(name))).second;
        return now.ru_utime.tv_sec + now.ru_utime.tv_usec * 1e-6;
      }
      
    private: // private members:
      typedef std::map<std::string,timeval> Timers;
      Timers timers;              // a collection of pairs (k,v) where
                                  // k is the timer name and v is the
                                  // (started) timer associated with k
    };
  } // end of namespace Approximate_min_ellipsoid_d_impl
} // end of namespace CGAL
#endif // CGAL_APPEL_TIMER_MODE

namespace CGAL {
  namespace Approximate_min_ellipsoid_d_impl {

    template <typename T>
    std::string tostr(const T& t) {
      std::stringstream strm;
      strm << t;
      return strm.str();
    }
    
    class Eps_export_2 {
    // An instance of the following class accepts circles and ellipses
    // and procudes an Enhanced-PostScript figure.
    public:
      enum Stroke_mode { Solid=0, Solid_filled=1, Dashed=2 };
      enum Label_mode { None, Angle, Random_angle };
      
    private:
      std::vector<double> bb;            // bounding box
      Stroke_mode bm;                    // current stroke mode
      Label_mode lm;                     // current label mode
      double next_angle;                 // next angle to use (for Angle mode)
      int count;                         // counts the objects (starting at 1)
      std::string next_label;            // next label to use
      std::ostringstream body;           // buffered output
      bool adjust_bb;
      double zoom;
      std::string filename;
      boost::rand48 rng;  

    public: // construction and destruction:
      
      // Begins exporting to file filename.  Sets the current stroke mode
      // to Solid and the current label-mode to Random_angle.
      Eps_export_2(std::string filename, double zoom, int /* seed */ = 0) 
	: bb(4,0.0), bm(Solid), lm(Random_angle), count(1), 
	  next_label(default_label(1)), adjust_bb(true),
	  zoom(zoom), filename(filename)
      {}
      
      // Ends exporting and closes the file.
      ~Eps_export_2()
      {
	// open output file:
	using std::endl;
	std::ofstream f(filename.c_str());
	if (!f)
	  std::cerr << "Couldn't open file " << filename << "." << endl;
	
	// write header:
	const double margin = 10.0;
	f << "%!PS-Adobe-2.0 EPSF-2.0" << endl
	  << "%%Title: Some balls" << endl
	  << "%%Creator: a program by kf." << endl
	  << "%%CreationDate: now" << endl
	  << "%%For: you, the unknown" << endl
	  << "%%Pages: 1" << endl
	  << "%%DocumentFonts: Helvetica" << endl
	  << "%%BoundingBox: "
	  << bb[0]*zoom-margin << " "
	  << bb[1]*zoom-margin << " "
	  << bb[2]*zoom+margin << " "
	  << bb[3]*zoom+margin << endl
	  << "%%Magnification: 1.0000" << endl
	  << "%%EndComments" << endl
	  << "%%BeginProlog" << endl
	  << "%" << endl
	  << "/zoom " << zoom << " def" << endl
	  << "zoom zoom scale" << endl
	  << "/Helvetica findfont 7 zoom div scalefont setfont" << endl
	  << "0.06299 7.500 mul zoom div setlinewidth" << endl
	  << "/dashlen 4.5 zoom div def" << endl
	  << "/dashlen2 1.5 zoom div def" << endl
	  << "/ptlen 3.5 zoom div def" << endl << "%" << endl
	  << "% center-x center-y radius bp" << endl
	  << "% creates a path representing a circle" << endl
	  << "/bp { /radius exch def /cy exch def /cx exch def" << endl
	  << "radius 0 eq" << endl
	  << "{ newpath cx ptlen sub cy moveto ptlen ptlen add 0" << endl
	  << "rlineto cx cy ptlen sub moveto 0 ptlen ptlen add" << endl
	  << "rlineto }" << endl
	  << "{ newpath cx radius add cy moveto cx cy radius" << endl
	  << "0 360 arc closepath } ifelse } def" << "%" << endl
	  << "%center-x center-y radius ball" << endl
	  << "% draws a circle, unfilled, no dash" << endl
	  << "/ball { bp stroke } def" << endl
	  << "% center-x center-y radius bball" << endl
	  << "% draws a circle, filled, no dash" << endl
	  << "/bball { /radius exch def /cy exch def /cx exch def" << endl
	  << "cx cy radius radius 0 eq" << endl
	  << "{ gsave bp 0.6 setgray stroke" << endl
	  << "newpath cx cy moveto closepath 1 setlinecap" << endl
	  << "0.0 setgray stroke grestore }" << endl
	  << "{ bp gsave 0.6 setgray fill grestore stroke }" << endl
	  << "ifelse } def" << endl
	  << "% center-x center-y radius dball" << endl
	  << "% draws a circle, unfilled, dashed" << endl
	  << "/dball { gsave bp [dashlen] 0 setdash stroke grestore" << endl
	  << "} def" << endl
	  << "/eball { gsave bp 0.6 setgray stroke grestore } def" << endl
	  << "%" << endl
	  << "% center-x center-y radius label dist angle drawl" << endl
	  << "% Draws the label label." << endl
	  << "/drawl { /angle exch def /dist exch def /label exch def" << endl
	  << "/radius exch def /cy exch def /cx exch def" << endl
	  << "cx radius dist zoom div add angle cos mul add" << endl
	  << "cy radius dist zoom div add angle sin mul add" << endl
	  << "moveto label show } def" << endl
	  << "%" << endl
	  << "% Usage a b c find-roots x1 x2" << endl
	  << "% Finds the two roots of a x^2 + b x + c = 0, assuming" << endl
	  << "% that a is nonzero." << endl
	  << "/find-roots {" << endl
	  << "    % compute discriminant:" << endl
	  << "    3 copy 3 2 roll" << endl
	  << "    mul 4.0 mul neg exch dup mul add" << endl
	  << "    % stack: a b c disc" << endl
	  << "    sqrt" << endl
	  << "    2 index 0 ge { neg } { } ifelse" << endl
	  << "    % stack: a b c sqrt(disc)" << endl
	  << "    % compute first solution (sqrt(disc)-b) / (2*a):" << endl
	  << "    dup 3 index sub 4 index 2.0 mul div" << endl
	  << "    5 1 roll" << endl
	  << "    % stack: x1 a b c sqrt(disc)" << endl
	  << "    % compute second solution 2*c / (sqrt(disc)-b):" << endl
	  << "    3 2 roll sub exch 2.0 mul exch div" << endl
	  << "    exch pop" << endl
	  << "} def" << endl
	  << "%" << endl
	  << "% Usage: u v normalize u' v'" << endl
	  << "% Takes the vector [u,v] and normalizes it to length 1." << endl
	  << "/normalize {" << endl
	  << "    dup dup mul" << endl
	  << "    % stack: u v v^2" << endl
	  << "    2 index dup mul add" << endl
	  << "    % stack: u v v^2+u^2" << endl
	  << "    sqrt 1.0 exch div dup 4 1 roll mul 3 1 roll mul exch" << endl
	  << "} def" << endl
	  << "%" << endl
	  << "% Usage: a b h cellipse" << endl
	  << "% Draws the ellipse x^T M x <= 1 with M = [[a,h],[h,b]]." << endl
	  << "%" << endl
	  << "/cellipse {" << endl
	  << "% compute Eigen decomposition of M:" << endl
	  << "%" << endl
	  << "% stack: a b h" << endl
	  << "dup 0.0 eq {" << endl
	  << "	% ellipse is a sphere:" << endl
	  << "	2 index sqrt" << endl
	  << "	2 index sqrt" << endl
	  << "	[ 1.0 0.0 0.0 1.0 0.0 0.0 ]" << endl
	  << "}" << endl
	  << "{" << endl
	  << "	1.0" << endl
	  << "	4 copy pop pop add neg" << endl
	  << "	5 copy pop pop dup mul neg 3 1 roll mul add" << endl
	  << "	find-roots" << endl
	  << "	% stack: a b h x1 x2" << endl
	  << "	%" << endl
	  << "	% build first Eigenvector:" << endl
	  << "	2 index 2 index 6 index sub normalize" << endl
	  << "	% build second Eigenvector:" << endl
	  << "	% stack: a b h x1 x2 ev1x ev1y" << endl
	  << "	4 index 3 index 8 index sub normalize" << endl
	  << "	% build matrix containing Eigenvectors:" << endl
	  << "	% stack: a b h x1 x2 ev1x ev1y ev2x ev2y" << endl
	  << "	6 array" << endl
	  << "	dup 0 6 index put" << endl
	  << "	dup 1 4 index put" << endl
	  << "	dup 2 5 index put" << endl
	  << "	dup 3 3 index put" << endl
	  << "	dup 4 0 put" << endl
	  << "	dup 5 0 put" << endl
	  << "	5 1 roll pop pop pop pop" << endl
	  << "} ifelse" << endl
	  << "    % stack:  a b h x1 x2 T" << endl
	  << "    %" << endl
	  << "    % We now draw an circle scaled by x1 (along the " << endl
	  << "    % x-axis) and x2 (along the y-axis):" << endl
	  << "    matrix currentmatrix    %  remember CTM" << endl
	  << "    % stack:  a b h x1 x2 T old-ctm" << endl
	  << "    exch matrix invertmatrix concat" << endl
	  << "    2 index sqrt 1.0 exch div 2 index sqrt 1.0" << endl
	  << "    exch div scale" << endl
	  << "    newpath" << endl
	  << "    0 0 1 0 360 arc closepath" << endl
	  << "    setmatrix               %  restore CTM" << endl
	  << "    pop pop pop pop pop" << endl
	  << "} def" << endl
	  << "%" << endl
	  << "% Usage: a b d u v mu ellipse" << endl
	  << "% Finds the path (i.e., border) of the ellipse" << endl
	  << "%" << endl
	  << "%   E = { x | x^T M x + x^T m + mu <= 0 }," << endl
	  << "%" << endl
	  << "% where" << endl
	  << "%" << endl
	  << "%       [ a  b ]        [ u ]" << endl
	  << "%   M = [      ],   m = [   ]" << endl
	  << "%       [ b  d ]        [ v ]" << endl
	  << "%" << endl
	  << "% and det(M) > 0." << endl
	  << "/ellipse {" << endl
	  << "  /emu exch def" << endl
	  << "  /ev exch def" << endl
	  << "  /eu exch def" << endl
	  << "  /ed exch def" << endl
	  << "  /eb exch def" << endl
	  << "  /ea exch def" << endl
	  << "  /ematrix [ ea eb eb ed 0 0 ] matrix invertmatrix def" << endl
	  << "  % compute z = M^{-1} m:" << endl
	  << "  eu ev ematrix transform" << endl
	  << "  /ezy exch def" << endl
	  << "  /ezx exch def" << endl
	  << "  % translate to center:" << endl
	  << "  matrix currentmatrix    %  remember CTM" << endl
	  << "  ezx neg 2 div ezy neg 2 div translate" << endl
	  << "  % compute matrix of (now centrally symmetric) ellipse:" << endl
	  << "  /efactor ezx eu mul ezy ev mul add 4.0 div emu sub def" << endl
	  << "  ea efactor div ed efactor div eb efactor div cellipse" << endl
	  << "  setmatrix               % restore CTM" << endl
	  << "} def" << endl
	  << "/setcol { 193 255 div 152 255 div 159 255 div" << endl
	  << "          setrgbcolor" << endl
	  << "} def" << endl
	  << "%" << endl
	  << "%" << endl
	  << "%%EndProlog" << endl
	  << "%%Page: 1 1" << endl;
	  
	// write content:
	f << "gsave" << endl
	  << body.str()
	  << "grestore" << endl
	  << "showpage" << endl
	  << "%%Trailer" << endl;
	
	// write footer:
	f.close();
      }
      
    public: // modifiers:
      // Sets the stroke mode for the next objects.
      void set_stroke_mode(Stroke_mode m)
      {
	bm = m;
      }
      
      // Sets the labelmode for the next objects.
      void set_label_mode(Label_mode m)
      {
	lm = m;
      }
      
      // Sets the label for the next object.  (This only makes sense if
      // the label mode for the next object will be different from None.)
      void set_label(const std::string& l)
      {
	next_label = l;
      }
      
      // Sets the angle for the next object drawn in mode Angle or
      // Random_angle.  (So you can temporarily overwrite the random
      // angle.)
      void set_angle(double angle)
      {
	next_angle = angle;
      }
      
      // Enable or disable automatic adjusting of the bounding box.
      void enlarge_bounding_box(bool active)
      {
	adjust_bb = active;
      }
      
      // Renders the line from (x,y) to (u,v) in the current stroke
      // and label mode.
      //
      // Implementation note: stroke and label mode are not yet implemented.
      void write_line(double x, double y, double u, double v,
		      bool colored = false)
      {
	// set color:
	if (colored)
	  body << "gsave setcol" << std::endl;
	
	// output ball in appropriate mode:
	adjust_bounding_box(x,y);
	adjust_bounding_box(u,v);
	body << "newpath\n"
	     << x << ' ' << y << " moveto\n"
	     << u << ' ' << v << " lineto\n"
	     << "stroke\n";

	// restore color:
	if (colored)
	  body << "grestore" << std::endl;
      }
      
      // Renders the circle with center (x,y) and radius r in
      // the current stroke and label mode.  All coordinates must be of
      // type double.
      void write_circle(double x, double y, double r, bool colored = false)
      {
	// set color:
	if (colored)
	  body << "gsave setcol" << std::endl;
	
	// output ball in appropriate mode:
	adjust_bounding_box(x-r,y-r);
	adjust_bounding_box(x+r,y+r);
	body << x << ' ' << y << ' ' << r << ' ';
	const char *mode[] = { "ball", "bball", "dball" };
	body << mode[static_cast<int>(bm)] << std::endl;
	
	// output label:
	if (lm != None) {
	  const double dist = 5.0;
	  const double lx = x + (r+dist/zoom)*std::cos(next_angle);
	  const double ly = y + (r+dist/zoom)*std::sin(next_angle);
	  adjust_bounding_box(lx,ly);
	  body << x << ' ' << y << ' ' << r << " (" << next_label << ") "
	       << dist << ' ' << next_angle << " drawl" << std::endl;
	}
	
	// restore color:
	if (colored)
	  body << "grestore" << std::endl;

	// generate next angle (if needed):
	if (lm == Random_angle)
	  next_angle = static_cast<double>(rng());
	
	// update counter and label:
	next_label = default_label(++count);
      }
      
      // Renders the centrally symmetric ellipse x^T M x <= 1 with M =
      // [[a,h],[h,b]] in the current stroke and label mode.
      void write_cs_ellipse(double a,double b,double h)
      {
	using std::endl;

	// Adjust the bounding box: For this, we use the fact that the
	// ellipse has semi-axes of lengths (1/sqrt(c),1/sqrt(d)) where
	//
        //   (c,d) = find_roots(1.0,-(a+b),a b - h^2)
	//
	// (See Sven's thesis, code on page 87.)  So it it enclosed in
	// some (rotated) rectangle of side lengths 2*c,2*d, centered
	// at the origin.  Consequently, the circle of radius
	// r=sqrt(c^2+d^2) encloses the ellipse, and we use this to
	// find the bounding box:
	const std::pair<double,double> semi = find_roots(1.0,-(a+b),a*b-h*h);
	const double r = std::sqrt(1/semi.first+1/semi.second);
	adjust_bounding_box(-r,-r);
	adjust_bounding_box(r,r);

	// draw the ellipse:
	body << "gsave" << endl 
	     << "  " << a << ' ' << b << ' ' << h << " cellipse" << endl;
	if (bm == Solid_filled)
	  body << "  gsave 0.6 setgray eofill grestore" << endl;
	if (bm == Dashed)
	  body << "  [dashlen] 0 setdash" << endl;
	body << "  stroke" << endl << "grestore" << endl;
	
	// output label:
	// Todo: Not implemented yet.
	if (false && lm != None) {
	  const double distance = 5.0;
	  body << a << ' ' << b << ' ' << distance << ' ' << next_angle
	       << " cellipse-pt moveto (" << next_label << ") show" << endl;
	}
	
	// generate next angle (if needed):
	if (lm == Random_angle)
	  next_angle = static_cast<double>(rng());
	
	// update counter and label:
	next_label = default_label(++count);
      }

      // Renders the ellipse
      // 
      //    E = { x | x^T M x + x^T m + mu <= 0 }
      //
      // where
      //
      //       [ a  h ]        [ u ]
      //   M = [      ],   m = [   ]
      //       [ h  b ]        [ v ]
      //
      // and det(M) > 0 in the current stroke and label mode.
      void write_ellipse(double a,double b,double h,
			 double u,double v,double mu)
      {
	using std::endl;
	
	// adjust bounding box (see file boundingbox.mw and my file note):
	const double tmp1 = a*b-h*h;
	const double tmp2 = ((u*(u*b-2.0*v*h)+v*v*a)/tmp1-4.0*mu)/tmp1;
	const double hw   = 0.5*std::sqrt(b*tmp2),  vw = 0.5*std::sqrt(a*tmp2);
	const double hoff = -0.5*(u*b-v*h)/tmp1, voff = 0.5*(u*h-v*a)/tmp1;
	adjust_bounding_box((std::min)(hw+hoff,-hw+hoff),
			    (std::min)(vw+voff,-vw+voff));
	adjust_bounding_box((std::max)(hw+hoff,-hw+hoff),
			    (std::max)(vw+voff,-vw+voff));

	#if 0
	// draw bounding box:
	body << "newpath " << (std::min)(hw+hoff,-hw+hoff) << " "
	     << (std::min)(vw+voff,-vw+voff) << " moveto "
	     << (std::max)(hw+hoff,-hw+hoff) << " " 
	     << (std::min)(vw+voff,-vw+voff) << " lineto "
	     << (std::max)(hw+hoff,-hw+hoff) << " "
	     << (std::max)(vw+voff,-vw+voff) << " lineto "
	     << (std::min)(hw+hoff,-hw+hoff)  << " "
	     << (std::max)(vw+voff,-vw+voff) << " lineto closepath stroke\n";
	#endif

	// begin drawing:
	body << "gsave" << endl;

	// draw the ellipse:
	body << "  " << a << ' ' << h << ' ' << b
	     << ' ' << u << ' ' << v << ' ' << mu
	     << " ellipse" << endl;
	if (bm == Solid_filled)
	  body << "  gsave 0.6 setgray eofill grestore" << endl;
	if (bm == Dashed)
	  body << "  [dashlen] 0 setdash" << endl;
	body << "  stroke" << endl;

	// output label:
	// Todo: Not implemented yet.

	// end drawing:
	body << "grestore" << endl;

	// generate next angle (if needed):
	if (lm == Random_angle)
	  next_angle = static_cast<double>(rng());
	
	// update counter and label:
	next_label = default_label(++count);
      }

      void adjust_bounding_box(double x,double y)
      // Make sure the bounding box is large enough to contain the point (x,y).
      {
	if (adjust_bb) {
	  bb[0] = (std::min)(x,bb[0]);
	  bb[2] = (std::max)(x,bb[2]);
	  bb[1] = (std::min)(y,bb[1]);
	  bb[3] = (std::max)(y,bb[3]);
	}
      }
      
    private: // utilities:

      std::pair<double,double> find_roots(double a,double b,double c)
      // Assuming a is nonzero, finds the roots of a x^2 + b x + c = 0.
      {
	double sd = std::sqrt(b*b-4.0*a*c);
	if (b >= 0.0)
	  sd = -sd;
	return std::pair<double,double>( (sd-b)/(2*a), 2*c/(sd-b) );
      }
      
      static std::string default_label(int count)
      {
	return tostr(count);
      }
    };
    
  } // namespace Approximate_min_ellipsoid_d_impl
  
} // namespace CGAL

#endif // CGAL_APPROX_MIN_ELL_D_DEBUG_H
