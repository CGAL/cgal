// Copyright (c) 1999-2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri, Herve Bronnimann, Sylvain Pion

#include <CGAL/basic.h>

#ifdef CGAL_USE_GEOMVIEW

#include <sstream>
#include <csignal>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

#include <sys/types.h> // kill() on SunPro requires these 2 #includes.
#include <signal.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/binary_file_io.h>

namespace CGAL {

Geomview_stream::Geomview_stream(const Bbox_3 &bbox,
				 const char *machine,
				 const char *login)
    : bb(bbox), vertex_color(BLACK), edge_color(BLACK), face_color(BLACK),
      wired_flag(false), echo_flag(true), raw_flag(false),
      trace_flag(false), binary_flag(false),
      line_width(1)
{
    setup_geomview(machine, login);
    frame(bbox);
    pickplane(bbox);
    set_vertex_radius((bbox.xmax() - bbox.xmin())/100.0);
}

Geomview_stream::~Geomview_stream()
{
    kill(pid, SIGKILL);  // kills geomview
}

void Geomview_stream::setup_geomview(const char *machine, const char *login)
{
    int pipe_out[2], pipe_in[2];

    // Communication between CGAL and geomview should be possible
    // in two directions. To achieve this we open two pipes

    std::cout << "Starting Geomview..." << std::flush;
    if (pipe(pipe_out) < 0) {
        CGAL_error_msg( "out pipe failed" );
    }

    if (pipe(pipe_in) < 0) {
        CGAL_error_msg(  "in pipe failed" );
    }

    switch (pid = fork()){
    case -1:
        CGAL_error_msg( "fork failed" );
    case 0:               // The child process
        close(pipe_out[1]); // does not write to the out pipe,
        close(pipe_in[0]);  // does not read from the in pipe.

	if (dup2(pipe_out[0], 0) != 0)
	    std::cerr << "Connect pipe to stdin failed." << std::endl;
	if (dup2(pipe_in[1], 1) != 1)
	    std::cerr << "Connect pipe to stdout failed." << std::endl;

        if (machine && (std::strlen(machine)>0)) {
	    std::string s (" rgeomview ");
	    s += machine;
	    s += ":0.0";
            execlp("rsh", "rsh", machine, "-l", login, s.data(),
                   static_cast<void *>(NULL)); // cast to stop gcc warning
        } else {
            execlp("geomview", "geomview", "-c", "-",
                   static_cast<void *>(NULL)); // cast to stop gcc warning
        }

        // if we get to this point something went wrong.
        std::cerr << "execl geomview failed" << std::endl;
        switch(errno) {
        case EACCES:
            std::cerr << "please check your environment variable PATH"
		      << std::endl;
            std::cerr << "make sure the file `geomview' is contained in it"
		      << std::endl;
            std::cerr << "and is executable" << std::endl;
            break;
        case ELOOP:
            std::cerr << "too many links for filename `geomview'" << std::endl;
            break;
        default:
            std::cerr << "error number " << errno << " (check `man execlp')"
		      << std::endl;
        };
        CGAL_error();
    default:              // The parent process
        close(pipe_out[0]); // does not read from the out pipe,
        close(pipe_in[1]);  // does not write to the in pipe.

        in = pipe_in[0];
        out = pipe_out[1];

	// Necessary to wait a little bit for Geomview,
        // otherwise you won't be able to ask for points...
        sleep(1);

#if 1
        // We want to get rid of the requirement in the CGAL doc about
	// (echo "started").  But we want to be backward compatible, that is,
	// people who have this echo in their .geomview must still have CGAL
	// working, at least for a few public releases.
        // So the plan is to send, from CGAL, the command : (echo "CGAL-3D")
        // It's the same length as "started", 7.
        // Then we read 7 chars from Geomview, and test which string it is.
        // If it's "CGAL-3D", then fine, the user doesn't have .geomview with
        // the back-compatible echo command.
        // In the very long run, we'll be able to get rid of all this code as
        // well.
	// Maybe we should simply read the pipe, till we find "CGAL-3D" ?

        *this << "(echo \"CGAL-3D\")";

        char inbuf[10];
        std::size_t retread=::read(in, inbuf, 7);
        (void)retread;

        if (std::strncmp(inbuf, "started", 7) == 0)
        {
            // std::cerr << "You still have a .geomview file with the\n"
                   // << "(echo \"started\") command. Note that this is not\n"
                   // << "compulsory anymore, since CGAL 2.3" << std::endl;

            // Then the next one is supposed to be CGAL-3D.
            retread=::read(in, inbuf, 7);
            (void)retread;
            if (std::strncmp(inbuf, "CGAL-3D", 7) != 0)
                std::cerr << "Unexpected string from Geomview !" << std::endl;
        }
        else if (std::strncmp(inbuf, "CGAL-3D", 7) == 0)
        {
            // std::cerr << "Good, you don't have a .geomview file with the\n"
                      // << "(echo \"started\") command" << std::endl;
        }
        else
        {
            std::cerr << "Unexcepted string from Geomview at initialization!\n"
                      << "Going on nevertheless !" << std::endl;
        }
#else
        // Old original version
        char inbuf[10];
        // Waits for "started" from the .geomview file.
        retread=::read(in, inbuf, 7);
        (void)retread;
#endif

        std::cout << "done." << std::endl;

        (*this) << "(normalization g* none)(bbox-draw g* no)";
    }
}

void
Geomview_stream::pickplane(const Bbox_3 &bbox)
{
    bool bin_bak = set_binary_mode();
    (*this) << "(geometry pickplane {QUAD BINARY\n"
            << 1
    // here are the four corners
            << bbox.xmin() << bbox.ymin() << bbox.zmin()
            << bbox.xmin() << bbox.ymax() << bbox.zmin()
            << bbox.xmax() << bbox.ymax() << bbox.zmin()
            << bbox.xmax() << bbox.ymin() << bbox.zmin()

    // close the text bracket
            << "}) (pickable pickplane no)";
    set_ascii_mode(bin_bak);
}

void
Geomview_stream::clear()
{
    (*this) << "(delete World)";
    id.clear();
}

void
Geomview_stream::look_recenter()
{
    (*this) << "(look-recenter World)";
}

Geomview_stream&
Geomview_stream::operator<<(const std::string & s)
{
    if ((int)s.length() != ::write(out, s.data(), s.length())) {
        CGAL_error_msg( "write problem in the pipe while sending data to geomview" );
    }
    trace(s);

    return *this;
}

Geomview_stream&
Geomview_stream::operator<<(int i)
{
    // Depending on the mode chosen
    if (get_binary_mode()) {
        // we write raw binary data to the stream.
        int num = i;
        I_swap_to_big_endian(num);
        std::size_t retwrite=::write(out, (char*)&num, sizeof(num));
        (void)retwrite;
        trace(i);
    } else {
        // transform the int in a character sequence and put whitespace around
        std::ostringstream str;
        str << i << ' ' << std::ends;
        *this << str.str().c_str();
    }

    return *this;
}

Geomview_stream&
Geomview_stream::operator<<(unsigned int i)
{
    // Depending on the mode chosen
    if (get_binary_mode()) {
        // we write raw binary data to the stream.
        unsigned int num = i;
        I_swap_to_big_endian(num);
        std::size_t retwrite=::write(out, (char*)&num, sizeof(num));
        (void)retwrite;
        trace(i);
    } else {
        // transform the int in a character sequence and put whitespace around
        std::ostringstream str;
        str << i << ' ' << std::ends;
        *this << str.str().c_str();
    }

    return *this;
}

Geomview_stream&
Geomview_stream::operator<<(long i)
{
    return operator<<((int) i);
}

Geomview_stream&
Geomview_stream::operator<<(unsigned long i)
{
    return operator<<((unsigned int) i);
}

Geomview_stream&
Geomview_stream::operator<<(double d)
{
    float f = d;
    if (get_binary_mode()) {
        float num = d;
        I_swap_to_big_endian(num);
        std::size_t retwrite= ::write(out, (char*)&num, sizeof(num));
        (void)retwrite;
        trace(f);
    } else {
        // 'copy' the float in a string and append a blank
        std::ostringstream str;
        str << f << ' ' << std::ends;
        *this << str.str().c_str();
    }
    return *this;
}

Geomview_stream&
operator<<(Geomview_stream &gv, const Bbox_2 &bbox)
{
    bool ascii_bak = gv.set_ascii_mode();
    gv << "(geometry " << gv.get_new_id("Bbox")
       << " {VECT 1 5 0 5 0 ";
    // here are the four corners

    gv << bbox.xmin() << bbox.ymin() << 0.0
       << bbox.xmin() << bbox.ymax() << 0.0
       << bbox.xmax() << bbox.ymax() << 0.0
       << bbox.xmax() << bbox.ymin() << 0.0
       << bbox.xmin() << bbox.ymin() << 0.0;

    // close the text bracket
    gv << "})";
    gv.set_ascii_mode(ascii_bak);

    return gv;
}

Geomview_stream&
operator<<(Geomview_stream &gv, const Bbox_3 &bbox)
{
    bool ascii_bak = gv.set_ascii_mode();
    gv << "(geometry " << gv.get_new_id("Bbox")
       << " {appearance {material {edgecolor "
       << gv.ecr() << gv.ecg() << gv.ecb() <<  "}}{SKEL 8 4 "
    // here are the corners
       << bbox.xmin() << bbox.ymin() << bbox.zmin()
       << bbox.xmin() << bbox.ymax() << bbox.zmin()
       << bbox.xmax() << bbox.ymax() << bbox.zmin()
       << bbox.xmax() << bbox.ymin() << bbox.zmin()
       << bbox.xmax() << bbox.ymin() << bbox.zmax()
       << bbox.xmax() << bbox.ymax() << bbox.zmax()
       << bbox.xmin() << bbox.ymax() << bbox.zmax()
       << bbox.xmin() << bbox.ymin() << bbox.zmax()

       << "10 0 1 2 3 4 5 6 7 0 3\n"
       << "2 1 6\n"
       << "2 2 5\n"
       << "2 4 7\n"

    // close the text bracket
       << "}})";
    gv.set_ascii_mode(ascii_bak);

    return gv;
}

void
Geomview_stream::set_bg_color(const Color &c)
{
    bool ascii_bak = set_ascii_mode();
    *this << "(backcolor \"Camera\" "
          << double(c.r())/255.0
          << double(c.g())/255.0
          << double(c.b())/255.0
          << ")";
    set_ascii_mode(ascii_bak);
}

Geomview_stream&
Geomview_stream::operator<<(const Color &c)
{
    vertex_color = edge_color = face_color = c;
    return (*this);
}

Color
Geomview_stream::get_vertex_color() const
{
    return vertex_color;
}

Color
Geomview_stream::get_edge_color() const
{
    return edge_color;
}

Color
Geomview_stream::get_face_color() const
{
    return face_color;
}

Color
Geomview_stream::set_vertex_color(const Color &c)
{
    Color old = vertex_color;
    vertex_color = c;
    return old;
}

Color
Geomview_stream::set_edge_color(const Color &c)
{
    Color old = edge_color;
    edge_color = c;
    return old;
}

Color
Geomview_stream::set_face_color(const Color &c)
{
    Color old = face_color;
    face_color = c;
    return old;
}

double
Geomview_stream::vcr() const
{
    return double(vertex_color.r())/255.0;
}

double
Geomview_stream::vcg() const
{
    return double(vertex_color.g())/255.0;
}

double
Geomview_stream::vcb() const
{
    return double(vertex_color.b())/255.0;
}

double
Geomview_stream::ecr() const
{
    return double(edge_color.r())/255.0;
}

double
Geomview_stream::ecg() const
{
    return double(edge_color.g())/255.0;
}

double
Geomview_stream::ecb() const
{
    return double(edge_color.b())/255.0;
}

double
Geomview_stream::fcr() const
{
    return double(face_color.r())/255.0;
}

double
Geomview_stream::fcg() const
{
    return double(face_color.g())/255.0;
}

double
Geomview_stream::fcb() const
{
    return double(face_color.b())/255.0;
}

void
Geomview_stream::frame(const Bbox_3 &bbox)
{
    (*this) << bbox << "(look-recenter g0 c0)";
}

Geomview_stream&
Geomview_stream::operator>>(char *expr)
{
    // Skip whitespaces
    do {
      std::size_t retread=::read(in, expr, 1);
      (void)retread;
    } while (expr[0] != '(');

    int pcount = 1;
    int i = 1;
    while (1) {
        std::size_t retread=::read(in, &expr[i], 1);
        (void)retread;
        if (expr[i] == ')'){
            pcount--;
        } else if (expr[i] == '('){
            pcount++;
        }
        if (pcount == 0){
            expr[i+1]='\0';
            break;  // we encountered a balanced number of parantheses
        }
        i++;
    }
    return *this;
}

// Parse a Lisp expression, return a pointer to the beginning of the
// nth subexpression, and terminate it by '\0'.
// It's either a word terminated by ' ' or ')', or a well parenthesed
// expression, or a quoted "string".
char*
Geomview_stream::nth(char* s, int count)
{
    s++; // skip first character (always a parenthesis)

    // Skip "count" words.
    for(; count != 0; count--) {
        while (*s == ' ')       // skip whitespaces
            s++;
        s++;
        while (*s != ' ')       // skip a word
            s++;
    }
    while (*s == ' ')           // skip whitespaces
        s++;

    // Now we have the beginning of the searched sub-expression.
    int j = 1;
    if (*s == '(')              // Case of a well-parenthesed expression.
        for (int pcount = 1; pcount != 0;) {
            if (s[j] == ')') pcount--;
            if (s[j] == '(') pcount++;
            j++;
        }
    else if (*s == '"') {       // Case of a quoted "string".
        while (s[j] != '"')
            j++;
        j++;
    }
    else                        // Case of a word terminated by ' ' or ')'.
        while (s[j] != ' ' && s[j] != ')')
            j++;

    s[j] = '\0';
    return s;
}

void
Geomview_stream::parse_point(const char* pickpoint,
		     double &x, double &y, double &z, double &w)
{
    std::stringstream ss;
    ss << pickpoint << std::ends;

    char parenthesis;
    ss >> parenthesis >> x >> y >> z >> w;
}

std::string
Geomview_stream::get_new_id(const std::string & s)
{
    std::ostringstream str;
    str << s << id[s]++ << std::ends;
    return str.str();
}

} //namespace CGAL

#else

// Add a dummy symbol to prevent warnings of empty translation unit.
namespace CGAL {
namespace {
int dummy;
}
} //namespace CGAL

#endif // CGAL_USE_GEOMVIEW
