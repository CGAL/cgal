// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//

// release       :
// release_date  :
//
// file          : src/Geomview_stream.C
// source        : web/Geomview_stream.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri and Herve Bronnimann
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// ============================================================================


#include <CGAL/IO/Geomview_stream.h>

#include <unistd.h>
#include <errno.h>
#include <stdlib.h>


CGAL_Geomview_stream::CGAL_Geomview_stream(const CGAL_Bbox_3 &bbox,
                                           const char *machine,
                                           const char *login)
    : _line_width(1)
{
    setup_geomview(machine,login);
    frame(bbox);
    pickplane(bbox);
    set_vertex_radius((bbox.xmax() - bbox.xmin())/100.0);
}

CGAL_Geomview_stream::CGAL_Geomview_stream(const char *machine,
                                           const char *login,
                                           const CGAL_Bbox_3 &bbox)
    : _line_width(1)
{
    cerr << "Warning: This constructor is going to disappear" << endl
         << "The bounding box should come as first argument" << endl
         << "machine and login default to NULL" << endl;
    setup_geomview(machine,login);
    frame(bbox);
    pickplane(bbox);
    set_vertex_radius((bbox.xmax() - bbox.xmin())/100.0);
}

CGAL_Geomview_stream::~CGAL_Geomview_stream()
{
    kill(pid, SIGKILL);  // kills geomview
}


void CGAL_Geomview_stream::setup_geomview(const char *machine,
                                          const char *login)
{
    bflag = 0;
    _trace = false;
    col = CGAL_BLACK;
    vertex_color = CGAL_BLACK;
    edge_color = CGAL_BLACK;
    face_color = CGAL_BLACK;
    int pipe_out[2], pipe_in[2];

    // Communication between CGAL and geomview should be possible
    // in two directions. To achieve this we open two pipes

    cout << "Starting Geomview..." << flush ;
    if (pipe(pipe_out) < 0) {
        cerr << "out pipe failed" << endl ;
        exit(-1);
    }

    if (pipe(pipe_in) < 0) {
        cerr << "in pipe failed" << endl ;
        exit(-1);
    }

    switch (pid = fork()){
    case -1:
        cerr << "fork failed" << endl ;
        exit(-1);
    case 0:               // The child process
        close(pipe_out[1]); // does not write to the out pipe,
        close(pipe_in[0]);  // does not read from the in pipe.


        close (0);          // this is the file descriptor of cin
        dup(pipe_out[0]);   // we connect it to the pipe
        close (1);          // this is the file descriptor of cout
        dup(pipe_in[1]);    // we connect it to the pipe
        if (machine && (strlen(machine)>0)) {
            ostrstream os;
            os << " rgeomview " << machine << ":0.0" << ends ;
            ostrstream logos;
            execlp("rsh", "rsh", machine, "-l", login, os.str(), (char *)0);
        } else {
            execlp("geomview", "geomview", "-c", "-", (char *)0);
        }

        // if we get to this point something went wrong.
        cerr << "execl geomview failed" << endl ;
        switch(errno) {
        case EACCES:
            cerr << "please check your environment variable PATH" << endl;
            cerr << "make sure the file `geomview' is contained in it" << endl;
            cerr << "and is executable" << endl;
            break;
        case ELOOP:
            cerr << "too many links for filename `geomview'" << endl;
            break;
        default:
            cerr << "error number " << errno << " (check `man execlp')" << endl;
        };
        exit(-1);
    default:              // The parent process
        close(pipe_out[0]); // does not read from the out pipe,
        close(pipe_in[1]);  // does not write to the in pipe.

        in = pipe_in[0];
        out = pipe_out[1];

        char inbuf[10];
        read(in, inbuf, 7);

        cout << "done." << endl;

        bbox_count = 0;
        triangle_count = 0;
        segment_count = 0;
        point_count = 0;
        tetrahedron_count = 0;
        (*this) << "(normalization g* none)(bbox-draw g* no)" ;

        break;
    }
}


void
CGAL_Geomview_stream::pickplane(const CGAL_Bbox_3 &bbox)
{
    (*this) << binary
            << "(geometry pickplane {QUAD BINARY\n"
            << 1
    // here are the four corners
            << bbox.xmin() << bbox.ymin() << bbox.zmin()
            << bbox.xmin() << bbox.ymax() << bbox.zmin()
            << bbox.xmax() << bbox.ymax() << bbox.zmin()
            << bbox.xmax() << bbox.ymin() << bbox.zmin()

    // close the text bracket
            << "}) (pickable pickplane no)"
            << ascii ;
}


void
CGAL_Geomview_stream::set_binary_mode()
{
    bflag = 1;
}

void
CGAL_Geomview_stream::set_ascii_mode()
{
    bflag = 0;
}

bool
CGAL_Geomview_stream::in_binary_mode() const
{
    return bflag;
}

bool
CGAL_Geomview_stream::in_ascii_mode() const
{
    return ! bflag;
}
CGAL_Geomview_stream&
CGAL_Geomview_stream::operator<<
(CGAL_Geomview_stream&(*fct)(CGAL_Geomview_stream&))
{
  (*fct)(*this);
  return *this;
}
bool
CGAL_Geomview_stream::get_trace() const
{
    return _trace;
}

bool
CGAL_Geomview_stream::set_trace(bool b)
{
    bool old = _trace;
    _trace = b;
    return old;
}


void
CGAL_Geomview_stream::trace(const char *cptr) const
{
    if(_trace){
        cerr << cptr;
    }
}

void
CGAL_Geomview_stream::trace(double d) const
{
    if(_trace){
        cerr << d << ' ';
    }
}

void
CGAL_Geomview_stream::trace(int i) const
{
    if(_trace){
        cerr << i << ' ';
    }
}
double
CGAL_Geomview_stream::get_vertex_radius() const
{
    return _radius;
}


double
CGAL_Geomview_stream::set_vertex_radius(double r)
{
    double old = _radius;
    _radius = r;
    return old;
}
int
CGAL_Geomview_stream::get_line_width() const
{
    return _line_width;
}


int
CGAL_Geomview_stream::set_line_width(int w)
{
    int old = _line_width;
    _line_width = w;
    return old;
}
void
CGAL_Geomview_stream::clear()
{
    (*this) << "(delete World)";
}

void
CGAL_Geomview_stream::look_recenter() const
{
    CGAL_Geomview_stream* ncthis = (CGAL_Geomview_stream*)this;
    (*ncthis) << "(look-recenter World)";
}
CGAL_Geomview_stream&
CGAL_Geomview_stream::operator<<(const char *cptr)
{
    int length = strlen(cptr);
    if (length != write(out, cptr, length)) {
        cerr << "write problem in the pipe while sending data to geomview"
             << endl;
        exit(-1);
    }
    trace(cptr);

    return *this;
}
CGAL_Geomview_stream&
CGAL_Geomview_stream::operator<<(int i)
{
    // Depending on the mode chosen
    if (in_binary_mode()) {
        // we write raw binary data to the stream.
        write(out, (char*)&i, sizeof(i));
    } else {
        // transform the int in a character sequence and put whitespace around
        ostrstream str;
        str << i << ' ' << ends;
        char *bptr = str.str();
        write(out, bptr, int(strlen(bptr)));
    }
    trace(i);

    return *this;
}
CGAL_Geomview_stream&
CGAL_Geomview_stream::operator<<(double d)
{
    float f = d;

    if (in_binary_mode()) {
        write(out, (char*)&f, sizeof(f));
    } else {
        // 'copy' the float in a string and append a blank
        ostrstream str;
        str << f << " " << ends ;
        char *bptr = str.str();

        write(out, bptr, int(strlen(bptr)));
    }
    trace(f);
    return *this;
}


CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Bbox_2 &bbox)
{
    ostrstream os;
    os << "bbox" << gv.bbox_count++ << ends ;
    char *id = os.str();

    gv << ascii
       << "(geometry " << id << " {VECT 1 5 0 5 0 " ;
    // here are the four corners

    gv << bbox.xmin() << bbox.ymin() << 0.0
       << bbox.xmin() << bbox.ymax() << 0.0
       << bbox.xmax() << bbox.ymax() << 0.0
       << bbox.xmax() << bbox.ymin() << 0.0
       << bbox.xmin() << bbox.ymin() << 0.0 ;

    // close the text bracket
    gv << "})" ;

    return gv;
}
CGAL_Geomview_stream&
operator<<(CGAL_Geomview_stream &gv,
           const CGAL_Bbox_3 &bbox)
{
    ostrstream os;
    os << "bbox" << gv.bbox_count++ << ends ;
    char *id = os.str();


    gv << ascii
       << "(geometry " << id << " {appearance {material {edgecolor "
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
       << "}})" ;

    return gv;
}
void
CGAL_Geomview_stream::set_bg_color(const CGAL_Color &c)
{
    *this << ascii
          << "(backcolor \"Camera\" "
          << double(c.r())/255.0
          << double(c.g())/255.0
          << double(c.b())/255.0
          << ")";
}

CGAL_Geomview_stream&
CGAL_Geomview_stream::operator<<(const CGAL_Color &c)
{
    col = c;
    vertex_color = c;
    edge_color = c;
    face_color = c;
    return (*this);
}


CGAL_Color
CGAL_Geomview_stream::get_vertex_color() const
{
    return vertex_color;
}


CGAL_Color
CGAL_Geomview_stream::get_edge_color() const
{
    return edge_color;
}


CGAL_Color
CGAL_Geomview_stream::get_face_color() const
{
    return face_color;
}


CGAL_Color
CGAL_Geomview_stream::set_vertex_color(const CGAL_Color &c)
{
    CGAL_Color old = vertex_color;
    vertex_color = c;
    return old;
}


CGAL_Color
CGAL_Geomview_stream::set_edge_color(const CGAL_Color &c)
{
    CGAL_Color old = edge_color;
    edge_color = c;
    return old;
}


CGAL_Color
CGAL_Geomview_stream::set_face_color(const CGAL_Color &c)
{
    CGAL_Color old = face_color;
    face_color = c;
    return old;
}


double
CGAL_Geomview_stream::vcr() const
{
    return double(vertex_color.r())/255.0;
}

double
CGAL_Geomview_stream::vcg() const
{
    return double(vertex_color.g())/255.0;
}

double
CGAL_Geomview_stream::vcb() const
{
    return double(vertex_color.b())/255.0;
}


double
CGAL_Geomview_stream::ecr() const
{
    return double(edge_color.r())/255.0;
}

double
CGAL_Geomview_stream::ecg() const
{
    return double(edge_color.g())/255.0;
}

double
CGAL_Geomview_stream::ecb() const
{
    return double(edge_color.b())/255.0;
}

double
CGAL_Geomview_stream::fcr() const
{
    return double(face_color.r())/255.0;
}

double
CGAL_Geomview_stream::fcg() const
{
    return double(face_color.g())/255.0;
}

double
CGAL_Geomview_stream::fcb() const
{
    return double(face_color.b())/255.0;
}


void
CGAL_Geomview_stream::frame(const CGAL_Bbox_3 &bbox)
{
    (*this) << bbox
            << ascii
            << "(look-recenter g0 c0)(delete bbox0)" ;
}


CGAL_Geomview_stream&
CGAL_Geomview_stream::operator>>(char *expr)
{
    // skip whitespace
    read(in, expr, 1);
    while(expr[0] != '('){
        read(in, expr, 1);
    }
    int pcount = 1;
    int i = 1;
    while (1) {
        read(in, &expr[i], 1);
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
    return *this ;
}
char*
CGAL_nth(char* s, int count)
{
    int length = strlen(s);

    // the first character is always a parenthesis
    int i = 1,  // char count
        wc = 0; // word count

    while (1) {
        // skip whitespace
        while(s[i]==' '){
            i++;
        }
        int j = 0;
        int pcount = 1;
        switch (s[i]){
        case '(':
            j=i+1;
            while (1) {
                if ( s[j]==')' ) {
                    pcount--;
                } else if (s[j]=='(' ) {
                    pcount++;
                }
                if ((pcount == 0 )&& (wc == count)) {
                    s[j+1] = '\0';
                    return s + i;
                }
                j++;
            }
            i = j+1;
            break;
        case '"':
            j = i+1;
            while (1) {
                if ( (s[j]=='"') && (wc == count) ) {
                    s[j+1] = '\0';
                    return s + i;
                }
                j++;
            }
            i = j+1;
            break;
        default:
            j=i+1;
            while( ( s[j]!=' ' ) && ( s[j]!=')' ) ){
                j++;
            }
            if (wc == count){
                s[j] = '\0';
                return s + i;
            }
            i = j;
            break;
        }
        wc++;
    }

    return s;
}
bool
CGAL_is_prefix(const char* p, const char* w)
{
    while((*p != '\0') && (*w != '\0')){
        if(*p != *w){
            return false;
        }
        ++p;
        ++w;
    }
    if((*w == '\0') && (*p != '\0')){
        return false;
    }
    return true;
}


