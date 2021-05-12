
#include <CGAL/basic.h>
#include<CGAL/Handle_for.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/MP_Float.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>

#include <CGAL/intersections.h>

#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_arc_traits.h>
#include <CGAL/Circular_arc_traits_tracer.h>

#include <CGAL/Lazy_circular_kernel_2.h>

#include <CGAL/Filtered_hexagon_circular_kernel_2.h>

#include <CGAL/Filtered_bbox_circular_kernel_2.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Variant_traits.h>


#include "IntStr.h"
#include "benchmark_parser.h"
#include <iostream>
#include <sstream>
template <typename K,class ArcContainer>
struct Breader : public Benchmark_visitor {
        typedef typename CGAL::Quotient<CGAL::MP_Float>                       NT1;
        typedef typename K::Point_2 Point_2;
          typedef typename K::Circle_2 Circle_2;
          typedef std::list<Circle_2> Circles;
        typedef std::list<Point_2> Points;
        typedef std::list<NT1> Numbers;
          typedef std::list<std::pair<Point_2, double> > Polygon;
          typedef std::list<Polygon> Polygons;
        Points points;
        Numbers numbers;
        Circles circles;
        Polygons polygons;
        int x_;
        int y_;
    Breader() {}
    virtual void token_not_handled( std::string s) {}
    virtual void accept_benchmark_name( std::string s) {
        Benchmark_visitor::accept_benchmark_name(s);
        std::cerr << "name '" << s << "', ";
    }

    virtual void accept_classification( std::string problem,
                                        std::string geom,
                                        std::string clas,
                                        std::string family,
                                        std::string instance,
                                        std::string release) {

        if ((problem != "Arrangement") && (problem != "CSG")
            && (problem != " "))
            error_handler
                ( "classification error");

        if ((geom != "Lines") && (geom != "Circles") && (geom != "Conics")
            && (geom != "Cubics") && (geom != "Quartics")
            && (geom != "ArbitraryDegreeCurves") && ( geom != "Quadrics")
            && (geom != "Tori") && (geom != "Planes") && (geom != " "))
            error_handler ( "classification error" );

        if ((clas != "FullCurves") && (clas != "Ellipses")
            && (clas != "BoundedArcs") && (clas != "UnboundedArcs")
            && (clas != "FullSurfaces") && (clas != "BoundedPatches")
            && (clas != "UnboundedPatches") && (clas != " "))
            error_handler ( "classification error" );
    }
        virtual void begin_circle_2() {

                std::cerr << "circle begin" << std::endl;
                x_=0;
                y_=0;
        }

        virtual void end_circle_2() {

                std::cerr << "circle end" << std::endl;
                //Circle_2 cyrc = typename K::Construct_circle_2()(Point_2(x_,y_), integer*integer);
                NT1 r= numbers.back();
                Circle_2 circ = typename K::Construct_circle_2()(points.back(), r*r);
                circles.push_back(circ);
        }

        virtual void accept_point_2( std::string x, std::string y) {

                std::cerr << "accept_point_2 :"<<x << " ; "<< y << std::endl;

                NT1 qx,qy;
                from_string(qx,x);
                from_string(qy,y);
                points.push_back(typename K::Construct_point_2()(qx,qy));
                std::cerr << "accept_point_2 :"<<qx << " ; "<< qy << std::endl;
        }
        virtual void accept_rational( std::string num, std::string denom) {


                CGAL::MP_Float _num,_denom;
                from_string(_num,num);
                from_string(_denom,denom);
                NT1 x(_num,_denom);
                std::cerr << "accept rat:"<<x<< std::endl;
                numbers.push_back(x);
        }

        virtual void accept_point_2( std::string x_num, std::string x_denom,
                                 std::string y_num, std::string y_denom) {


                CGAL::MP_Float num,denom;
                from_string(num,x_num);
                from_string(denom,x_denom);
                NT1 x(num,denom);
                from_string(num,y_num);
                from_string(denom,y_denom);
                NT1 y(num,denom);
                std::cerr << "accept_point_2(rational) :"<< x<<","<< y<< std::endl;
                tnh( "Point_2(x_num, x_denom, y_num, y_denom)");
                points.push_back(typename K::Construct_point_2()(x,y));
            }
        virtual void accept_integer( std::string s) {

                NT1 n;
                from_string(n,s);
                std::cerr <<"accept_integer!!!!"<<n <<std::endl;
                numbers.push_back(n);

        }
        virtual void begin_line_segment_2() {
                std::cerr <<"Begin_line_segment_2"<<std::endl; }

        virtual void end_line_segment_2() {
                std::cerr <<"end_line_segment_2"<<std::endl;
                Polygon p1;
                Polygon p2;
                typename Points::iterator pit = points.end();
                p2.push_back(std::make_pair(*pit, 0));
                pit--;
                p1.push_back(std::make_pair(*pit, 0));
                polygons.push_back(p1);
                polygons.push_back(p2);
        }
            virtual void begin_CircularArc_2(){
                std::cerr <<"begin_Circular_arc_2"<<std::endl;
        }

            virtual void end_CircularArc_2(){
                std::cerr <<"end_Circular_arc_2"<<std::endl;
        }
   virtual void begin_LineArc_2(){
        tnh("begin_line_arc_2");
        std::cerr <<"begin_line_arc_2"<<std::endl;
        }

    virtual void end_LineArc_2(){
        tnh("end_line_arc_2");
        std::cerr <<"end_line_arc_2"<<std::endl;
        }
   virtual void begin_CircularPoint_2(){
        tnh("begin_circular_arc_poin");
        std::cerr <<"begin_circular_arc_poin"<<std::endl;
        }

    virtual void end_CircularPoint_2(){
        tnh("end_circular_arc_poin");
        std::cerr <<"end_circular_arc_poin"<<std::endl;
        }

};
////
int str_to_int(const char* ch){
int i;
sscanf(ch,"%i",&i);
return i;
}

int main( int argc, char* argv[] ) {

  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_circles_2_2<NT1>                      Algebraic_k1;
  typedef CGAL::Circular_kernel_2<Linear_k1, Algebraic_k1>         CK;
  typedef CK::Circular_arc_2                                  Circular_arc_2;
  typedef CK::Line_arc_2                                      Line_arc_2;
  typedef CGAL::Variant_traits<CK,Line_arc_2,Circular_arc_2>  CircularK_Variant_Traits;
  typedef boost::variant< Circular_arc_2, Line_arc_2 >        CircularKVarArc;
  typedef std::vector<CircularKVarArc>                        CircularKVarArcContainer;
        CircularKVarArcContainer arc;
        typedef CGAL::Arrangement_2<CircularK_Variant_Traits>                 Pmwx;
          typedef CGAL::Arr_naive_point_location<Pmwx>        Point_location;
        Pmwx _pm;
          Point_location _pl(_pm);
        typedef CK::Circle_2 Circle_2;
        typedef std::list<Circle_2> Circles;

int exit_status = 0;
    if ( argc < 2) {
        return 0;
    } else {
        for ( int i = 1; i < argc; ++i) {
           Breader<CK,CircularKVarArcContainer> check;
            std::cerr << "File '" << argv[i] << "', "<<std::endl;
            if ( benchmark_parse_file( argv[i], &check)) {
                std::cerr<<check.circles.size()<<std::endl;
                std::cerr<<check.polygons.size()<<std::endl;
            } else {
                std::cerr << "is malformed." << std::endl;
                exit_status = 1;

            }
        }
    }


    return exit_status;

}

