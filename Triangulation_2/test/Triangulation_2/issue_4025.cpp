#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
typedef CGAL::Polygon_2<K>                                                Polygon_2;
typedef CGAL::Exact_intersections_tag                                     Itag_;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default, Itag_> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>                       CDTP;

typedef CDTP::Point                                                       Point;
typedef CDTP::Constraint_id                                               Cid;
typedef CDTP::Vertex_handle                                               Vertex_handle;
typedef CDTP::Constraint_id                                               Constraint_id;
typedef CDTP::Vertices_in_constraint_iterator                             Vertices_in_constraint_iterator;

int countVertex(CDTP &cdtp, CDTP::Constraint_id id)
{
    Vertices_in_constraint_iterator v=cdtp.vertices_in_constraint_begin(id);

    int count=0;
    while(v!=cdtp.vertices_in_constraint_end(id))
    {
        count++;
        v++;
    }

    return count;
}


int main()
{
    CDTP cdtp;

    std::list<Point> pointsListCollinear;

    pointsListCollinear.push_back(Point(0,0));
    pointsListCollinear.push_back(Point(0,1));
    pointsListCollinear.push_back(Point(0,2));
    pointsListCollinear.push_back(Point(0,3));
    pointsListCollinear.push_back(Point(0,4));
    pointsListCollinear.push_back(Point(0,5));

    std::list<Point> pointsListNoCollinear;

    pointsListNoCollinear.push_back(Point(1,0));
    pointsListNoCollinear.push_back(Point(2,1));
    pointsListNoCollinear.push_back(Point(4,2));
    pointsListNoCollinear.push_back(Point(2,3));
    pointsListNoCollinear.push_back(Point(4,4));
    pointsListNoCollinear.push_back(Point(1,5));


    Constraint_id ctIdCollinear=cdtp.insert_constraint(pointsListCollinear.begin(),pointsListCollinear.end());
    Constraint_id ctIdNoCollinear=cdtp.insert_constraint(pointsListNoCollinear.begin(),pointsListNoCollinear.end());


    //******************************* attempt with the collinear constraint
    Vertices_in_constraint_iterator vertexToRemoveCollinear=cdtp.vertices_in_constraint_begin(ctIdCollinear);
    vertexToRemoveCollinear++;
    vertexToRemoveCollinear++;


        std::cout<<"attempt to remove vertex "<<(*vertexToRemoveCollinear)->point().x()<<" , "<<(*vertexToRemoveCollinear)->point().y() <<std::endl;
    cdtp.remove_vertex_from_constraint(ctIdCollinear,vertexToRemoveCollinear);

    std::cout<<"number of subconstraints "<<cdtp.number_of_subconstraints()<<std::endl; //--> 5, expected 4
    std::cout<<"number of constraints "<<cdtp.number_of_constraints()<<std::endl; //--> 1
    std::cout<<"number of vertex in constraint "<<countVertex(cdtp,ctIdCollinear)<<std::endl; //--> 6, expected 5


    //******************************* attempt with the collinear constraint
    Vertices_in_constraint_iterator vertexToRemoveNoCollinear=cdtp.vertices_in_constraint_begin(ctIdNoCollinear);
    vertexToRemoveNoCollinear++;
    vertexToRemoveNoCollinear++;

        std::cout<<"attempt to remove vertex "<<(*vertexToRemoveNoCollinear)->point().x()<<" , "<<(*vertexToRemoveNoCollinear)->point().y() << std::endl;
    cdtp.remove_vertex_from_constraint(ctIdNoCollinear,vertexToRemoveNoCollinear);

    std::cout<<"number of subconstraints "<<cdtp.number_of_subconstraints()<<std::endl; //--> 4, ok
    std::cout<<"number of constraints "<<cdtp.number_of_constraints()<<std::endl;  //--> 1
    std::cout<<"number of vertex in constraint "<<countVertex(cdtp,ctIdNoCollinear)<<std::endl;  //--> 5, ok

    return 0;

}
