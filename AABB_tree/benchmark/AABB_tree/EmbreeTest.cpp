#include <iostream>
#include <fstream>
#include <limits>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <embree3/rtcore.h> 

#include <CGAL/Timer.h>

#include "RaysGenerate.h"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

struct Vertex   { float x,y,z,r; }; 
struct Triangle { int v0, v1, v2; };

int main(int argc, char *argv[])
{   
    const char* filename = (argc > 1) ? argv[1] : "data/data.ply";
    std::ifstream input(filename);

    Mesh surfaceMesh;
    CGAL::read_ply(input, surfaceMesh);

    RTCDevice device = rtcNewDevice("verbose=0");
    RTCScene scene = rtcNewScene(device);

    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),surfaceMesh.number_of_vertices());

    for(vertex_descriptor vd : surfaceMesh.vertices()){
        Point data = surfaceMesh.point(vd);
        vertices[vd.idx()].x = data.x();       
        vertices[vd.idx()].y = data.y();       
        vertices[vd.idx()].z = data.z();       
    }

    Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),surfaceMesh.number_of_faces());
    
    for (face_descriptor fd : surfaceMesh.faces()){
        Mesh::Halfedge_index hf = surfaceMesh.halfedge(fd);
        int temp[3]; int i=0;
        for(Mesh::Halfedge_index hi : halfedges_around_face(hf, surfaceMesh)){
            Mesh::Vertex_index vi = target(hi, surfaceMesh);
            temp[i] = vi.idx();
            i++;
        }
        triangles[fd.idx()].v0 = temp[0];
        triangles[fd.idx()].v1 = temp[1];
        triangles[fd.idx()].v2 = temp[2];
    }

    rtcCommitGeometry(mesh);
    unsigned int geomID = rtcAttachGeometry(scene, mesh);
    rtcReleaseGeometry(mesh);

    rtcCommitScene(scene);
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRayHit rayhit;
    rayhit.ray.org_x = -2.0; /*POINT.X*/ 
    rayhit.ray.org_y =  0.0; /*POINT.Y*/
    rayhit.ray.org_z =  0.0; /*POINT.Z*/

    rayhit.ray.tnear = 0.0;
    rayhit.ray.tfar = std::numeric_limits<double>::infinity();
    rayhit.ray.flags = 0;

    // rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;
    // rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    
    int numberOfRays = 50000; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays); 
    CGAL::Timer time;
    time.start();
    for(size_t n=0; n!=numberOfRays; ++n){
        rayhit.ray.dir_x = rg.rayDirections[n]._x;
        rayhit.ray.dir_x = rg.rayDirections[n]._y;
        rayhit.ray.dir_x = rg.rayDirections[n]._z;

        rtcIntersect1(scene, &context, &rayhit);

    }
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;
        
    rtcReleaseDevice(device);
    return 0;
}
