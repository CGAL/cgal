#include <iostream>
#include <limits>

#include <embree3/rtcore.h> 

struct Vertex   { float x,y,z,r;  }; 
struct Triangle { int v0, v1, v2; };

int main(int argc, char const *argv[])
{   
      
    RTCDevice device = rtcNewDevice("verbose=0");
    RTCScene scene = rtcNewScene(device);

    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),8);
    vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1;
    vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1;
    vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1;
    vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1;
    vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1;
    vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1;
    vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1;
    vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1;

    int tri = 0;
    Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),12);

    triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 2; tri++;
    triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 2; tri++;

    triangles[tri].v0 = 4; triangles[tri].v1 = 6; triangles[tri].v2 = 5; tri++;
    triangles[tri].v0 = 5; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;

    triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 1; tri++;
    triangles[tri].v0 = 1; triangles[tri].v1 = 4; triangles[tri].v2 = 5; tri++;

    triangles[tri].v0 = 2; triangles[tri].v1 = 3; triangles[tri].v2 = 6; tri++;
    triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;

    triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 4; tri++;
    triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 4; tri++;

    triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 3; tri++;
    triangles[tri].v0 = 3; triangles[tri].v1 = 5; triangles[tri].v2 = 7; tri++;

    rtcCommitGeometry(mesh);
    unsigned int geomID = rtcAttachGeometry(scene, mesh);
    rtcReleaseGeometry(mesh);

    rtcCommitScene(scene);
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRayHit rayhit;
    rayhit.ray.org_x = -2.0;
    rayhit.ray.org_y = 0;
    rayhit.ray.org_z = 0;
    rayhit.ray.tnear = 0.0;

    rayhit.ray.dir_x = 1;
    rayhit.ray.dir_x = 0;
    rayhit.ray.dir_x = 0;

    rayhit.ray.tfar = std::numeric_limits<double>::infinity();
    rayhit.ray.flags = 0;

    rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(scene, &context, &rayhit);

    if(rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
            std::cout<<"Intersection."<<std::endl;

    rtcReleaseDevice(device);
    return 0;
}
