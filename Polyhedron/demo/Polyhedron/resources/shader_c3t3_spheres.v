#version 120

attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec4 colors;
attribute highp vec3 barycenter;
uniform highp vec4 cutplane;
attribute highp float radius;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 color;
varying highp float dist[6];

void main(void)
{
  color = vec4(colors.xyz, barycenter.x * cutplane.x  + barycenter.y * cutplane.y  + barycenter.z * cutplane.z  +  cutplane.w);
  for(int i=0; i<6; ++i)
    dist[i] = 1;
  vec4 my_vertex =
  vec4(radius*vertex.x + barycenter.x, radius* vertex.y + barycenter.y, radius*vertex.z + barycenter.z, 1.0) ;
  fP = mv_matrix * my_vertex;
  fN = mat3(mv_matrix)* normals;
  gl_Position =  mvp_matrix * my_vertex ;
}
