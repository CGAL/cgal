#version 120
attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec3 colors;
attribute highp vec3 center;
attribute highp float radius;
uniform highp vec4 cutplane;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 color;


void main(void)
{
  color = vec4(colors, center.x * cutplane.x  + center.y * cutplane.y  + center.z * cutplane.z  +  cutplane.w);
  fP = mv_matrix * vertex;
  fN = mat3(mv_matrix)* normals;
  gl_Position =  mvp_matrix *
  vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, radius*vertex.z + center.z, 1.0) ;
}
