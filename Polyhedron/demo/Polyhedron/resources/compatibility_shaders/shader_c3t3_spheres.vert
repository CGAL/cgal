//#version 100 
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
uniform highp float point_size;
void main(void)
{
  gl_PointSize = point_size;
  color = vec4(colors, center.x * cutplane.x  + center.y * cutplane.y  + center.z * cutplane.z  +  cutplane.w);
  highp vec4 my_vertex = vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, radius*vertex.z + center.z, 1.0) ;
  fP = mv_matrix * my_vertex;
  mat3 mv_matrix_3;                    
  mv_matrix_3[0] = mv_matrix[0].xyz;   
  mv_matrix_3[1] = mv_matrix[1].xyz;   
  mv_matrix_3[2] = mv_matrix[2].xyz;   
  fN = mv_matrix_3* normals;           
  gl_Position =  mvp_matrix * my_vertex;
}
