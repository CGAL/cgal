#version 150
in vec4 vertex;
in vec3 normals;
in vec3 colors;
in vec3 center;
in float radius;
uniform vec4 cutplane;
uniform mat4 mvp_matrix;
uniform mat4 mv_matrix;
out vec4 fP;
out vec3 fN;
out vec4 color;
uniform float point_size;
void main(void)
{
  gl_PointSize = point_size;
  color = vec4(colors, center.x * cutplane.x  + center.y * cutplane.y  + center.z * cutplane.z  +  cutplane.w);
  vec4 my_vertex = vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, radius*vertex.z + center.z, 1.0) ;
  fP = mv_matrix * my_vertex;
  mat3 mv_matrix_3;                    
  mv_matrix_3[0] = mv_matrix[0].xyz;   
  mv_matrix_3[1] = mv_matrix[1].xyz;   
  mv_matrix_3[2] = mv_matrix[2].xyz;   
  fN = mv_matrix_3* normals;           
  gl_Position =  mvp_matrix * my_vertex;
}
