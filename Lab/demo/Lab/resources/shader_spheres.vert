#version 150
in vec4 vertex;
in vec3 normals;
in vec3 colors;
in vec3 center;
in float radius;
uniform mat4 mvp_matrix;
uniform mat4 mv_matrix;
uniform mat4 norm_matrix;
uniform mat4 f_matrix;
out vec4 fP;
out vec3 fN;
out vec4 color;
out float dist[6];
uniform float point_size;

void main(void)
{
  gl_PointSize = point_size;
  for(int i=0; i<6; ++i)
    dist[i] = 1.0;
  color = vec4(colors, 1.0);
  vec4 transformed_vertex =
      vec4(radius*vertex.x + center.x, radius* vertex.y + center.y, radius*vertex.z + center.z, 1.0);
  fP = mv_matrix * transformed_vertex;
  mat3 norm_matrix_3;
  norm_matrix_3[0] = norm_matrix[0].xyz;
  norm_matrix_3[1] = norm_matrix[1].xyz;
  norm_matrix_3[2] = norm_matrix[2].xyz;
  fN = norm_matrix_3* normals;
  gl_Position =  mvp_matrix * f_matrix * transformed_vertex;
}
