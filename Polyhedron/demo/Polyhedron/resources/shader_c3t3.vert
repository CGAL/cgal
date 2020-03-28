#version 150
in vec4 vertex;
in vec3 normals;
in vec3 colors;
in vec3 center;
uniform  mat4 mvp_matrix;
uniform  mat4 mv_matrix;
uniform  vec4 cutplane;
uniform  float shrink_factor;
out vec4 fP; 
out vec3 fN; 
out vec4 color; 
uniform  float point_size;
void main(void)
{
  gl_PointSize = point_size;
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  fP = mv_matrix * vertex; 

  mat3 mv_matrix_3;                    
  mv_matrix_3[0] = mv_matrix[0].xyz;   
  mv_matrix_3[1] = mv_matrix[1].xyz;   
  mv_matrix_3[2] = mv_matrix[2].xyz;   
  fN = mv_matrix_3* normals;           
  
   mat4 transOB = mat4(1, 0, 0, 0, // first column
   0, 1, 0, 0, // second column
   0, 0, 1, 0, // third column
   center.x, center.y, center.z, 1); // fourth column
   mat4 transBO = mat4(1, 0, 0, 0, // first column
    0, 1, 0, 0, // second column
    0, 0, 1, 0, // third column
    -center.x, -center.y, -center.z, 1); // fourth column
    mat4 scaling = mat4(shrink_factor, 0, 0, 0,
    0, shrink_factor, 0, 0,
    0, 0, shrink_factor, 0,
    0, 0, 0, 1);
  gl_Position = mvp_matrix *transOB * scaling * transBO * vertex;
}
