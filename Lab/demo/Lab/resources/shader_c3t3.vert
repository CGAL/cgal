#version 150
in vec4 vertex;
in vec3 normals;
in vec3 colors;
in vec3 center;
in vec2 subdomain_in;
uniform  mat4 mvp_matrix;
uniform  mat4 mv_matrix;
uniform mat4 norm_matrix;
uniform  vec4 cutplane;
uniform  float shrink_factor;
out vec4 fP;
out vec3 fN;
out vec4 color;
out float dist[6];
flat out vec2 subdomain_out;
uniform bool is_clipbox_on;
uniform mat4 clipbox1;
uniform mat4 clipbox2;
uniform  float point_size;

void compute_distances(void)
{
  for(int i=0; i<3; ++i)
  {
    dist[i]=
    clipbox1[i][0]*vertex.x+
    clipbox1[i][1]*vertex.y+
    clipbox1[i][2]*vertex.z +
    clipbox1[i][3];
    dist[i+3]=
    clipbox2[i][0]*vertex.x+
    clipbox2[i][1]*vertex.y+
    clipbox2[i][2]*vertex.z +
    clipbox2[i][3];
  }
}

void main(void)
{
  subdomain_out = subdomain_in;
  gl_PointSize = point_size;
  color = vec4(colors, vertex.x * cutplane.x  + vertex.y * cutplane.y  + vertex.z * cutplane.z  +  cutplane.w);
  if(is_clipbox_on)
   compute_distances();
  fP = mv_matrix * vertex;

  mat3 norm_matrix_3;
  norm_matrix_3[0] = norm_matrix[0].xyz;
  norm_matrix_3[1] = norm_matrix[1].xyz;
  norm_matrix_3[2] = norm_matrix[2].xyz;
  fN = norm_matrix_3* normals;

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
