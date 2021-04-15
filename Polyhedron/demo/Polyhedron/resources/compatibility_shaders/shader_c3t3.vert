//#version 100
attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec3 colors;
attribute highp vec3 center;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp vec4 cutplane;
uniform highp float shrink_factor;
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 color;
uniform highp float point_size;
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

  highp mat4 transOB = mat4(1., 0., 0., 0., // first column
   0., 1., 0., 0., // second column
   0., 0., 1., 0., // third column
   center.x, center.y, center.z, 1); // fourth column
  highp mat4 transBO = mat4(1., 0., 0., 0., // first column
    0., 1., 0., 0., // second column
    0., 0., 1., 0., // third column
    -center.x, -center.y, -center.z, 1); // fourth column
   highp mat4 scaling = mat4(shrink_factor, 0., 0., 0.,
    0., shrink_factor, 0., 0.,
    0., 0., shrink_factor, 0.,
    0., 0., 0., 1.);
  gl_Position = mvp_matrix *transOB * scaling * transBO * vertex;
}
