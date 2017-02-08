#version 120
attribute highp vec4 vertex;
attribute highp vec3 normals;
attribute highp vec3 colors;
attribute highp vec3 barycenter;
uniform highp mat4 mvp_matrix;
uniform highp mat4 mv_matrix;
uniform highp float shrink_factor;
varying highp vec4 fP;
varying highp vec3 fN;
varying highp vec4 color;
void main(void)
{
  color = vec4(colors, 1.0);
  fP = mv_matrix * vertex;
  fN = mat3(mv_matrix)* normals;
  highp mat4 transOB = mat4(1, 0, 0, 0, // first column
   0, 1, 0, 0, // second column
   0, 0, 1, 0, // third column
   barycenter.x, barycenter.y, barycenter.z, 1); // fourth column
  highp mat4 transBO = mat4(1, 0, 0, 0, // first column
    0, 1, 0, 0, // second column
    0, 0, 1, 0, // third column
    -barycenter.x, -barycenter.y, -barycenter.z, 1); // fourth column
   highp mat4 scaling = mat4(shrink_factor, 0, 0, 0,
    0, shrink_factor, 0, 0,
    0, 0, shrink_factor, 0,
    0, 0, 0, 1);
  gl_Position = mvp_matrix *transOB * scaling * transBO * vertex;
}
