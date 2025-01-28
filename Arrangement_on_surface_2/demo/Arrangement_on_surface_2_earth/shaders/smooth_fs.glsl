
#version 330

uniform vec4 u_color;
uniform vec4 u_plane;

in vec3 v_pos;
in vec3 v_normal;

out vec4 out_color;

void main() {
  if (dot(u_plane, vec4(v_pos, 1)) < 0) discard;

  const vec3 lightDir = normalize(vec3(0,0,-1));

  //float c = clamp(dot(lightDir,triNormal), 0, 1);
  vec3 n = normalize(v_normal);
  float c = abs(dot(lightDir, n) );
  out_color = u_color * (vec4(.2, .2,.2,1) + 0.8 * vec4(c,c,c,1));

  //color = vec4(1,1,0,1);
  //color = vCol;
}
