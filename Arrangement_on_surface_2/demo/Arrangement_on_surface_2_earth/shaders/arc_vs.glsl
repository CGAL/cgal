
#version 330

layout (location = 0) in vec3 a_pos;

uniform mat4 u_mvp;

out vec3 v_pos;

void main() {
  v_pos = a_pos;
  gl_Position = u_mvp * vec4(a_pos.xyz, 1);
}
