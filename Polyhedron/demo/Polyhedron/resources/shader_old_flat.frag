#version 150
in vec4 color;
in vec4 fP;
uniform vec4 light_pos;
uniform vec4 light_diff;
uniform vec4 light_spec;
uniform vec4 light_amb;
uniform float spec_power ;
uniform int is_two_side;
uniform bool is_selected;
out vec4 out_color;
void main(void) {
   vec3 L = light_pos.xyz - fP.xyz;
   vec3 V = -fP.xyz;
   vec3 N;
   vec3 X = dFdx(fP.xyz);
   vec3 Y = dFdy(fP.xyz);
   vec3 normal=normalize(cross(X,Y));

   if(normal == vec3(0.0,0.0,0.0))
       N = vec3(0.0,0.0,0.0);
   else
       N = normalize(normal);
   L = normalize(L);
   V = normalize(V);
   vec3 R = reflect(-L, N);
  vec4 diffuse;
   if(is_two_side == 1)
       diffuse = abs(dot(N,L)) * light_diff * color;
   else
       diffuse = max(dot(N,L), 0.0) * light_diff * color;
   vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;
   vec4 ret_color = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1.0);
   if(is_selected)
       out_color = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, 1.0);
   else
       out_color = ret_color;
}
