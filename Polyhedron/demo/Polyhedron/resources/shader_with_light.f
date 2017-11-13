#version 120
varying highp vec4 color;
varying highp vec4 fP; 
varying highp vec3 fN;
varying highp float dist[6];
uniform highp vec4 light_pos;  
uniform highp vec4 light_diff; 
uniform highp vec4 light_spec; 
uniform highp vec4 light_amb;
uniform highp float spec_power ; 
uniform int is_two_side; 
uniform bool is_selected;
uniform bool is_clipbox_on;
void main(void) {

if(is_clipbox_on)
  if(dist[0]>0 ||
     dist[1]>0 ||
     dist[2]>0 ||
     dist[3]>0 ||
     dist[4]>0 ||
     dist[5]>0)
    discard;

   highp vec3 L = light_pos.xyz - fP.xyz;
   highp vec3 V = -fP.xyz;
   highp vec3 N;
   if(fN == highp vec3(0.0,0.0,0.0))
       N = highp vec3(0.0,0.0,0.0);
   else
       N = normalize(fN);
   L = normalize(L);
   V = normalize(V);
   highp vec3 R = reflect(-L, N);
  vec4 diffuse;
   if(is_two_side == 1)
       diffuse = abs(dot(N,L)) * light_diff * color;
   else
       diffuse = max(dot(N,L), 0.0) * light_diff * color;
   highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec;
   vec4 ret_color = vec4((color*light_amb).xyz + diffuse.xyz + specular.xyz,1);
   if(is_selected)
       gl_FragColor = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, 1.0);
   else
       gl_FragColor = ret_color;
}


