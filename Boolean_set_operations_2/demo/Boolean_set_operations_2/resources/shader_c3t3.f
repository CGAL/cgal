#version 120
varying highp vec4 color;
varying highp vec4 fP; 
varying highp vec3 fN; 
uniform highp vec4 light_pos;  
uniform highp vec4 light_diff; 
uniform highp vec4 light_spec; 
uniform highp vec4 light_amb;  
uniform highp float spec_power ;
uniform int is_two_side; 
uniform bool is_selected;
void main(void) {
 if(color.w<0)
 {
    vec4 my_color = vec4(color.xyz, 1.);
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
        diffuse = abs(dot(N,L)) * light_diff * my_color; 
    else 
        diffuse = max(dot(N,L), 0.0) * light_diff * my_color; 
    highp vec4 specular = pow(max(dot(R,V), 0.0), spec_power) * light_spec; 
    vec4 ret_color = vec4((my_color*light_amb).xyz + diffuse.xyz + specular.xyz,1); 
    if(is_selected)
        gl_FragColor = vec4(ret_color.r+70.0/255.0, ret_color.g+70.0/255.0, ret_color.b+70.0/255.0, 1.0);
    else
        gl_FragColor = ret_color;
 }
 else
   discard;
}


