%%%%%%%%%%%%%%%%%%%%
%  	 Ryan Heyser   %    
%  Model Rendering %
%%%%%%%%%%%%%%%%%%%%               

%%%%%%%%%%%%  Initial Setup  %%%%%%%%%%%%%%

%refresh to initial state
clear;
clc;
clf;

%lightsource [x y z]
ls_pos = [10 50 1];

%lightsource [r g b]
ls_rgb = [100 50 100];

%field of view or aspect ratio
fov = 70;

%ambient light RGB
am_rgb = [50 50 50];
%emissive material RGB
em_rgb = [200 255 100];
%diffuse material RGB
dm_rgb = [150 140 200];
%specular material RGB
sm_rgb = [240 240 255];
sp_pow = 6; %specular intensity
%material RGB
material_rgb = [20 100 110];

%object's world position
world_obj_coord = [0 0 0];

%camera location
camera = [15 0 0];

%camera lookat
lookat = [0 0 0];

%camera's up
up = [0 1 0];

%aspect ratio
aspect = 1;

%setup axis
axis([-1 1 -1 1]);

%world-space position and orientation of the object 
%represented by rotations around x,y,z axis in degrees
%object = [45 45 90] where each value is the degree
%of rotation
obj_rot = [0 45 10];

%%%%%%%%%%%%  Read Model  %%%%%%%%%%%%%%

%read txt file
object = load('shuttle.txt');

%calculate min and max of object
min = min(min(object));
max = max(max(object));
diff = max - min;

%use min and max to 
%zn = abs(min)/2 + norm(lookat)/2;
%zf = abs(max) + 2*norm(world_obj_coord);
zn=1;
zf=20;

%setup
xdata = [];
ydata = [];
cdata = [];
zdata = [];
%ztemp = [];

%%%%%%%%%%%%  Transformations  %%%%%%%%%%%%%%
for index = 1:length(object)
    triangle = object(index,:);
    subvert1 = [triangle(1,1:3) 1];
    subvert2 = [triangle(1,4:6) 1];
    subvert3 = [triangle(1,7:9) 1];

   %object transformations 

%rotation !!!!!!!!!!!!USING LEFT HANDED SYSTEM!!!!!!!!!!!!
%rotate z
    rotz_mat1 = subvert1*[cosd(obj_rot(3)), sind(obj_rot(3)), 0, 0; ...
        -sind(obj_rot(3)), cosd(obj_rot(3)), 0, 0; 0,0,1,0; 0,0,0,1];
    rotz_mat2 = subvert2*[cosd(obj_rot(3)), sind(obj_rot(3)), 0, 0; ...
        -sind(obj_rot(3)), cosd(obj_rot(3)), 0, 0; 0,0,1,0; 0,0,0,1];
    rotz_mat3 = subvert3*[cosd(obj_rot(3)), sind(obj_rot(3)), 0, 0; ...
        -sind(obj_rot(3)), cosd(obj_rot(3)), 0, 0; 0,0,1,0; 0,0,0,1];

%rotate y
    roty_mat1 = rotz_mat1*[cosd(obj_rot(2)),0,-sind(obj_rot(2)),0; ...
        0,1,0,0; sind(obj_rot(2)),0,cosd(obj_rot(2)),0; 0,0,0,1];
    roty_mat2 = rotz_mat2*[cosd(obj_rot(2)),0,-sind(obj_rot(2)),0; ...
        0,1,0,0; sind(obj_rot(2)),0,cosd(obj_rot(2)),0; 0,0,0,1];
    roty_mat3 = rotz_mat3*[cosd(obj_rot(2)),0,-sind(obj_rot(2)),0; ...
        0,1,0,0; sind(obj_rot(2)),0,cosd(obj_rot(2)),0; 0,0,0,1];

%rotate x
    rotx_mat1 = roty_mat1*[1,0,0,0; 0, cosd(obj_rot(1)),sind(obj_rot(1)),0; ...
        0,-sind(obj_rot(1)),cosd(obj_rot(1)),0; 0,0,0,1];
    rotx_mat2 = roty_mat2*[1,0,0,0; 0, cosd(obj_rot(1)),sind(obj_rot(1)),0; ...
        0,-sind(obj_rot(1)),cosd(obj_rot(1)),0; 0,0,0,1];
    rotx_mat3 = roty_mat3*[1,0,0,0; 0, cosd(obj_rot(1)),sind(obj_rot(1)),0; ...
        0,-sind(obj_rot(1)),cosd(obj_rot(1)),0; 0,0,0,1];

    
%apply world transformation to each of the vertices 
%of the model to get appropriate position and 
%orientation in world coordinates.
    world_mat1 = [world_obj_coord - rotx_mat1(1:3) 1];
    world_mat2 = [world_obj_coord - rotx_mat2(1:3) 1];
    world_mat3 = [world_obj_coord - rotx_mat3(1:3) 1];

    
%%%%%%%%%%%%% Lighting %%%%%%%%%%%%%%%%

%ambient + diffuse + specular + emissive
    %surf_norm = 
    surf_norm = -1*cross(world_mat2(1:3)-world_mat1(1:3),world_mat3(1:3)-world_mat1(1:3));
    surf_norm = surf_norm/norm(surf_norm);
    center = [mean([world_mat1(1),world_mat2(1),world_mat3(1)]), mean([world_mat1(2),world_mat2(2),world_mat3(2)]), mean([world_mat1(3),world_mat2(3),world_mat3(3)]) ];
    E = center - camera;
    E = E/norm(E);
    
    %backface culling
    if (dot(E, surf_norm) > 0)
        %sprintf('backface culled')
        continue;
    end
    
    L = center - ls_pos;
    L = L/norm(L);
    h = (E + L)/norm(E+L);
    
    am_clr = material_rgb.*am_rgb;
    %since the version of matlab on my desktop hates the max function
    %dm_clr = max(dot(-L,surf_norm),0) * (material_rgb .* dm_rgb);
    %sp_clr = (max(dot(-surf_norm,h),0))^sp_pow * (material_rgb .* sm_rgb);
    if (dot(-L,surf_norm) > 0)
        dm_clr = dot(-L,surf_norm) * (material_rgb .* dm_rgb);
    else
        dm_clr = 0;
    end
    
    if (dot(-surf_norm,h) > 0)
        sp_clr = (dot(-surf_norm,h))^sp_pow * (material_rgb .* sm_rgb);
    else
        sp_clr = 0;
    end
    
    light_rgb = am_clr + em_rgb + dm_clr + sp_clr;     

%apply view transformation to get into eyespace
% zaxis = normal(At - Eye)
% xaxis = normal(cross(Up, zaxis))
% yaxis = cross(zaxis, xaxis)
%     
%  xaxis.x           yaxis.x           zaxis.x          0
%  xaxis.y           yaxis.y           zaxis.y          0
%  xaxis.z           yaxis.z           zaxis.z          0
% -dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1
    zaxis = (lookat - camera)/norm(lookat - camera);
    xaxis = cross(up,zaxis)/norm(cross(up,zaxis));
    yaxis = cross(zaxis,xaxis);
    viewspaceLH = [xaxis(1), yaxis(1), zaxis(1), 0; 
        xaxis(2), yaxis(2), zaxis(2), 0;
        xaxis(3),   yaxis(3),   zaxis(3),   0; 
        dot(-xaxis,camera), dot(-yaxis,camera), dot(-zaxis,camera), 1];
    view_mat1 = world_mat1*viewspaceLH;
    view_mat2 = world_mat2*viewspaceLH;
    view_mat3 = world_mat3*viewspaceLH;
    

    
%apply projection transformation to get into
%normalized coordinates
% fovy = field of view in the y direction, in rad
% aspect = aspect ratio = 1 for this project
% zn = z-value of the near view-plane
% zf = z-value of the far view-plane
% xScale     0          0               0
% 0        yScale       0               0
% 0          0       zf/(zf-zn)         1
% 0          0       -zn*zf/(zf-zn)     0
% where:
% yScale = cot(fovY/2)
% 
% xScale = yScale / aspect ratio
    yscale = cotd(fov/2);
    xscale = yscale/aspect;
    perspectiveLH = [xscale,0,0,0;0,yscale,0,0;
        0,0,zf/(zf-zn),1; 0,0,-zn*zf/(zf-zn),0];
    pers_mat1 = view_mat1*perspectiveLH;
    pers_mat2 = view_mat2*perspectiveLH;
    pers_mat3 = view_mat3*perspectiveLH;
    
    
%divide the x,y,z coordinates by w to impliment
%the perspective effect.
    final_mat1 = pers_mat1(1:3)/pers_mat1(4);
    final_mat2 = pers_mat2(1:3)/pers_mat2(4);
    final_mat3 = pers_mat3(1:3)/pers_mat3(4);

%convert to patch
    zmean = mean([final_mat1(3), final_mat2(3), final_mat3(3)],2);
    
    %z-clipping
       
    if(zmean > (1))
        %sprintf('clipped >1 zmean = %f',zmean)
        continue;
    elseif (zmean < 0) 
        %sprintf('clipped <0 zmean = %f',zmean)
        continue;
    end
        
    xdata = [xdata; final_mat1(1), final_mat2(1), final_mat3(1)];
    ydata = [ydata; final_mat1(2), final_mat2(2), final_mat3(2)];
    zdata = [zdata; zmean];
    cdata = [cdata; light_rgb]; 
    
end

[zsort index] = sort(zdata,1,'descend');
xsort = xdata(index,:);
ysort = ydata(index,:);
csort = cdata(index,:);

xsort = rot90(xsort);
ysort = rot90(ysort);

colorfix = 81225;
csort = csort/colorfix;

for index = 1:length(csort)
    patch(xsort(:,index),ysort(:,index),csort(index,:));
    pause(.01); %show triangles being drawn
end









