function [x,y,z]=threeDcoord(p,L,C)

%Give the 3D coordinates of points from the index p within a 3D matrix with
%L lines and C columns

z=floor((p-1)/(L*C))+1;  %%%depth position
y=floor((p-(z-1)*L*C-1)/L)+1;  %%% line position
x=p-(z-1)*L*C-(y-1)*L; %%% column position