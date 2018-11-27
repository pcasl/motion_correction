% simulation of the motion artifacts

clear
clc
close all

N = 128;
Nframe = 50;

motion(:,1) = sin(2*pi*2*(1:Nframe)/Nframe) + 0.2*sin(2*pi*4*(1:Nframe)/Nframe); % displacement in x
motion(:,2) = cos(2*pi*2*(1:Nframe)/Nframe) + 0.3*cos(2*pi*3*(1:Nframe)/Nframe); % displacement in y
motion(:,3) = cos(2*pi*(1:Nframe)/Nframe) + (1:Nframe)/Nframe;

% displacement 
max_displacement = 20;
I_square = ones(round(N/2));
[x_length,y_length] = size(I_square);
x_orgin = round(N/2);
y_orgin = round(N/2);
for frame = 1:Nframe
    x_start = round(x_orgin + max_displacement*motion(frame,1));
    y_start = round(y_orgin + max_displacement*motion(frame,2));
    x_end = x_start+x_length-1;
    y_end = y_start+y_length-1;
    
    I1(x_start:x_end,y_start:y_end,frame) = I_square;
end

% room-in
diameter = 0.9*(smooth(motion(:,1))-min(motion(:,1))+2) ./ (max(motion(:,1))-min(motion(:,1))+2);
for frame = 1:50
    I2(:,:,frame) = CircleMask(N,diameter(frame));
end

% rotation
