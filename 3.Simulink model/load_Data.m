ref_pos   = [out.tout, out.JointAngle];
ref_vel   = [out.tout, out.JointVel];
ref_accel = [out.tout, out.JointAccel];
% Set the filter parameters
window_size = 51; 
poly_order = 3; 

% 1. Smooth the raw data
smooth_pos = sgolayfilt(out.JointAngle, poly_order, window_size);
smooth_vel = sgolayfilt(out.JointVel, poly_order, window_size);
smooth_acc = sgolayfilt(out.JointAccel, poly_order, window_size);

% 2. Stitch the time vector and save as BRAND NEW variables
filtered_pos   = [out.tout, smooth_pos];
filtered_vel   = [out.tout, smooth_vel];
filtered_accel = [out.tout, smooth_acc];