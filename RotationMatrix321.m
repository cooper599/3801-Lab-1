function [DCM] = RotationMatrix321(angles)
mat1 = [1 0 0; 0 cos(angles(1)) sin(angles(1)); 0 -sin(angles(1)) cos(angles(1))];
mat2 = [cos(angles(2)) 0 -sin(angles(2)); 0 1 0; sin(angles(2)) 0 cos(angles(2))];
mat3 = [cos(angles(3)) sin(angles(3)) 0; -sin(angles(3)) cos(angles(3)) 0; 0 0 1];
DCM = mat1 * mat2 * mat3;
end