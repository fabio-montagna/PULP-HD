function [M32] = compress_hypervectors(M)
%
% DESCRIPTION   : to compress a matrix/vector into a matrix/vector composed
%                 by 32-bit unsigned integer variables  
%
% INPUTS:
%   M           : input matrix/vector
% OUTPUTS:
%   M32         : compressed matrix/vector
%    
    [r_M,c_M] = size(M);
    
    dim = floor(c_M/32);
    if (mod(dim, 32)) 
        dim=dim+1;
    end
    
    M32 = [];
 
    temp = uint32(0); 
  
    for z = 1 : r_M
        for j = 1 : dim - 1
            for i = 1 : 32
                temp = uint32(temp + uint32(bitsll(M(z, (i + (32 * (j - 1)))), (32 - i))));
            end  
            M32(z, j) = temp;
            temp = 0; 
        end   
        for i = 1 : 16
            temp = temp + uint32(bitsll(M(z, ((i + (32 * (dim - 1))))), (32 - i)));
        end
        M32(z, dim) = temp;
        temp = 0; 
    end
    
end