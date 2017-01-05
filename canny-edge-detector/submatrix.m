function [ mat ] = submatrix( matrix, i, j, radius )
    
   dim = radius*2+1;
   mat = zeros(dim);
   
   minCol = i-radius;
   maxCol = i+radius;
   
   minRow = j-radius;
   maxRow = j+radius;
   
   cols = [minCol:maxCol];
   rows = [minRow:maxRow];
   
   mat = matrix(cols,rows);
   

end

