% set dimensions of matrix to height and width by filling with zeros

function newmatrix=fillzeros(matrix, height, width)

if size(matrix,1) > height || size(matrix,2) > width
    fprintf('[error - fillzeros] dimensions do not match')
end

newmatrix = zeros(height, width);
newmatrix(1:size(matrix,1),1:size(matrix,2))=matrix;

    