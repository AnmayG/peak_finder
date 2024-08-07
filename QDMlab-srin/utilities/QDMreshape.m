function data = QDMreshape(dataStack, imgNumRows, imgNumCols) 
%[data] = QDMreshape(dataStack, imgNumRows, imgNumCols)
% reshapes the 51xnumber_of_pixel data from the QDM into Y x X x freq
    data = reshape(dataStack, [], imgNumCols, imgNumRows); % reshape into freq x col x rows
    data = permute(data,[3 2 1]); % permute the axis to rows x cols x freq
end