function res = sat(blah)
    thesize = size(blah);
    numrows = thesize(1);
    numcols = thesize(2);
    res(1:numrows,1:numcols) = 0;
    for rowindex = 1:numrows
        for colindex = 1:numcols
            if abs(blah(rowindex,colindex))<=1
                res(rowindex,colindex) = blah(rowindex,colindex);
            else
                res(rowindex,colindex) = sign(blah(rowindex,colindex));
            end
        end
    end
end