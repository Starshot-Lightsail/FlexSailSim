%colormapexpander

mycmp = bjet;

newcmp = [];
ef = 3;  %these are powers of two dumbass... 1 = 2x, 2 = 4x, 3 = 8x...
efn = 2^ef;
firstpart = ceil((efn-1)/2);
secondpart = efn-firstpart;
for n=1:firstpart
    newcmp(end+1,:) = mycmp(1,:);
end
for n=1:(length(mycmp)-1)
    for nn = 1:efn
        newcmp(end+1,:) = (nn/efn)*(mycmp(n+1,:)) + ((efn-nn)/efn)*(mycmp(n,:));
    end
end

for n=1:secondpart
    newcmp(end+1,:) = mycmp(end,:);
end
