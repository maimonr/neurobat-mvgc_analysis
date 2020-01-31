function idx = slidingWin(N, winSize, overlap)

nWins = floor((N)/(winSize - overlap));
idx = zeros(nWins,winSize);

for w = 1:nWins
   idx(w,:) = (winSize-overlap)*(w-1) +( 1:winSize); 
end

badWins = find(idx>N);
[I,~] = ind2sub(size(idx),badWins);
idx(unique(I),:) = [];


end