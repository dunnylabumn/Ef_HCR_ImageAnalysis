idx = ( [I_cells.Area] < 30);

BW2 = ismember(L,find(idx));
BW3 = ismember(L,find(~idx));

figure
subplot(131), imshow(I)
subplot(132), imshow(BW2)
subplot(133), imshow(BW3)

