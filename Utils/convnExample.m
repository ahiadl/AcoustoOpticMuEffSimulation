

dx = 1;

spaceMat = zeros(101,101,201);
xVec = -10:dx:190;
yVec = -50:dx:50;
zVec = -50:dx:50;

[Z,Y,X] = meshgrid(yVec, zVec, xVec);
rMat = sqrt(Z.^2+Y.^2+X.^2);
spaceMat(rMat<10) = rMat(rMat<10);

kernel = zeros(21,21,201);
kernel(11,11,1:10) = 1:-0.1:0.1;
kernel(11,11,end) = 1;

figure()
subplot(2,2,2)
imagesc(squeeze(rMat(:,51,:)))
title("rMat")
colorbar
subplot(2,2,3)
imagesc(squeeze(kernel(:,11,:)))
title("Kernel")
colorbar
subplot(2,2,1)
imagesc(squeeze(spaceMat(:,51,:)))
title("Space")
colorbar


A = convn(spaceMat, kernel, 'valid');

figure();
imagesc(squeeze(A(:,61,:)))

figure();
imagesc(squeeze(A))
%%

