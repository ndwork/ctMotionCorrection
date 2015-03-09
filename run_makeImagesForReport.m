function run_makeImagesForReport
close all;
    img = phantom();
    translations = [0.5,0]; %meters
    pixsize = 0.01;
    translations = translations/pixsize; %pixels
    rotations = -60; %degrees
    img = padImgForRadon(img,50,50,1);
    figure;
    imshow(img,[])
    makeImagesForReport(img,rotations,translations)
end