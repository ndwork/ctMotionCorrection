function run_makeImagesForReport
close all;
    img = phantom();
    translations = [0.5,0]; %meters
    pixsize = 0.01;
    translations = translations/pixsize; %pixels
    rotations = -120; %degrees
    img = imrotate(img,-rotations/2,'crop');
    img = padImgForRadon(img,50,50,1);
%     figure;
%     imshow(img,[])
    numFrames = 4;
    images = makeImagesForReport(img,rotations,translations,numFrames);
end