function [images] = makeImagesForReport(img,rotations,translations,numFrames)
%rotations: max rotation in degrees
%translations: [max x translation, max y translation]

thetas = [0:180];
lines = zeros(size(img));

lines(size(img,1)-300:size(img,1)-50,100:102) = 1;
lines(size(img,1)-300:size(img,1)-50,125:127) = 1;
lines(size(img,1)-300:size(img,1)-50,150:152) = 1;
lines(size(img,1)-300:size(img,1)-50,175:177) = 1;
lines(size(img,1)-300:size(img,1)-50,200:202) = 1;
lines(size(img,1)-300:size(img,1)-50,225:227) = 1;
lines(size(img,1)-300:size(img,1)-50,250:252) = 1;

r = linspace(0,rotations,numel(thetas));
x = linspace(-translations(1),translations(1),numel(thetas));
y = linspace(-translations(2),translations(2),numel(thetas));

frames = linspace(1,numel(thetas),numFrames);
frames = [floor(frames) 0];
k = 1;

movie = figure;
images = zeros(size(img,1), size(img,2), numFrames);
for t = 1:numel(thetas)
    line1 = imrotate(lines,thetas(t),'crop');
    img1 = imrotate(img,r(t),'crop');
    img2 = translateImg(img1,[y(t),x(t)]); 
    img3 = img2+line1;
    index = find(img3>1);
    img3(index) = 1;
    figure(movie); 
    clf; 
    imshow(img3,[]);
    if t == frames(k)
      images(:,:,k) = img3;
      k = k + 1;
    end
end

end