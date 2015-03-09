function makeImagesForReport(img,rotations,translations)
%rotations: max rotation in degrees
%translations: [max x translation, max y translation]

thetas = [0:180];
line = zeros(size(img));

line(100:250,10) = 1;

r = linspace(0,rotations,numel(thetas));
x = linspace(0,translations(1),numel(thetas));
y = linspace(0,translations(2),numel(thetas));

movie = figure;
for t = 1:numel(thetas)
    line1 = imrotate(line,thetas(t),'crop');
    img1 = imrotate(img,r(t),'crop');
    img2 = translateImg(img1,[y(t),x(t)]); 
    figure(movie); clf; imshow(img2+line1,[]);
end

end