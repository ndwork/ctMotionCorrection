function makeBlurryPhantoms
    img0 = phantom();
    y = 0;
    horizontaltrans = [5,1,0];
    for i = horizontaltrans
       figure;
       x = i;
       img1 = translateImg(img0,[0,x]);
       img2 = translateImg(img0,[0,2*x]);
       img3 = translateImg(img0,[0,3*x]);
       img4 = translateImg(img0,[0,x/2]);
       img5 = translateImg(img0,[0,x/3]);
       img6 = translateImg(img0,[0,x/4]);
       img7 = translateImg(img0,[0,x/5]);
       img8 = translateImg(img0,[0,x/6]);

       img = img0+img1+img2+img3+img4+img5+img6+img7+img8;
       clf;
       imshow(img,[]); 
    end
end