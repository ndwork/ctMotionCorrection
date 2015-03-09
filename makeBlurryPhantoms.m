function [img] = makeBlurryPhantoms(horizontaltrans)
    img0 = phantom();
    y = 0;
%     horizontaltrans = [5,1,0];
    
    img = zeros(size(phantom,1), size(phantom,2), 3);
    
    for k = 1:numel(horizontaltrans)
    
%     for i = horizontaltrans
%        figure;
       x = horizontaltrans(k);
       img1 = translateImg(img0,[0,x]);
       img2 = translateImg(img0,[0,2*x]);
       img3 = translateImg(img0,[0,3*x]);
       img4 = translateImg(img0,[0,x/2]);
       img5 = translateImg(img0,[0,x/3]);
       img6 = translateImg(img0,[0,x/4]);
       img7 = translateImg(img0,[0,x/5]);
       img8 = translateImg(img0,[0,x/6]);

       img(:,:,k) = img0+img1+img2+img3+img4+img5+img6+img7+img8;
%        clf;
%        imshow(img(:,:,k),[]); 
       
    end
end