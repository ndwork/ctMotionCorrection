function [x, metric_normgrad, metric_laplacian, metric_histogram, ...
  metric_variance] = makeMetricImages()
close all;
    img0 = phantom();
    y = 0;  
    x = -5:.1:5;
    metric_normgrad = zeros(1,numel(x));
    metric_laplacian = zeros(1,numel(x));
    metric_histogram = zeros(1,numel(x));
    metric_variance = zeros(1,numel(x));
    %figure;
    for i = 1:numel(x)
       img1 = translateImg(img0,[0,x(i)]);
       img2 = translateImg(img0,[0,2*x(i)]);
       img3 = translateImg(img0,[0,3*x(i)]);
       img4 = translateImg(img0,[0,x(i)/2]);
       img5 = translateImg(img0,[0,x(i)/3]);
       img6 = translateImg(img0,[0,x(i)/4]);
       img7 = translateImg(img0,[0,x(i)/5]);
       img8 = translateImg(img0,[0,x(i)/6]);

       img = img0+img1+img2+img3+img4+img5+img6+img7+img8;
       
       metric_normgrad(i) = normgrad(img);
       metric_laplacian(i) = laplacian(img);
       metric_histogram(i) = histogram(img);
       metric_variance(i) = variance(img);       
       %clf;
       %imshow(img,[]); 
    end
    
%     lw = 2;
%     fs = 16; % font size for labels
%     ts = 14; % font size for tick marks
%     
%     figure;
%     set(gca,'FontSize',ts)
%     plot(x,metric_normgrad,'LineWidth',lw);
%     xlabel('Image Shift (pixels)','FontSize',fs)
%     ylabel('Normalized Gradient Squared','FontSize',fs)
%     set(gca,'XTick',linspace(-5,5,3))
%     L = get(gca,'YLim');
%     set(gca,'YTick',linspace(L(1),L(2),3))
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     
%     figure;
%     set(gca,'FontSize',ts)
%     plot(x,metric_laplacian,'LineWidth',lw);
%     xlabel('Image Shift (pixels)','FontSize',fs)
%     ylabel('Laplacian','FontSize',fs)
%     set(gca,'XTick',linspace(-5,5,3))
%     L = get(gca,'YLim');
%     set(gca,'YTick',linspace(L(1),L(2),3))
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     
%     figure;
%     set(gca,'FontSize',ts)
%     plot(x,metric_histogram,'LineWidth',lw);
%     xlabel('Image Shift (pixels)','FontSize',fs)
%     ylabel('Histogram Energy','FontSize',fs)
%     set(gca,'XTick',linspace(-5,5,3))
%     L = get(gca,'YLim');
%     set(gca,'YTick',linspace(L(1),L(2),3))
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     
%     figure;
%     set(gca,'FontSize',ts)
%     plot(x,metric_variance,'LineWidth',lw); 
%     xlabel('Image Shift (pixels)','FontSize',fs)
%     ylabel('Variance','FontSize',fs)
%     set(gca,'XTick',linspace(-5,5,3))
%     L = get(gca,'YLim');
%     set(gca,'YTick',linspace(L(1),L(2),3))
%     set(gca,'XMinorTick','on','YMinorTick','on')
end

function [out] = normgrad(img)
    imgv = img(:);
    out = sum((abs(conv([1;-1],imgv))./sum(abs(conv([1;-1],imgv)))).^2);    
end

function [out] = laplacian(img)
    out = sum(sum(abs(conv2([0,-1,0;-1,4,-1;0,-1,0],img))));    
end

function [out] = histogram(img)
    imgv = img(:);
    count = hist(imgv);
    out = sum(abs(count).^2);    
end

function [out] = variance(img)
    imgv = img(:);
    out = sum(var(imgv));    
end