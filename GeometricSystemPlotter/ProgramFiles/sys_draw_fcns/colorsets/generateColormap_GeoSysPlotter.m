% Colormap for surfaces and other "intense color" plots
CUB_contour = flipud([79.6875,49.8015,31.7156;
          159.3750,99.6030,63.4312;
          239.0625,149.4045,95.1469;
          255.0000,199.2060,126.8625;
            247,247,247;
            204,204,204;
            150,150,150;
            99,99,99;
            37,37,37])*(1/255);

CUB_contour = interp1(linspace(0,100,size(CUB_contour,1)), CUB_contour, linspace(0,100,251));

% Colormap for contour plots
CUB = flipud([79.6875,49.8015,31.7156;
          159.3750,99.6030,63.4312;
          239.0625,149.4045,95.1469;
          255.0000,199.2060,126.8625;
            204,204,204;
            150,150,150;
            99,99,99;
            37,37,37])*(1/255);

CUB = interp1(linspace(0,100,size(CUB,1)), CUB, linspace(0,100,251));

% Save these as separate files to call in the 'color_**..' file:
save('CUBColormap.mat','CUB');
save('CUB_contourColormap.mat','CUB_contour');




