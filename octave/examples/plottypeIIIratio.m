     clf;
     x= -0.2:0.1:0.2;
     fid = fopen('data/typeIIIratio.txt');
     [ratios] = fscanf(fid, '&%f ',[5,8])'
     hold on;
     for i = 1:2:7
        plot (x, ratios (i,:), x, ratios (i+1,:));
     endfor
     #title ("sin(x) for x = -10:0.1:10");
     xlabel ("displacement/gram increment");
     ylabel ("Ratio Type IIIL/Type IIIR");
     text (0.15, 1.25, "l=4");
     text (0.15, 2.0, "l=5");
     text (0.15, 2.5, "l=6");
     text (0.17, 4.7, "l=7");
     text (0.15, 4.8, "l=8");
     text (0.085, ratios (6,4), "l=9");
     text (0.0525, 4.0, "l=10");
     text (0.025, 4.25, "l=11");
     # legend ("sin ..(x)","cos");
     axis([x(1,1),x(1,5), 0, 7]);
     hold off;
 