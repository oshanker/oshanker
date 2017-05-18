     clf;
     x= -0.2:0.1:0.2;
     fid = fopen('data/typeIIratio.txt');
     [ratios] = fscanf(fid, '&%f ',[5,10])'
     hold on;
     for i = 1:2:9
        plot (x, ratios (i,:), x, ratios (i+1,:));
     endfor
     #title ("sin(x) for x = -10:0.1:10");
     xlabel ("displacement/gram increment");
     ylabel ("Ratio Type II/Type I");
     text (-0.1, ratios (1,2), "l=2");
     text (-0.1, ratios (2,2), "l=3");
     text (-0.1, ratios (3,2), "l=4");
     text (-0.1, ratios (4,2), "l=5");
     text (-0.095, ratios (5,2), "l=6");
     text (-0.095, ratios (6,2), "l=7");
     text (-0.0625, 4.0, "l=8");
     text (-0.045, 4.25, "l=9");
     text (-0.033, 4.65, "l=10");
     text (-0.01, 4.85, "l=11");
     # legend ("sin ..(x)","cos");
     axis([x(1,1),x(1,5), 0, 5]);
     hold off;
 