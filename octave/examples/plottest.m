     x = -10:0.1:10;
     plot (x, sin (x), x, cos(x));
     title ("sin(x) for x = -10:0.1:10");
     xlabel ("x");
     ylabel ("f (x)");
     text (0, sin(0), "sin");
     text (1, cos(1), "cos");
     legend ("sin ..(x)","cos");
