df <- read.csv("../oldriemann/out/gzetaE28/values.csv", header = FALSE);
matplot(df[1],df[2],type = "b", tck=1)