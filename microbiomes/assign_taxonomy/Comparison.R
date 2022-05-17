library("arsenal")
library("stringr")
library("dplyr")

taxa1 <- read.table("file.path", 
                    sep = "\t" , 
                    header = T,
                    row.names = 1)

taxa2 <- read.table("file.path",
                    sep = "\t",
                    header = T,
                    row.names = 1)

taxa1a <- taxa1[,c('Class',
                'Order',
                'Family',
                'Genus',
                'Species')]

taxa2a <- taxa2[c('class',
                'order',
                'family',
                'genus',
                'species')]




s <- summary(comparedf(x = taxa1a, y = taxa2a, 
                  control = comparedf.control(tol.vars = "case")))

capture.output(summ, file = "file.path")



