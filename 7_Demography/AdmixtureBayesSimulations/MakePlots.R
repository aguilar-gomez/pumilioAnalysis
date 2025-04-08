pdf(file = "Plot.pdf",     width = 13*0.65,  height = 7*0.65)

A = read.table("~/desktop/PumilioSimulations/0events/AdBayesTop.txt", header = F)
B = read.table("~/desktop/PumilioSimulations/1events/AdBayesTop.txt", header = F)

df = data.frame(A[[1]] , A[[2]] , B[[1]], B[[2]])
names(df) = c("0 admixtures (Mode)", "0 admixtures (Mean)",  
              "1 admixture (Mode)", "1 admixture (Mean)")

boxplot(df, ylab = "Topology Equality",   ylim = c(0,1.01), range=0,  
        col = c("red", "green"  ,"steelblue3" , "lightskyblue"))

dev.off()
