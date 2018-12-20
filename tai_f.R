#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 


library('tAI')
library(doParallel)

setwd(args[1])

corder = read.table('../codOrder.txt', header=F, sep=" ")
trna = read.table('trna.txt', header=T, row.names = NULL, sep=" ")
trna <- trna[match(corder$V2,trna$codon),]
hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
f.counts <- read.table("freq.counts", header=FALSE, row.names = 1, sep=" ")
colnames(f.counts) <- c(as.character(unlist(hdata[,2:60])), "AUG", "UGG")
corder2 <- corder[-c(11,12,15),]
t<-t(f.counts)
t <- t[match(corder2$V2,row.names(t)),]
f.o <- t(t)


indexes <- read.table("C1.table", header=T, row.names = NULL, sep="\t")

fi <- as.data.frame(cbind(f.o,indexes))
if (length(grep("[*]",indexes$Nc)) > 0 ){fi <- fi[-grep("[*]",indexes$Nc),]}

f.o <- fi[,1:61]
fi <- fi[,62:66]
fi$Nc <- as.numeric(as.character(fi$Nc))
fi$GC3s <- as.numeric(as.character(fi$GC3s))

if (trna$trna[35] == 0){sk <- 1 } else {sk <- 0}
tai_nc_corr <- function(sv, tRNA, sking, Nc, GC3s)
{
fixed.s=c(0,0,0,0)
s=append(fixed.s,sv)

if (min(s)>=0 & max(s)<=1){
    ws <- get.ws(tRNA=tRNA, s=s, sking=sk)
   # cat("new data iteration, s = ",s, "\n")
    tai <- get.tai(f.o[,-33], ws)
    tai <- as.numeric(tai)
    ts <- get.s(tai, fi$Nc, fi$GC3s)
    return(1-ts)
    }
else{
    return(1000)
    }
}

#s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)
#remember that the function returns 1-ts
fixed.s=c(0,0,0,0)



#single thread

#s.values <- NULL
#for (i in 1:2){
#sv <- runif(5, 0, 1)
#o <- optim(par = sv, fn = tai_nc_corr, tRNA = trna$trna, sking = 1, Nc = indexes$Nc, GC3s = indexes$GC3s, method = "Nelder-Mead", control = list(maxit=10000))
#svo <- round(o$par, 2)
#so <- t(as.data.frame(append(fixed.s,svo)))
#colnames(so) <- c(1:9)
#rownames(so) <- NULL
#so <- cbind(so, TScor=(1-o$value))
#s.values <- rbind(s.values, so)
#print(s.values[i,])
#}


#merge function

merge_data <- function(a,b){
    rbind(a, b)
}


#parallel
#system.time(block)
ncores <- detectCores()
registerDoParallel(cores=ncores)
#getDoParWorkers()

s.values = foreach(i=1:args[2], .export=c('tai_nc_corr'), .packages='tAI', .combine=merge_data) %dopar% {
    sv <- runif(5, 0, 1)
    o <- optim(par = sv, fn = tai_nc_corr, tRNA = trna$trna, sking = 1, Nc = indexes$Nc, GC3s = indexes$GC3s, method = "Nelder-Mead", control = list(maxit=10000))
    so <- t(as.data.frame(append(fixed.s,o$par)))
    colnames(so) <- c(1:9)
    rownames(so) <- NULL
    so <- cbind(so, TScor=(1-o$value))
}

#results summary and graphics

Sopt <- s.values[which(s.values[,10] == (max(s.values[,10]))),1:9]
ws <- get.ws(tRNA=trna$trna,s=Sopt, sking=sk)
tai <- get.tai(f.o[,-33], ws)
#hist(tai)
sv=Sopt[5:9]
tai_nc_corr(sv, trna$trna, sking=sk, fi$Nc, fi$GC3s)

nc.aj <- nc.adj(fi$Nc, fi$GC3s)
#plot(tai,nc.aj)
data <- cbind(tai,nc.aj)

TS_table <- s.values[order(s.values[,10], decreasing=T),]

write.table(TS_table, file = "s_opts.txt", sep = "\t", dec = ".", row.names = TRUE)
write.table(data, file = "tai_vs_ncaj.txt", sep = "\t", dec = ".", row.names = TRUE)
#round(o$par, 2) if Sopt values rounded to 2 dec, correlation changes...


library('ggplot2')
library('ggrepel')
library('ggthemes')
library('ggpmisc')

data <- as.data.frame(cbind(tai,nc.aj))
my.formula <- y ~ x

svg(filename="tai_vs_NcAj.svg", width=8, height=6, pointsize=10)

ggplot(data, aes(x=nc.aj, y=tai)) +
geom_point(show.legend = FALSE) +
geom_smooth(method=lm) +
stat_poly_eq(formula = my.formula, eq.with.lhs = "italic(hat(y))~`=`~",
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE) +
ggtitle("Adaptation Index tAI vs Nc.Adjusted") + 
xlab("Nc.Aj") +
ylab("tAI") + 
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()




