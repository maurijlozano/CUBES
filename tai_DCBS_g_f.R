#!/usr/bin/env Rscript
# run with Rscript --vanilla scriptName.R
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (Input file, Iterations).\n", call.=FALSE)
} 

library('tAI')
library(doParallel)

setwd(args[1])

corder = read.table('../codOrder.txt', header=F, sep=" ")
trna = read.table('trna.txt', header=T, row.names = NULL, sep=" ")
trna <- trna[match(corder$V2,trna$codon),]
dcbs <- read.table("GNM_DCBS.txt", header=T, row.names = 1, sep=" ")
dcbs <- dcbs[,c(65)]

hdata <- read.table("../freq_head.txt", header=FALSE, sep=" ")
f.counts <- read.table("genome.freq.counts", header=FALSE, row.names = 1, sep=" ")
colnames(f.counts) <- c(as.character(unlist(hdata[,2:60])), "AUG", "UGG")
corder2 <- corder[-c(11,12,15),]
t<-t(f.counts)
t <- t[match(corder2$V2,row.names(t)),]
f.o <- t(t)

fi <- as.data.frame(cbind(f.o, DCBS = dcbs))
fi$DCBS <- as.numeric(as.character(fi$DCBS))

if (trna$trna[35] == 0){sk <- 1 } else {sk <- 0}

tai_DCBS_corr <- function(sv, tRNA, sking, DCBS)
{
fixed.s=c(0,0,0,0)
s=append(fixed.s,sv)

if (min(s)>=0 & max(s)<=1){
    ws <- get.ws(tRNA=tRNA, s=s, sking=sk)
   # cat("new data iteration, s = ",s, "\n")
    tai <- get.tai(f.o[,-33], ws)
    tai <- as.numeric(tai)
    ts <- cor(tai, DCBS, method="spearman")
    return(1-ts)
    }
else{
    return(1000)
    }
}
#s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)
#remember that the function returns 1-ts
fixed.s=c(0,0,0,0)

#merge function

merge_data <- function(a,b){
    rbind(a, b)
}


#parallel
#system.time(block)
ncores <- detectCores()
registerDoParallel(cores=ncores)
#getDoParWorkers()

s.values = foreach(i=1:args[2], .export=c('tai_DCBS_corr'), .packages='tAI', .combine=merge_data) %dopar% {
    sv <- runif(5, 0, 1)
    o <- optim(par = sv, fn = tai_DCBS_corr, tRNA = trna$trna, sking = 1, DCBS = fi$DCBS, method = "Nelder-Mead", control = list(maxit=10000))
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
tai_DCBS_corr(sv, trna$trna, sking=sk, fi$DCBS)

#write tables
data <- cbind(tai, fi$DCBS)
TS_table <- s.values[order(s.values[,10], decreasing=T),]
write.table(TS_table, file = "s_opts_DCBS_GNM.txt", sep = "\t", dec = ".", row.names = TRUE)
write.table(data, file = "tai_vs_DCBS_GNM.txt", sep = "\t", dec = ".", row.names = TRUE)
#round(o$par, 2) if Sopt values rounded to 2 dec, correlation changes...


library('ggplot2')
library('ggrepel')
library('ggthemes')
library('ggpmisc')

data <- as.data.frame(cbind(tai,DCBS=fi$DCBS))
my.formula <- y ~ x

svg(filename="Tai_vs_DCBS_GNM.svg", width=8, height=6, pointsize=10)

ggplot(data, aes(x=DCBS, y=tai)) +
geom_point(show.legend = FALSE) +
geom_smooth(method=lm) +
stat_poly_eq(formula = my.formula, eq.with.lhs = "italic(hat(y))~`=`~",
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE) +
ggtitle("Adaptation Index tAI vs DCBS") + 
xlab("DCBS") +
ylab("tAI") + 
theme_calc() +
theme(axis.text.x = element_text(size=10),
	  axis.text.y = element_text(size=10),
	  axis.title.x = element_text(size=12, face="bold"),
	  axis.title.y = element_text(size=12, face="bold")
) 

dev.off()




