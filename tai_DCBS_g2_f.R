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

get.ws <- function (tRNA, s = NULL, sking) 
{
    if (is.null(s)) 
        s <- c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68, 0.89, 0.5)
    p = 1 - s
    W = NULL
	box4 <- c(5,17,21,29,33,37,49,53,61)
	box2 <- c(1,9,13,25,41,45,57)
	codonWorder <- c(1,2,3,4,9,10,11,12,13,14,15,16,25,26,27,28,41,42,43,44,45,46,47,48,57,58,59,60,5,6,7,8,17,18,19,20,21,22,23,24,29,30,31,32,33,34,35,36,37,38,39,40,49,50,51,52,53,54,55,56,61,62,63,64)
   
	for (i in box2) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])
	for (i in box4) W = c(W, p[1] * tRNA[i] + p[5] * tRNA[i + 1] + p[10] * tRNA[i + 2], p[2] * tRNA[i + 1] + p[6] * tRNA[i], p[3] * tRNA[i + 2] + p[7] * tRNA[i], p[4] * tRNA[i + 3] + p[8] * tRNA[i + 2])

	W <- W[match(1:64,codonWorder)]
	W[36] = p[4] * tRNA[36]
    if (sking == 1) 
		W[35] = p[9]
    W = W[-c(11, 12, 15, 36)]
    w = W/max(W)
    if (sum(w == 0) > 0) {
        ws <- w[w != 0]
        gm <- exp(sum(log(ws))/length(ws))
        w[w == 0] = gm
    }
    return(w)
}



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
    sv <- runif(6, 0, 1)
    o <- optim(par = sv, fn = tai_DCBS_corr, tRNA = trna$trna, sking = 1, DCBS = fi$DCBS, method = "Nelder-Mead", control = list(maxit=10000))
    so <- t(as.data.frame(append(fixed.s,o$par)))
    colnames(so) <- c(1:10)
    rownames(so) <- NULL
    so <- cbind(so, TScor=(1-o$value))
}

#results summary and graphics

Sopt <- s.values[which(s.values[,11] == (max(s.values[,11]))),1:10]

ws <- get.ws(tRNA=trna$trna,s=Sopt, sking=sk)
tai <- get.tai(f.o[,-33], ws)
#hist(tai)
sv=Sopt[5:10]
tai_DCBS_corr(sv, trna$trna, sking=sk, fi$DCBS)

#write tables
data <- cbind(tai, fi$DCBS)
TS_table <- s.values[order(s.values[,11], decreasing=T),]
write.table(TS_table, file = "s_opts_DCBS_GNM2.txt", sep = "\t", dec = ".", row.names = TRUE)
write.table(data, file = "tai_vs_DCBS_GNM2.txt", sep = "\t", dec = ".", row.names = TRUE)
#round(o$par, 2) if Sopt values rounded to 2 dec, correlation changes...


library('ggplot2')
library('ggrepel')
library('ggthemes')
library('ggpmisc')

data <- as.data.frame(cbind(tai,DCBS=fi$DCBS))
my.formula <- y ~ x

svg(filename="Tai_vs_DCBS_GNM2.svg", width=8, height=6, pointsize=10)

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




