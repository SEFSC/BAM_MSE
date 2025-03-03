##################################################################################
###   Sample R program for testing and demonstrating the "FishGraph" collection
###     of R graphics functions for analysis of stock-assessment results.
##################################################################################

##### Start fresh ###########
rm(list=ls(all=TRUE))
graphics.off()
.SavedPlots <- NULL

library(FishGraph)

##### Read in the data from the ASCII .rdat file: #####
spp <- dget("BAM-sim.rdat") 

##### Common arguments in FishGraph functions. Convenient to define them once.
ptype=NULL #(no quotes) for no plots saved; other options: "pdf", "wmf", "eps"
dtype=TRUE  #draft type
ctype=TRUE  #color type 
########## Open a graphics device ###########
windows(width = 8, height = 8, record = TRUE)
#pdf(file="BamSim1.pdf", width = 8, height = 8)
########## Call the functions #########################

Landings.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, 
               L.units=c("mt"))

Comp.yearly.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.neff=FALSE, print.neff=TRUE,
                  compact = TRUE, print.n=TRUE)

Selectivity.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.points=T, compact=TRUE)

Comp.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, c.min=0.25, corr=FALSE)
#Comp.plots(spp, use.color=ctype, graphics.type=ptype, c.min=0.25, corr=FALSE)

Growth.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, plot.all = TRUE)

PerRec.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
             user.PR = list("SPR", "ypr.kg"), F.references=list("F35"))

Eq.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
         F.references=list("F35"), user.Eq=list("L.eq.mt", "SSB.eq", "R.eq"))

StockRec.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype,
               draw.lowess = FALSE, start.drop = 0, units.rec="number age-1 fish")

CLD.total.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, first.year = "1965",
                units.CLD.w="mt", CLD.w.references=list(NULL, "msy.mt", NULL),
                CLD.n.references=NULL, plot.proportion = TRUE)

NFZ.age.plots(spp,draft=dtype, use.color=ctype, graphics.type=ptype,
              user.plots="N.age.mdyr", start.drop=0, plot.CLD=F)

BSR.time.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, legend.pos="top",
               BSR.references=list("Bmsy", "SSBmsy", "Rmsy"))
  
Index.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, resid.analysis=T)

# ylim.temp=c(0,max(na.omit(spp$t.series[,21:22])))
# plot(spp$t.series[spp$t.series[,21]>0,1],spp$t.series[spp$t.series[,21]>0,21]
#      ,xlab="Year",ylab="Survey 1",ylim=ylim.temp,cex=2,lwd=2,col="darkblue")
# lines(spp$t.series[spp$t.series[,22]>0,1],spp$t.series[spp$t.series[,22]>0,22]
#       ,type="b",cex=1.5,pch=16,col="orange",lwd=1,lty=2)
# legend(x="bottom",legend=c("Observed","Predicted"),
#        lty=c(0,2),lwd=c(2,1),pch=c(1,16),col=c("darkblue","orange"),pt.cex=c(2,1.5))
# 
# ylim.temp=c(0,max(na.omit(spp$t.series[,24:25])))
# plot(spp$t.series[spp$t.series[,24]>0,1],spp$t.series[spp$t.series[,24]>0,24]
#      ,xlab="Year",ylab="Survey 2",ylim=ylim.temp,cex=2,lwd=2,col="darkblue")
# lines(spp$t.series[spp$t.series[,25]>0,1],spp$t.series[spp$t.series[,25]>0,25]
#       ,type="b",cex=1.5,pch=16,col="orange",lwd=1,lty=2)
# legend(x="bottom",legend=c("Observed","Predicted"),
#        lty=c(0,2),lwd=c(2,1),pch=c(1,16),col=c("darkblue","orange"),pt.cex=c(2,1.5))
# 
# ylim.temp=c(0,max(na.omit(spp$t.series[,27:28])))
# plot(spp$t.series[spp$t.series[,27]>0,1],spp$t.series[spp$t.series[,27]>0,27]
#      ,xlab="Year",ylab="Survey 3",ylim=ylim.temp,cex=2,lwd=2,col="darkblue")
# lines(spp$t.series[spp$t.series[,28]>0,1],spp$t.series[spp$t.series[,28]>0,28]
#       ,type="b",cex=1.5,pch=16,col="orange",lwd=1,lty=2)
# legend(x="bottom",legend=c("Observed","Predicted"),
#        lty=c(0,2),lwd=c(2,1),pch=c(1,16),col=c("darkblue","orange"),pt.cex=c(2,1.5))
# 
# ylim.temp=c(0,max(na.omit(spp$t.series[,30:31])))
# plot(spp$t.series[spp$t.series[,30]>0,1],spp$t.series[spp$t.series[,30]>0,30]
#      ,xlab="Year",ylab="Survey 4",ylim=ylim.temp,cex=2,lwd=2,col="darkblue")
# lines(spp$t.series[spp$t.series[,31]>0,1],spp$t.series[spp$t.series[,31]>0,31]
#       ,type="b",cex=1.5,pch=16,col="orange",lwd=1,lty=2)
# legend(x="bottom",legend=c("Observed","Predicted"),
#        lty=c(0,2),lwd=c(2,1),pch=c(1,16),col=c("darkblue","orange"),pt.cex=c(2,1.5))

F.time.plots(spp, draft=dtype, use.color=ctype, graphics.type=ptype, 
            F.references=list(a="F35"), F.additional=c("F.F35.ratio"))

Cohort.plots(spp,draft=dtype, graphics.type=ptype)

Parm.plots(spp, graphics.type=ptype)

Bound.vec.plots(spp, draft=dtype, graphics.type=ptype)

Phase.plots(spp, start.drop=0, draft=dtype, use.color=ctype, graphics.type=ptype, Xaxis.F=F, year.pos=3,F.B.references=list("Fmsy","msst"))

# windows(width = 8, height = 9, record = TRUE)
# par(mfrow=c(3,2))
# plot(spp$projection$year, spp$projection$F.proj, type="o", lwd=2, col="blue", xlab="", ylab="Projected F")
# plot(spp$projection$year, spp$projection$SSB, type="o", lwd=2, col="blue", xlab="", ylab="Projected SSB (mt)")
# plot(spp$projection$year, spp$projection$L.knum.proj, type="o", lwd=2, col="blue", xlab="Year", ylab="Landings (1000s)")
# plot(spp$projection$year, spp$projection$L.mt.proj, type="o", lwd=2, col="blue", xlab="Year", ylab="Landings (mt)")
# plot(spp$projection$year, spp$projection$Ddead.knum.proj, type="o", lwd=2, col="blue", xlab="Year", ylab="Dead discards (1000s)")
# plot(spp$projection$year, spp$projection$Ddead.mt.proj, type="o", lwd=2, col="blue", xlab="Year", ylab="Dead discards (mt)")

#dev.off()

#out=data.frame(year=spp$t.series$year, ssb=spp$t.series$SSB, recruits=spp$t.series$recruits, F=spp$t.series$F.full)
#write.csv(out, file="BAM.output.csv", quote=F)

