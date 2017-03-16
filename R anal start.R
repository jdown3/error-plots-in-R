#learn how to do things with nts results
ntcount <- read.table ("C:/Users/Julia/Documents/GitHub/error-profiles/results/output.nts.txt", header=T) #txt files in notepad, values separated by tabs 
#called it ntcount, \ is escape so need double or /

max <- signif(max(ntcount$Number, na.rm=T), 2) #signif rounds the first arg to the number of sigfig in the 2nd arg
#ceil <- ceiling(max(ntcount$Number, na.rm=T)) #want to round UP to 1 sig fig...

#install.packages("ggplot2")
library("ggplot2") #have to say at start to open package, like import in py

#PLOT 1 - number of each nt
ntplot <- ggplot(ntcount, aes(x = reorder(Nucleotide, -Number), y = Number + geom_bar(stat="identity") + geom_text(aes(label=Number)) + labs(title="Frequency of Nucleotides in Sequences", x="Nucleotide", y="Frequency") + ylim(0, max)
#plot from ntcount, nucleotide against number, make it a bar graph, order nucleotides by their number

#PLOT 1.5 - no of nts as %
#get individually, then put into 2nd col FIX
nts <- c("A", "C", "G", "T")
acount <- subset(ntcount, Nucleotide=="A", select=Number) #faster way than filtering?, only after one cell value
ccount <- subset(ntcount, Nucleotide=="C", select=Number) 
gcount <- subset(ntcount, Nucleotide=="G", select=Number)
tcount <- subset(ntcount, Nucleotide=="T", select=Number)
numnts <- sum(acount, ccount, gcount, tcount)
apct <- acount/numnts *100 
cpct <- ccount/numnts *100 
gpct <- gcount/numnts *100 
tpct <- tcount/numnts *100
pcts <- as.numeric(c(apct, cpct, gpct, tpct)) #number vector
ntpcts <- data.frame(nts, pcts, stringsAsFactors = FALSE)
ntpctP <-ggplot(ntpcts, aes(x = reorder(nts, -pcts), y = pcts)) + geom_bar(stat="identity", fill="#996699") + geom_text(aes(label=signif(pcts, 2), y=pcts + 1)) + labs(title="Nucleotide composition of Sequences", x="Nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(ntpcts$pcts)+2))

#errors
errors.df <- read.table ("C:/Users/Julia/Documents/GitHub/error-profiles/results/output.errors.txt", header=T)

emax <- signif(max(errors.df$Number, na.rm=T), 2)

errors.df$newcol <- paste(errors.df$Difference, errors.df$Consensus) #combines diff and cons columns to get change for x axis

#PLOT 2 - number of each error type 
eplot <- ggplot(errors.df, aes(x=reorder(newcol,-Number), y=Number)) + geom_bar(stat="identity", width=0.8, fill='#6666CC') + geom_text(aes(label=Number, y=Number+10000)) + labs(title="Count of error types", x="Change", y="Count") + theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, emax+20000))

#plot raw no. of each nt being cons and being diffs
#cons
aconslist <- subset(errors.df, Consensus=="A")
acons <- sum(aconslist$Number)
cconslist <- subset(errors.df, Consensus=="C")
ccons <- sum(cconslist$Number)
gconslist <- subset(errors.df, Consensus=="G")
gcons <- sum(gconslist$Number)
tconslist <- subset(errors.df, Consensus=="T")
tcons <- sum(tconslist$Number)

cons <- as.numeric(c(acons,ccons,gcons,tcons))
ntsAScons <- data.frame(nts,cons)
uplim <- max(ntsAScons$cons)
lowlim <- min(ntsAScons$cons)
#PLOT 3 - number of times each nt is cons
consp <- ggplot(ntsAScons, aes(x=reorder(nts,-cons), y=cons)) + geom_bar(stat="identity") + geom_text(aes(cons))
consp2 <- consp + labs(title="Frequency of each nucleotide as the Consensus", x="Nucleotide", y="Frequency") + coord_cartesian(ylim=c(lowlim,uplim))

#percentages 
ntsAScons$Pct <- c(ntsAScons$cons/sum(ntsAScons$cons)*100)
#PLOT 4 - % of each nt that will be cons
consprop <- ggplot(ntsAScons, aes(x=reorder(nts,-Pct), y=Pct)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(Pct, 2))) + labs(title="Probability of Consensus nucleotide", x="Nucleotide", y="Percentage (%)") 

#diff
adiffslist <- subset(errors.df, Difference=="A")
adiffs <- sum(adiffslist$Number) #DOESN'T USE DIFFNUM SO WRONG!
cdiffslist <- subset(errors.df, Difference=="C")
cdiffs <- sum(cdiffslist$Number)
gdiffslist <- subset(errors.df, Difference=="G")
gdiffs <- sum(gdiffslist$Number)
tdiffslist <- subset(errors.df, Difference=="T")
tdiffs <- sum(tdiffslist$Number)

diffs <- as.numeric(c(adiffs,cdiffs,gdiffs,tdiffs))
ntsASdiffs <- data.frame(nts,diffs)
uplim <- max(ntsASdiffs$diffs)
lowlim <- min(ntsASdiffs$diffs)
#PLOT 5 - number of times each nt is diff
diffsp <- ggplot(ntsASdiffs, aes(x=reorder(nts,-diffs), y=diffs)) + geom_bar(stat="identity") + geom_text(aes(label=diffs))
diffsp2 <- diffsp + labs(title="Frequency of each nucleotide as the Difference", x="Nucleotide", y="Frequency") + coord_cartesian(ylim=c(lowlim,uplim))

#percentages 
ntsASdiffs$Pct <- c(ntsASdiffs$diffs/sum(ntsASdiffs$diffs)*100)
#PLOT 6 - % of each nt that will be diffs
diffsprop <- ggplot(ntsASdiffs, aes(x=reorder(nts,-Pct), y=Pct)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(Pct,2))) + labs(title="Probability of Difference nucleotide", x="Nucleotide", y="Percentage (%)") #+ coord_cartesian(ylim=c(lowlim,uplim))

#raw no. of single letter being in error (cons or diff)
aerrors <- adiffs + acons
cerrors <- cdiffs + ccons
gerrors <- gdiffs + gcons
terrors <- tdiffs + tcons

errors <- as.numeric(c(aerrors, cerrors, gerrors, terrors)) #values with L at end, not nums from data frame
ntsASerrs <- data.frame(nts,errors) #FIX writing out four times
uplim <- max(ntsASerrs$errors)
lowlim <- min(ntsASerrs$errors)
#PLOT 7 - number of times each nt is in errors (diff or cons)
errsp <- ggplot(ntsASerrs, aes(x=reorder(nts,-errors), y=errors)) + geom_bar(stat="identity") + geom_text(aes(label=errors))
errsp2 <- errsp + labs(title="Frequency of each nucleotide in a change (diff or cons)", x="Nucleotide", y="Frequency") + coord_cartesian(ylim=c(lowlim,uplim))

#percentages
#make new % col in data frame
numerrors <- sum(errors.df$Number)
ntsASerrs$Pct <- c(aerrors/numerrors *100, cerrors/numerrors*100, gerrors/numerrors *100, terrors/numerrors*100)
#PLOT 8  -  % nt makeup of errors
errsprop <- ggplot(ntsASerrs, aes(x=reorder(nts,-Pct), y=Pct)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(Pct,2))) + labs(title="Probability of nucleotide occurance in a change (diff or cons)", x="Nucleotide", y="Percentage (%)") #+ coord_cartesian(ylim=c(lowlim,uplim))

consPctMakeup <- c(acons/numerrors *100, ccons/numerrors*100, gcons/numerrors *100, tcons/numerrors*100)
consmakeup <- data.frame(nts, consPctMakeup)
#PLOT 8c - % makeup of cons, num cons is num errors
cmakeupp <- ggplot(consmakeup, aes(x=reorder(nts,-consPctMakeup), y=consPctMakeup)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(consPctMakeup,2))) + labs(title="Probability of nucleotide in Consensus", x="Nucleotide", y="Percentage (%)")

diffsPctMakeup <- c(adiffs/numerrors *100, cdiffs/numerrors*100, gdiffs/numerrors *100, tdiffs/numerrors*100)
diffsmakeup <- data.frame(nts, diffsPctMakeup)
#PLOT 8d - % makeup of diffs
dmakeupp <- ggplot(diffsmakeup, aes(x=reorder(nts,-diffsPctMakeup), y=diffsPctMakeup)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(diffsPctMakeup,2))) + labs(title="Probability of nucleotide in Difference", x="Nucleotide", y="Percentage (%)")

#PLOT 8s - stacked % makeup of errors 
#make data frame
aecol <- c(adiffs/numerrors *100, acons/numerrors *100) #e for normalised w respect to errors
cecol <- c(cdiffs/numerrors *100, ccons/numerrors *100) #NUMBERS DON'T MAKE SENSE
gecol <- c(gdiffs/numerrors *100, gcons/numerrors *100)
tecol <- c(tdiffs/numerrors *100, tcons/numerrors *100)
dc <- c("Difference", "Consensus")
stackPctErrs <- data.frame( dc, aecol, cecol, gecol, tecol) #cols of data frames are variables

stackPctErrs2 <- plyr::rename(stackPctErrs, c("aecol"="A", "cecol"="C", "gecol"="G", "tecol"="T")) #so legend has right labels

stackPctErrslong <- melt(stackPctErrs2, id.var="dc") #"melt" into long table form for plotting

stackPctErrsP <- ggplot(stackPctErrslong, aes(x = dc, y = value, fill = variable, label=signif(value,2))) + geom_bar(stat = "identity") + geom_text(position=position_stack(vjust=0.5)) + labs(title="Probability of nucleotide in Errors", x="Error character", y="Probability (%)", fill="Nucleotide") + xlim(rev(levels(stackPctErrslong$dc)))+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0)) 
#how label top with total??

#normalised by nt count plot
#make data frame, div by..
ancol <- rbind(adiffs/acount*100, acons/acount*100) #normalised with respect to how many of itself in all seqs
cncol <- rbind(cdiffs/ccount*100, ccons/ccount*100)
gncol <- rbind(gdiffs/gcount*100, gcons/gcount*100)
tncol <- rbind(tdiffs/tcount*100, tcons/tcount*100)
stackPctErrsNN <- data.frame(dc, ancol, cncol, gncol, tncol) #NOT WORK SAME AS ABOVE
stackPctErrsNN <- plyr::rename(stackPctErrsNN, c("Number"="A", "Number.1"="C", "Number.2"="G", "Number.3"="T"))
stackPctErrsNNlong <- melt(stackPctErrsNN, id.var="dc") 
stackPctErrsNNP <- ggplot(stackPctErrsNNlong, aes(x = dc, y = value, fill = variable, label=signif(value,2))) + geom_bar(stat = "identity") + geom_text(position=position_stack(vjust=0.5)) + labs(title="Proportion of nucleotide in Errors", x="Error character", y="Percentage (%)", fill="Nucleotide") + xlim(rev(levels(stackPctErrsNNlong$dc)))+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0)) 

#normalised errors, calc then show
#in errors, add (+) numbers with g, divide (/) by g number in ntcount
#A
#acount is up top now #faster way than filtering?, only after one cell value
AinErrs <- aerrors/acount *100
AinDiff <- adiffs/acount *100
AinCons <- acons/acount *100

#C
CinErrs <- cerrors/ccount *100
CinDiff <- cdiffs/ccount *100
CinCons <- ccons/ccount *100

#G
#error number with g
#gerrorslist <- subset(errors, Difference=="G" | Concensus=="G", select=Number) #get me number where G is either 
#gerrors <- sum(gerrorslist$Number) #add into one int, SEPARATE? ADD DIFF AND CONS LISTS?
#must use == in subset fnc
GinErrs <- gerrors/gcount *100
GinDiff <- gdiffs/gcount *100
GinCons <- gcons/gcount *100

#T
TinErrs <- terrors/tcount *100
TinDiff <- tdiffs/tcount *100
TinCons <- tcons/tcount *100

#plot, plot fractions, 3 plots for errors, diffs, cons
#errors
#make data frame to plot
inErrs <- as.numeric(c(AinErrs, CinErrs, GinErrs, TinErrs)) #makes list, as.numeric changes it to number vector
ntinErrs <- data.frame(nts, inErrs)
#PLOT 9 - % of each nt that is an error
neplot <- ggplot(ntinErrs, aes(x=reorder(nts, -inErrs), y=inErrs)) + geom_bar(stat="identity") + geom_text(aes(label=signif(inErrs,2)))
neplot2 <- neplot + labs(title="Normalised single nucleotide changes", x="Nucleotide", y="Percentage of nucleotide present in Errors (%)")
uplim <- round(max(ntinErrs$inErrs), 2)
lowlim <- round(min(ntinErrs$inErrs), 2) 
neplot2 + coord_cartesian(ylim=c(lowlim,uplim)) #rounds up and down to 2dp, cause fractions will be between 1 and 10%?
#neplot2 + scale_y_continuous(limits=c(0.06,0.07)) #rids bars
#neplot2 + ylim(0.06,0.07) #rids bars

#diffs
inDiffs <- as.numeric(c(AinDiff, CinDiff, GinDiff, TinDiff)) #makes list, as.numeric changes it to number vector
ntinDiffs <- data.frame(nts, inDiffs)
#PLOT 10 - % of each nt that is a diff
ndplot <- ggplot(ntinDiffs, aes(x=reorder(nts, -inDiffs), y=inDiffs)) + geom_bar(stat="identity") + geom_text(aes(label=signif(inDiffs,2)))
ndplot2 <- ndplot + labs(title="Normalised single nucleotide changes (differences)", x="Nucleotide", y="Percentage of nucleotide present in the Difference (%)")
uplim <- round(max(ntinDiffs$inDiffs)+0.005, 2)
lowlim <- round(min(ntinDiffs$inDiffs)-0.005, 2)
ndplot2 + coord_cartesian(ylim=c(lowlim,uplim))

#cons
inCons <- as.numeric(c(AinCons, CinCons, GinCons, TinCons)) #makes list, as.numeric changes it to number vector
ntinCons <- data.frame(nts, inCons)
#PLOT 11 - % of each nt that is cons
ncplot <- ggplot(ntinCons, aes(x=reorder(nts, -inCons), y=inCons)) + geom_bar(stat="identity") + geom_text(aes(label=signif(inCons, 2)))
ncplot2 <- ncplot + labs(title="Normalised single nucleotide changes (concensus)", x="Nucleotide", y="Percentage of nucleotide present in the Concensus (%)")
uplim <- round(max(ntinCons$inCons)+0.005, 2)
lowlim <- round(min(ntinCons$inCons)-0.005, 2) 
ncplot2 + coord_cartesian(ylim=c(lowlim,uplim))

#flankings, long x axis, what's best plot/way to present?? BREAK UP
flanks_many <- read.table("C:/Users/Julia/Documents/GitHub/error-profiles/results/output.flankings_counts.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)
flanks <- read.table("C:/Users/Julia/Documents/GitHub/error-profiles/results/output.flankings_counts_single.txt", sep='\t', header=TRUE, stringsAsFactors=FALSE)
#remove rows with no seqs
flanks_many <- subset(flanks_many, Prev.Seq != " " | Difference != " " | Consensus != " " | Post.Seq != " ")# getting rid of end errors

fmax <- signif(max(flanks$Number, na.rm=T), 2)
  
other <- read.table("C:/Users/Julia/Documents/GitHub/error-profiles/results/output.other.txt", sep='\t', header=TRUE)

#PLOT 12 - number of each flanking seq **CRAMMED
fplot <- ggplot(flanks, aes(x=reorder(all,-Number), y=Number)) + geom_bar(stat="identity") + geom_text(aes(label=Number))
fplot + labs(title="Frequency of flanking Sequences at error sites", x="Flanking Sequence", y="Frequency") + ylim(0,fmax)

#normalised flankings
#% of each letter in flank, post seq, prev seq
#calcs
#A
aprevseq <- flanks[grep("A", flanks$Prev.Seq),] #all the rows were prev seq has an A, using regex
aprevlist <- subset(aprevseq, select=Number) #take the number column
aprev <- sum(aprevlist$Number)
apostseq <- flanks[grep("A", flanks$Post.Seq),]
apostlist <- subset(apostseq, select=Number)
apost <- sum(apostlist$Number)

aflanktable <- flanks[grep("A", flanks$newcol),]
aflanklist <- subset(aflanktable, select=Number)
aflank <- sum(aflanklist$Number)
#aflank <- aprev + adiffs + acons + apost #doubling up, need to find A in whole row, so use all cols together and then grep

AinFlank <- aflank/acount *100
AinPrev <- aprev/acount *100
AinPost <- apost/acount *100

#C
cprevseq <- flanks[grep("C", flanks$Prev.Seq),] #all the rows where prev seq has an A, using regex
cprevlist <- subset(cprevseq, select=Number) #take the number column
cprev <- sum(cprevlist$Number)
cpostseq <- flanks[grep("C", flanks$Post.Seq),]
cpostlist <- subset(cpostseq, select=Number)
cpost <- sum(cpostlist$Number)

cflanktable <- flanks[grep("C", flanks$newcol),]
cflanklist <- subset(cflanktable, select=Number)
cflank <- sum(cflanklist$Number)

CinFlank <- cflank/ccount *100
CinPrev <- cprev/ccount *100
CinPost <- cpost/ccount *100

#G
gprevseq <- flanks[grep("G", flanks$Prev.Seq),] #all the rows were prev seq has an A, using regex
gprevlist <- subset(gprevseq, select=Number) #take the number column
gprev <- sum(gprevlist$Number)
gpostseq <- flanks[grep("G", flanks$Post.Seq),]
gpostlist <- subset(gpostseq, select=Number)
gpost <- sum(gpostlist$Number)

gflanktable <- flanks[grep("G", flanks$newcol),]
gflanklist <- subset(gflanktable, select=Number)
gflank <- sum(gflanklist$Number)

GinFlank <- gflank/gcount *100
GinPrev <- gprev/gcount *100
GinPost <- gpost/gcount *100

#T
tprevseq <- flanks[grep("T", flanks$Prev.Seq),] #all the rows were prev seq has an A, using regex
tprevlist <- subset(tprevseq, select=Number) #take the number column
tprev <- sum(tprevlist$Number)
tpostseq <- flanks[grep("T", flanks$Post.Seq),]
tpostlist <- subset(tpostseq, select=Number)
tpost <- sum(tpostlist$Number)

tflanktable <- flanks[grep("T", flanks$newcol),]
tflanklist <- subset(tflanktable, select=Number)
tflank <- sum(tflanklist$Number)

TinFlank <- tflank/tcount *100 
TinPrev <- tprev/tcount *100
TinPost <- tpost/tcount *100

#plot
#flank 
inFlank <- as.numeric(c(AinFlank, CinFlank, GinFlank, TinFlank)) #makes list, as.numeric changes it to number vector
ntinFlank <- data.frame(nts, inFlank)
#PLOT 13 - % of each nt that is in flanking seq
nfplot <- ggplot(ntinFlank, aes(x=reorder(nts, -inFlank), y=inFlank)) + geom_bar(stat="identity") + geom_text(aes(label=signif(inFlank,2)))
nfplot2 <- nfplot + labs(title="Normalised flanking sequence changes", x="Nucleotide", y="Percentage of nucleotide present in Flanking sequence (%)")
uplim <- round(max(ntinFlank$inFlank) + 0.5, 0) #these numbers bigger, greater than 10%
lowlim <- round(min(ntinFlank$inFlank) - 0.5, 0)
nfplot2 + coord_cartesian(ylim=c(lowlim,uplim))

#prev
inPrev <- as.numeric(c(AinPrev, CinPrev, GinPrev, TinPrev)) #makes list, as.numeric changes it to number vector
ntinPrev <- data.frame(nts, inPrev)
uplim <- round(max(ntinPrev$inPrev) + 0.5, 0)
lowlim <- round(min(ntinPrev$inPrev) - 0.5, 0) 
#PLOT 14 - % of each nt that is in prev seq
nprplot <- ggplot(ntinPrev, aes(x=reorder(nts, -inPrev), y=inPrev)) + geom_bar(stat="identity", width=0.5, fill="orchid4")  + geom_text(aes(label=signif(inPrev,3), y=inPrev+1)) + labs(title="Normalised pre-error sequence", x="Nucleotide", y="Percentage of nucleotide present in pre-sequence (%)")  + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, uplim+2)) #gets rid of fill

#PLOT 14p - % makeup of prev seqs, any flankno
#need to know flanking no, get list of combos for x axis, count occurance of each combo
library(data.table)
flankno <- other$flankno
flanklengthcombos <- do.call(CJ, replicate(flankno, nts, FALSE)) #flankno=2, an int, get all 2 letter combos of nts
#put columns together with apply
combos <- apply(flanklengthcombos, 1, paste, collapse="") 
combos <- data.frame(cbind(combos)) #makes combo vector a column

#loop, for each row, subset, sum??
preMakeup <- data.frame(stringsAsFactors=FALSE)
#col added later is count in  all seqs
prevSeqs <- subset(flanks, select = c(Prev.Seq, Number)) #remove ones with no prev seq/two spaces, numbers were doubled up
#unique NOT WORKING
#prevSeqs <- unique(transform(prevSeqs, Number=ave(Number, Prev.Seq, FUN=sum)))
prevSeqs$Prev.Seq <- gsub('\\s+', '', prevSeqs$Prev.Seq) #DOES THIS AFFECT FUTURE THINGS? no
prevSeqs <- ddply(prevSeqs,.(Prev.Seq),summarize,Number=sum(Number))

for(i in 1:nrow(combos)) {
  row <- combos[i,] #factor w levels
  row <- as.character(row)
  pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 #many of them
  newrow <- cbind(row, pct)
  #output row and pct to dataframe
  preMakeup <- rbind(preMakeup, newrow) #binds new row to previous preMakeup table
}
#wasn't working because c(,) for newrow made a list, not a row

#DO NO MORE, in sep graph
#count combos of less letters separately, then add to one dataframe NOT FULLY AUTOMATED
#count combos of length flankno-1 until that =0?
if (flankno-1 != 0) {
  #make combos
  shorter1combos <- do.call(CJ, replicate(flankno-1, nts, FALSE))
  s1combos <- apply(shorter1combos, 1, paste, collapse="") 
  s1combos <- data.frame(cbind(s1combos)) 
  #do loop
  for(i in 1:nrow(s1combos)) {
    s1row <- s1combos[i,] 
    s1row <- as.character(s1row)
    s1pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
    s1newrow <- cbind(row, pct) #not getting
    #output row and pct to dataframe
    preMakeup <- rbind(preMakeup, s1newrow) 
  }
  
  if (flankno-2 != 0) {
    shorter2combos <- do.call(CJ, replicate(flankno-2, nts, FALSE))
    s2combos <- apply(shorter2combos, 1, paste, collapse="") 
    s2combos <- data.frame(cbind(s2combos)) 
    for(i in 1:nrow(s2combos)) {
      s2row <- s2combos[i,]
      s2row <- as.character(s2row)
      s2pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
      s2newrow <- cbind(row, pct)
      #output row and pct to dataframe
      preMakeup <- rbind(preMakeup, newrow) 
    }
    
    if (flankno-3 != 0) {
      shorter3combos <- do.call(CJ, replicate(flankno-3, nts, FALSE))
      s3combos <- apply(shortercombos, 1, paste, collapse="") 
      s3combos <- data.frame(cbind(s3combos)) 
      for(i in 1:nrow(s3combos)) {
        s3row <- s3combos[i,]
        s3row <- as.character(s3row)
        pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
        newrow <- cbind(row, pct)
        #output row and pct to dataframe
        preMakeup <- rbind(preMakeup, newrow) 
      }
      
      if (flankno-4 != 0) {
        shorter4combos <- do.call(CJ, replicate(flankno-4, nts, FALSE))
        s4combos <- apply(shorter4combos, 1, paste, collapse="") 
        s4combos <- data.frame(cbind(s4combos)) 
        for(i in 1:nrow(s4combos)) {
          s4row <- s4combos[i,]
          s4row <- as.character(s4row)
          pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
          newrow <- cbind(row, pct)
          #output row and pct to dataframe
          preMakeup <- rbind(preMakeup, newrow) 
        }
      }
      
    }
    
  }
  
}
#14p
preMakeupP <- ggplot(preMakeup, aes(x=reorder(row,-Number), y=Number)) + geom_bar(stat="identity", width=0.8, fill="#3333CC") + geom_text(aes(label=signif(Number, 2), y=Number+0.3)) + labs(title="Probability of pre-error sequence", x="Sequence", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(preMakeup$Number)+1))

#WHEN FLANKNO = 1 #CHANGE TO APRE, CPRE...., DON'T NEED??
#prePctMakeup <- c(acons/numerrors *100, ccons/numerrors*100, gcons/numerrors *100, tcons/numerrors*100)
#premakeup <- data.frame(nts, prePctMakeup)
#premakeupp <- ggplot(premakeup, aes(x=reorder(nts,-prePctMakeup), y=prePctMakeup)) + geom_bar(stat="identity") +  geom_text(aes(label=signif(prePctMakeup,2))) + labs(title="Probability of pre-error sequence", x="Sequence", y="Percentage (%)")

#post
inPost <- as.numeric(c(AinPost, CinPost, GinPost, TinPost)) #makes list, as.numeric changes it to number vector
ntinPost <- data.frame(nts, inPost)
uplim <- round(max(ntinPost$inPost)+0.5, 0)
lowlim <- round(min(ntinPost$inPost)-0.5, 0)
#PLOT 15 - % of each nt that is in post seq
npoplot <- ggplot(ntinPost, aes(x=reorder(nts, -inPost), y=inPost)) + geom_bar(stat="identity", width=0.5, fill="#336699") + geom_text(aes(label=signif(inPost,3), y=inPost + 1)) + labs(title="Normalised post-error sequence", x="Nucleotide", y="Percentage of nucleotide present in post-sequence (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits=c(0, uplim+2))

#PLOT 15p - % makeup of post seqs
postMakeup <- data.frame(stringsAsFactors=FALSE)
postSeqs <- subset(flanks, select = c(Post.Seq, Number)) #remove ones with no prev seq/two spaces, numbers were doubled up
#postSeqs <- unique(transform(postSeqs, Number=ave(Number, Post.Seq, FUN=sum))) #ERROR NOW
#cbind(aggregate(Number~Post.Seq, sum, data=postSeqs), table(postSeqs$Post.Seq)) #ERROR
#dtpostSeqs <- data.table(postSeqs)
#postSeqs <- dtpostSeqs[,list(Number = sum(Number), freq = .N), by = "Post.Seq"] #double ups
postSeqs$Post.Seq <- gsub('\\s+', '', postSeqs$Post.Seq) 
postSeqs <- ddply(postSeqs,.(Post.Seq),summarize,Number=sum(Number)) #remove white space b4
for(i in 1:nrow(combos)) {
  row <- combos[i,]
  row <- as.character(row)
  pct <- subset(postSeqs, Post.Seq==row, select=Number)/numerrors *100 #not getting
  newrow <- cbind(row, pct)
  postMakeup <- rbind(postMakeup, newrow) #binds new row to previous preMakeup table
}
#DO NO MORE
if (flankno-1 != 0) {
  #don't have to make combos, made above
  for(i in 1:nrow(s1combos)) {
    s1row <- s1combos[i,] 
    s1row <- as.character(s1row)
    s1pct <- subset(postSeqs, Post.Seq==row, select=Number)/numerrors *100 
    s1newrow <- cbind(row, pct) #not getting
    #output row and pct to dataframe
    postMakeup <- rbind(postMakeup, s1newrow) 
  }
  
  if (flankno-2 != 0) {
    for(i in 1:nrow(s2combos)) {
      s2row <- s2combos[i,]
      s2row <- as.character(s2row)
      s2pct <- subset(postSeqs, Post.Seq==row, select=Number)/numerrors *100 
      s2newrow <- cbind(row, pct)
      #output row and pct to dataframe
      postMakeup <- rbind(postMakeup, newrow) 
    }
    
    if (flankno-3 != 0) {
      for(i in 1:nrow(s3combos)) {
        s3row <- s3combos[i,]
        s3row <- as.character(s3row)
        pct <- subset(postSeqs, Post.Seq==row, select=Number)/numerrors *100 
        newrow <- cbind(row, pct)
        #output row and pct to dataframe
        postMakeup <- rbind(postMakeup, newrow) 
      }
      
      if (flankno-4 != 0) {
        for(i in 1:nrow(s4combos)) {
          s4row <- s4combos[i,]
          s4row <- as.character(s4row)
          pct <- subset(postSeqs, Post.Seq==row, select=Number)/numerrors *100 
          newrow <- cbind(row, pct)
          #output row and pct to dataframe
          postMakeup <- rbind(postMakeup, newrow) 
        }
      }
      
    }
    
  }
  
}
#15p
postMakeupP <- ggplot(postMakeup, aes(x=reorder(row,-Number), y=Number)) + geom_bar(stat="identity", width=0.8, fill="#663399") + geom_text(aes(label=signif(Number, 2), y=Number+0.3)) + labs(title="Probability of post-error sequence", x="Sequence", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(postMakeup$Number)+1))

#stack 14 and 15
stackntinpp <- data.frame(nts, ntinPrev$inPrev, ntinPost$inPost)
stackntinpp <- plyr::rename(stackntinpp, c("ntinPrev.inPrev"="Pre-error seq", "ntinPost.inPost"="Post-error seq"))
stackntinpplong <- melt(stackntinpp, id.var="nts")

order <- with(stackntinpplong, order(nts, value)) #doesn't work?

#PLOT 14-15 - for each nt as diffs, probability of what the cons will be
stackntinPPplot <- ggplot(stackntinpplong, aes(x = reorder(nts, -value), y = value, fill = variable, label=signif(value,2))) + 
  geom_bar(stat = "identity") + geom_text(position=position_stack(vjust=0.5)) + labs(title="Proportion of nucleotides surrounding errors?", x="Nucleotide", y="Percentage (%)", fill="Fragment position")+ theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0))
#move over

#stacked plot - diffs
#if A is diffs, prob of others being cons, get from errors data file
adiffccons <- sum(subset(adiffslist, Consensus=="C", select=Number))
adiffgcons <- sum(subset(adiffslist, Consensus=="G", select=Number))
adifftcons <- sum(subset(adiffslist, Consensus=='T', select=Number))
#get %
adiffcconsPct <- adiffccons/adiffs *100
adiffgconsPct <- adiffgcons/adiffs * 100
adifftconsPct <- adifftcons/adiffs *100

#cdiffslist, cdiffs is the sum 
cdiffacons <- sum(subset(cdiffslist, Consensus=="A", select=Number))
cdiffgcons <- sum(subset(cdiffslist, Consensus=="G", select=Number))
cdifftcons <- sum(subset(cdiffslist, Consensus=='T', select=Number))
#get %
cdiffaconsPct <- cdiffacons/cdiffs *100
cdiffgconsPct <- cdiffgcons/cdiffs * 100
cdifftconsPct <- cdifftcons/cdiffs *100

#repeat for other nts
gdiffacons <- sum(subset(gdiffslist, Consensus=="A", select=Number))
gdiffccons <- sum(subset(gdiffslist, Consensus=="C", select=Number))
gdifftcons <- sum(subset(gdiffslist, Consensus=='T', select=Number))
#get %
gdiffaconsPct <- gdiffacons/gdiffs *100
gdiffcconsPct <- gdiffccons/gdiffs * 100
gdifftconsPct <- gdifftcons/gdiffs *100

tdiffacons <- sum(subset(tdiffslist, Consensus=="A", select=Number))
tdiffccons <- sum(subset(tdiffslist, Consensus=="C", select=Number))
tdiffgcons <- sum(subset(tdiffslist, Consensus=='G', select=Number))
#get %
tdiffaconsPct <- tdiffacons/tdiffs *100
tdiffcconsPct <- tdiffccons/tdiffs * 100
tdiffgconsPct <- tdiffgcons/tdiffs *100

#got values, plot
#make data frame
aconscol <- c(0, cdiffaconsPct, gdiffaconsPct, tdiffaconsPct) #these will be diff colours
cconscol <- c(adiffcconsPct, 0, gdiffcconsPct, tdiffcconsPct)
gconscol <- c(adiffgconsPct, cdiffgconsPct, 0, tdiffgconsPct)
tconscol <- c(adifftconsPct, cdifftconsPct, gdifftconsPct, 0)

stackdiffs <- data.frame(nts, aconscol, cconscol, gconscol, tconscol) #cols of data frames are variables
library(plyr)
stackdiffs2 <- rename(stackdiffs, c("aconscol"="A", "cconscol"="C", "gconscol"="G", "tconscol"="T")) #so legend has right labels
library(reshape2)
stackdiffslong <- melt(stackdiffs2, id.var="nts") #"melt" into long table form for plotting
#remove 0s
stackdiffslong2 <- subset(stackdiffslong, value != 0)

#order variable by value, reordered dataframe
dorder <- with(stackdiffslong2, order(nts, value))

#PLOT 16 - for each nt as diffs, probability of what the cons will be
stackdiffsplot <- ggplot(stackdiffslong2[dorder,], aes(x = nts, y = value, fill = variable, label=signif(value,2))) + 
  geom_bar(stat = "identity") + geom_text(position=position_stack(vjust=0.5)) + labs(title="Probability of the Consensus nucleotide per Difference nucleotide", x="Difference nucleotide", y="Probability (%)", fill="Consensus nucleotide")
#position stack gets value labels to be in the stacks, vjust 0.5 puts it halfway up the stack

#stacked plot - cons
aconscdiff <- sum(subset(aconslist, Difference=="C", select=Number))
aconsgdiff <- sum(subset(aconslist, Difference=="G", select=Number))
aconstdiff <- sum(subset(aconslist, Difference=='T', select=Number))
#get %
aconscdiffPct <- aconscdiff/acons *100
aconsgdiffPct <- aconsgdiff/acons *100
aconstdiffPct <- aconstdiff/acons *100

cconsadiff <- sum(subset(cconslist, Difference=="A", select=Number))
cconsgdiff <- sum(subset(cconslist, Difference=="G", select=Number))
cconstdiff <- sum(subset(cconslist, Difference=='T', select=Number))
#get %
cconsadiffPct <- cconsadiff/ccons *100
cconsgdiffPct <- cconsgdiff/ccons *100
cconstdiffPct <- cconstdiff/ccons *100

gconsadiff <- sum(subset(gconslist, Difference=="A", select=Number))
gconscdiff <- sum(subset(gconslist, Difference=="C", select=Number))
gconstdiff <- sum(subset(gconslist, Difference=='T', select=Number))
#get %
gconsadiffPct <- gconsadiff/gcons *100
gconscdiffPct <- gconscdiff/gcons *100
gconstdiffPct <- gconstdiff/gcons *100

tconsadiff <- sum(subset(tconslist, Difference=="A", select=Number))
tconscdiff <- sum(subset(tconslist, Difference=="C", select=Number))
tconsgdiff <- sum(subset(tconslist, Difference=='G', select=Number))
#get %
tconsadiffPct <- tconsadiff/tcons *100
tconscdiffPct <- tconscdiff/tcons * 100
tconsgdiffPct <- tconsgdiff/tcons *100

#plot - cons
adiffscol <- c(0, cconsadiffPct, gconsadiffPct, tconsadiffPct) #these will be diff colours
cdiffscol <- c(aconscdiffPct, 0, gconscdiffPct, tconscdiffPct)
gdiffscol <- c(aconsgdiffPct, cconsgdiffPct, 0, tconsgdiffPct)
tdiffscol <- c(aconstdiffPct, cconstdiffPct, gconstdiffPct, 0)

stackcons <- data.frame(nts, adiffscol, cdiffscol, gdiffscol, tdiffscol) #cols of data frames are variables
#library(plyr)
stackcons2 <- rename(stackcons, c("adiffscol"="A", "cdiffscol"="C", "gdiffscol"="G", "tdiffscol"="T")) #so legend has right labels
#library(reshape2)
stackconslong <- melt(stackcons2, id.var="nts") #"melt" into long table form for plotting
stackconslong2 <- subset(stackconslong, value != 0)

#order variable by value, reordered dataframe
corder <- with(stackconslong2, order(nts, value))

#PLOT 17 - for each nt as cons, probability of what the diff will be
stackconsplot <- ggplot(stackconslong2[corder,], aes(x = nts, y = value, fill = variable, label=signif(value,2))) + 
  geom_bar(stat = "identity") + geom_text(position=position_stack(vjust=0.5)) + labs(title="Probability of the Difference nucleotide per Consensus nucleotide", x="Consensus nucleotide", y="Probability (%)", fill="Difference nucleotide")

#count of prepost together
bflanklengthcombos <- do.call(CJ, replicate(flankno*2, nts, FALSE))
bcombos <- apply(bflanklengthcombos, 1, paste, collapse="") 
bcombos <- data.frame(cbind(bcombos))

flanks$pp <- paste(flanks$Prev.Seq, flanks$Post.Seq, sep="") #double ups!!!!!!! does matter?
#flanks$pp <- gsub('\\s+', '', flanks$pp) #rids whitespace from pp col CAN'T DO
#ppflanks <- subset(flanks, Prev.Seq != " "*flankno & Difference == " " & Consensus == " " & Post.Seq != " "*flankno)

both <- data.frame()
for(i in 1:nrow(bcombos)) {
  row <- bcombos[i,]
  pct <- subset(flanks, pp==row, select=Number)/numerrors *100 #this can be multiple rows, need sum
  pct <- sum(pct$Number)
  newrow <- cbind(row, pct)
  both <- rbind(both, newrow) #binds new row to previous preMakeup table
}
#JUST DO ABOVE FOR NOW.......not using in markdown...
#for both of less letters
if (flankno-1 != 0) { 
  #make combos
  bshort1combos <- do.call(CJ, replicate(flankno*2-1, nts, FALSE))
  bs1combos <- apply(bshort1combos, 1, paste, collapse="") 
  bs1combos <- data.frame(cbind(bs1combos)) 
  #do loop
  for(i in 1:nrow(bs1combos)) {
    bs1row <- bs1combos[i,] 
    bs1pct <- subset(ppflanks, pp==row, select=Number)/numerrors *100 
    bs1newrow <- cbind(row, pct) #not getting
    #output row and pct to dataframe
    preMakeup <- rbind(preMakeup, s1newrow) 
  }
  
  if (flankno-2 != 0) {
    shorter2combos <- do.call(CJ, replicate(flankno-2, nts, FALSE))
    s2combos <- apply(shorter2combos, 1, paste, collapse="") 
    s2combos <- data.frame(cbind(s2combos)) 
    for(i in 1:nrow(s2combos)) {
      s2row <- s2combos[i,]
      s2pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
      s2newrow <- cbind(row, pct)
      #output row and pct to dataframe
      preMakeup <- rbind(preMakeup, newrow) 
    }
    
    if (flankno-3 != 0) {
      shorter3combos <- do.call(CJ, replicate(flankno-3, nts, FALSE))
      s3combos <- apply(shortercombos, 1, paste, collapse="") 
      s3combos <- data.frame(cbind(s3combos)) 
      for(i in 1:nrow(s3combos)) {
        row <- s3combos[i,]
        pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
        newrow <- cbind(row, pct)
        #output row and pct to dataframe
        preMakeup <- rbind(preMakeup, newrow) 
      }
      
      if (flankno-4 != 0) {
        shorter4combos <- do.call(CJ, replicate(flankno-4, nts, FALSE))
        s4combos <- apply(shorter4combos, 1, paste, collapse="") 
        s4combos <- data.frame(cbind(s4combos)) 
        for(i in 1:nrow(s4combos)) {
          row <- s4combos[i,]
          pct <- subset(prevSeqs, Prev.Seq==row, select=Number)/numerrors *100 
          newrow <- cbind(row, pct)
          #output row and pct to dataframe
          preMakeup <- rbind(preMakeup, newrow) 
          
          5
          
          6
          
          7
          
          8
          
          9
        }
      }
      
    }
    
  }
  
}
#PLOT 23 - number of prepost seqs (both), percent makeup
bothp <- ggplot(both, aes(x=reorder(row, -Number), y=Number)) + geom_bar(stat="identity", width=0.8, fill="#336699") + geom_text(aes(label=round(Number,2), y=Number+0.1)) + labs(title="Probability of prepost error sequences", x= "Pre and Post error sequence", y="Probability (%)")+ theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(both$Number)+0.1))

#IN WHOLE/TOTAL FLANK SEQ

#install bioconduct, shortread
fastq <- "C:/Users/Julia/Documents/GitHub/error-profiles/data/Jurkat_only_S5.consensus.fastq"
seqs <- ShortRead::readFastq("C:/Users/Julia/Documents/GitHub/error-profiles/data/Jurkat_only_S5.consensus.fastq")
#have to flip slashes
#shortreadQ files have id, sread and quality slots
seqs <- as.character(seqs) #no can do

#every 4th line...
#ShortRead::FastqStreamer(seqs, 4) #chunks, ERROR in function, unable to find inherited method 
  
#PLOT 24 - normalised pre

#for each row in prevseqs..count in sread...
col <- data.frame()
for(i in 1:nrow(prevSeqs)) { #does prevSeqs have double ups? NO using flanks_many
  row <- prevSeqs[i,]
  print(row)
  #overlapping?
  count <- sum(stringr::_count(_count(ShortRead::sread(seqs), row$Prev.Seq)) #will give you count for each seq
  #add col?
  print(count)
  col <- rbind(col, count)
}
#add to df that gets plotted 
prevSeqs <- cbind(prevSeqs, col) #overwrite or add another col? -add

prevSeqs <- plyr::rename(prevSeqs, c("X2715934L"="col")) #same num each run?

#divide number by overall count
npreplot <- ggplot(prevSeqs, aes(x=reorder(Prev.Seq, -(Number/col)), y=Number/col*100)) + geom_bar(stat="identity", width=0.8, fill="#660033")+ geom_text(aes(label=round(Number/col*100,1), y=(Number/col*100)+0.1)) + labs(title="Proportion of sequence fragments that are pre-error", x= "Pre-error sequence", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(prevSeqs$Number/prevSeqs$col*100)+0.1))
  
#PLOT 25 - normalised post
postSeqs$col <- 0
for(i in 1:nrow(postSeqs)) {
  row <- postSeqs[i,]
  print(row)
  #overlapping?
  count <- sum(stringr::str_count(ShortRead::sread(seqs), row$Post.Seq)) #will give you count for each seq
  print(count)
  postSeqs$col[i] <- count
}

postSeqs <- subset(postSeqs, col != 0) #QUICK FIX FOR NOW, RID LESS LETTERS PROPERLY LATER

#divide number by overall count
npostplot <- ggplot(postSeqs, aes(x=reorder(Post.Seq, -(Number/col)), y=Number/col*100)) + geom_bar(stat="identity", width=0.8, fill="#660066")+ geom_text(aes(label=round(Number/col*100,1), y=(Number/col*100)+0.1)) + labs(title="Proportion of sequence fragments that are post-error", x= "Post-error sequence", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) +scale_y_continuous(expand=c(0,0), limits = c(0, max(postSeqs$Number/postSeqs$col*100)+0.5))
  
#PLOT 26 - normalised pdp
flanks$pdp <- paste(flanks$Prev.Seq, flanks$Difference, flanks$Post.Seq, sep='')
pdp <- as.character(flanks$pdp)

#errors with different consensus letters double up, should be less than 16384 (all poss 7 letter combos)
#subset flanks to get pdps, where no cons, cons =" " and diff = something
#REDO with this data, how fast 2day method now? ******* back to counter with enviro

pdp <- subset(flanks, select = c(Prev.Seq, Difference, Post.Seq, Number)) #removes cons
pdp$pdpcol <- paste(pdp$Prev.Seq, pdp$Difference, pdp$Post.Seq, sep='')
#remove white spaces
pdp$pdpcol <- gsub('\\s+', '', pdp$pdpcol) #replacing one or more spaces with no space
pdp$Prev.Seq <- NULL
pdp$Difference <- NULL
pdp$Post.Seq <- NULL
pdp <- unique(transform(pdp, Number=ave(Number, pdpcol, FUN=sum))) #unique pdpcols summed??? 22611
#remove ones of less letters
pdp$nchar <- nchar(pdp$pdpcol)
pdp <- subset(pdp, nchar == 7, select=c(pdpcol, Number)) #remove frags of less here at end or earlier?
#remove pdpcols with Ns
#pdp[grep("N", pdp$pdpcol),] #how remove not keep, none

#faster counting all pdps in seq at once
lines <- readLines(fastq) #remove seqs?

#THIS METHOD 2 DAYS, fixed not regex
flanks$count <- 0
for (i in 1:length(lines)){ #for each seq 
  time <- Sys.time()
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    line = lines[i]
    print(line)
    count <- stringr::str_count(line, stringr::fixed(pdp)) #default is regex, I only need it to look for fixed string e.g. ACA #vectorise only one of these, otherwise will match up instead of looking all
    print(count) #count is vector, sum counts
    flanks$count <- flanks$count + count #sum in loop, only one col?
  }
  print(Sys.time()-time)
}
#flanks_single$Total <- rowSums(flanks_single[ , 7:608862]) #TOO SLOW

#COUNTER OF COMBOS METHOD, might be slow cause 7 letter combos
#all flank seq length k-mers, here 7
p_pcombos <- do.call(CJ, replicate(flankno*2 + 1, nts, FALSE))
p_pcombos <- apply(p_pcombos, 1, paste, collapse="") #char array/vector
p_pcombos <- data.frame(combos=p_pcombos, stringsAsFactors = FALSE) #keys/index for counter

#TRY USING ENVIRONMENT TO MAKE DICT COUNTER
counters <- new.env(parent=emptyenv(), size=nrow(p_pcombos)*1.2) #reset enviro, removes counters
#assign all the letter combos, if not found they'll be zero #make combos into vector/array?

#loop through seqs and add to them 3HRS
for (i in 1:length(lines)){ #for each seq 
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    line = lines[i]
    len <- nchar(line)
    begin = 1
    endd = flankno*2+1 #reset for each seq
    #print(line)
    #print(len)
    #sliding window count, flankno*2+1 length
    #get subseq
    
    #substr
    #get substr to slide along by incrementing the start and stop,
    while (endd <= len){ #200ish loops #for faster than while, if only once
      frag <- substr(line, begin, endd)
      #print(frag)
      #then increment counter for whatever substr it gets
      counters[[frag]] <- mget(frag, counters, ifnotfound = 0)[[1]] + 1
      #print(counters[[frag]]) 
      begin = begin + 1
      endd = endd + 1
      #print(begin)
    } 
    #have to loop through combos?? NO MORE
    #for (i in 1:length(p_pcombos)){ #16000 loops
      #p_p <- p_pcombos[i] #current combo
      #count <- stringr::str_count(line, stringr::fixed(p_p)) #search #vectorise?? or substr for sliding window
      #print(count)
      #e[[p_p]] <- counters + count #add #ERROR
    #}
    if (i %% 1000 == 2){ #every 250 seqs, has to get to 2mill, change 4000
      print(i)
      print(Sys.time()-time) #to process a single seq
    }
  }
}
#unlist counters output to df to plot later...
p_pcombos$count <- unlist(mget(p_pcombos$combos, envir = counters, ifnotfound = 0))
#plot df

#just plot count of combos RUN
p_pcombosP <- ggplot(p_pcombos, aes(x=reorder(combos, -count), y=count)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(count,0), y=(p_pcombos$count)+0.1)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(p_pcombos$count)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#not a pct? of all 7-mers
#break up or zoom in with plotly?

total7mers <- sum(p_pcombos$count)

#web-based?
#devtools::install_github("ropensci/plotly") #4.5.6
library(plotly)
library(ggplot2) #2.2.1.900
ggplotly(p_pcombosP) #use CRAN version of plotly 
#try older plotly (dev), dev ggplot2?
#CRAN 4.5.2
devtools::install_github("ropensci/plotly", ref = "85128db")
#error #cannot open file 'C:/Users/Julia/Documents/R/win-library/3.3/plotly/R/sysdata.rdb': No such file or directory

#3-mers
combos3 <- do.call(CJ, replicate(3, nts, FALSE))
combos3 <- apply(combos3, 1, paste, collapse="") #char array/vector
combos3 <- data.frame(combos=combos3, stringsAsFactors = FALSE) #keys/index for counter

count3mers <- new.env(parent=emptyenv(), size=nrow(combos3)*1.2) 

#loop through seqs 1HR!!
for (i in 1:length(lines)){ #for each seq 
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    line = lines[i]
    len <- nchar(line)
    begin = 1
    endd = 3 #reset for each seq
    while (endd <= len){ #200ish loops #for faster than while, if only once
      frag <- substr(line, begin, endd)
      count3mers[[frag]] <- mget(frag, count3mers, ifnotfound = 0)[[1]] + 1
      #print(counters[[frag]]) 
      begin = begin + 1
      endd = endd + 1
      #print(begin)
    } 
    if (i %% 400000 == 2){ #every 100,000 seqs
      print(i)
      print(Sys.time()-time) #to process a single seq
    }
  }
}
combos3$count <- unlist(mget(combos3$combos, envir = count3mers, ifnotfound = 0)) #NOT FOUND
combos3P <- ggplot(combos3, aes(x=reorder(combos, -count), y=count)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(count,0), y=(count)+0.1)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(count)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

total3mers <- sum(combos3$count)

combos3$pct <- combos3$count/total3mers *100 #proportion of all 3mers, NOT error rate

#3mers without first one
count3mers2 <- new.env(parent=emptyenv(), size=nrow(combos3)*1.2) 

#loop through seqs 1HR!!
for (i in 1:length(lines)){ #for each seq 
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    line = lines[i]
    len <- nchar(line)
    begin = 2 #this is the change! so missing first 3mer
    endd = 4 #reset for each seq
    while (endd <= len){ #200ish loops #for faster than while, if only once
      frag <- substr(line, begin, endd)
      count3mers2[[frag]] <- mget(frag, count3mers2, ifnotfound = 0)[[1]] + 1
      #print(counters[[frag]]) 
      begin = begin + 1
      endd = endd + 1
      #print(begin)
    } 
    if (i %% 400000 == 2){ #every 100,000 seqs
      print(i)
      print(Sys.time()-time) #to process a single seq
    }
  }
}
combos3$count2 <- unlist(mget(combos3$combos, envir = count3mers2, ifnotfound = 0)) 
combos3P <- ggplot(combos3, aes(x=reorder(combos, -count2), y=count2)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(combos3$count2,0), y=(combos3$count2)+0.1)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(combos3$count2)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

total3mers2 <- sum(combos3$count2)

combos3$pct2 <- combos3$count2/total3mers2 *100


#loop through seqs and count them/inc counter
#set up counters to 0
#counter=defaultdict(int) #not in R

ctr=0 #counter is just an int var that you add to, will have to loop through pdpcombos

#NOT MUCH FASTER? cause 2 loops and str_count
#make counter a matrix?
counter=matrix(NA, nrow=length(p_pcombos), ncol=2) #worked??
rownames(counter)=p_pcombos
colnames = c("combos", "count")
colnames(counter)=colnames #didnt' work in one line

for (i in 1:length(lines)){ #for each seq 
  time <- Sys.time()
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    line = lines[i]
    print(line)
    #have to loop through combos?
    for (i in 1:length(p_pcombos)){
      p_p <- p_pcombos[i] #current combo
      count <- stringr::str_count(line, stringr::fixed(p_p))
      print(count)
      counter[p_p,2]=counter[p_p,2]+count #increment by count of times in that seq
    }
  }
  print(Sys.time()-time)
}

#INSTEAD, 2nd col of combos df is counter?
p_pcombos$pdpcounter <- 0

for (i in 1:length(lines)){ #for each seq 
  time <- Sys.time()
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    line = lines[i]
    print(line)
    
    #increment counter, sliding window? rollapply?
    counter[p_pcombos] + 1
    
    count <- stringr::str_count(line, stringr::fixed(pdp)) #default is regex, I only need it to look for fixed string e.g. ACA #vectorise only one of these, otherwise will match up instead of looking all
    print(count) #count is vector, sum counts
    flanks$count <- flanks$count + count #sum in loop, only one col?
  }
  print(Sys.time()-time)
}

###OTHER ATTEMPTS
##strsplit(readLines(seqs))

##fastq-grep('AAA', fastq)

##HTSeqGenie::FastQStreamer.getReads(fastq)

fastqchunk <- ShortRead::FastqStreamer(fastq, 1) #chunk is 4 reads or 4 lines?

while (length(fq <- ShortRead::yield(f))){
  print(read)
  #do stuff, count like b4
}
close(fastqchunk)

###strm <- ShortRead::FastqStreamer("C:/Users/Julia/Documents/GitHub/error-profiles/data/Jurkat_only_S5.consensus.fastq")
repeat {
  print(strm)
}


###for (i in 1:nrow(ShortRead::sread(seqs))){ #ARG OF LENGTH 0 ERROR
  seq <- ShortRead::sread(seqs)[i,]

  for(i in 1:nrow(flanks)) {
    row <- flanks[i,]
    print(row)
    count <- sum(stringr::str_count(seq, row$pdp)) #will give you count for each seq
    print(count)
    col <- rbind(col, count)
  }
}

#PLOT put number cols into same df
pdpP <- merge(pdp, p_pcombos, by.x = "pdpcol", by.y = "combos") #merge, these columns have the common row names 
pdpP$pct <- pdpP$Number/pdpP$count *100
pdpP <- pdpP[rev(order(pdpP$pct)),]
#divide Number by overall count
npdpplot <- ggplot(pdpP, aes(x=reorder(pdpcol, -(Number/count)), y=Number/count*100)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(Number/count*100,1), y=(Number/count*100)+0.1)) + labs(title="Proportion of sequence fragments that are the Difference flanking sequence", x= "Fragment (pre-diff-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max((pdpP$Number/pdpP$count)*100)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#x axis labels vertical
#can't see, break up.....

#plot 5000, third, 250 1/60
#order first
pdpPless <- pdpP[1:200,]
npdpPfront <- ggplot(pdpPless, aes(x=reorder(pdpcol, -(Number/count)), y=Number/count*100)) + geom_bar(stat="identity", width=0.8, fill="#CC99CC")+ geom_text(aes(label=round(Number/count*100,1), y=(Number/count*100)+8, angle=90)) + labs(title="Proportion of sequence fragments that are the Difference flanking sequence", x= "Fragment (pre-diff-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max((pdpPless$Number/pdpPless$count)*100)+15)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot last 500
pdpPend <- pdpP[16184:16384,]
npdpPback <- ggplot(pdpPend, aes(x=reorder(pdpcol, -(Number/count)), y=Number/count*100)) + geom_bar(stat="identity", width=0.8, fill="#CCCCFF")+ geom_text(aes(label=round(Number/count*100,2), y=(Number/count*100)+0.025, angle=90)) + labs(title="Proportion of sequence fragments that are the Difference flanking sequence", x= "Fragment (pre-diff-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max((pdpPend$Number/pdpPend$count)*100)+0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#PLOT 27 - normalised pcp
#get PCP, remove things
pcp <- subset(flanks, select = c(Prev.Seq, Difference, Consensus, Post.Seq, Number)) #removes cons
pcp$pcpcol <- paste(pcp$Prev.Seq, pcp$Consensus, pcp$Post.Seq, sep='')
pcp$Prev.Seq <- NULL
pcp$Consensus <- NULL
pcp$Post.Seq <- NULL
#remove white spaces
pcp$pcpcol <- gsub('\\s+', '', pcp$pcpcol) #replacing one or more spaces with no space
#remove ones of less letters
pcp$nchar <- nchar(pcp$pcpcol)
pcp <- subset(pcp, nchar == flankno*2+1, select=c(pcpcol, Number, Difference))
#pcp <- unique(transform(pcp, Number=ave(Number, pcpcol, FUN=sum))) #unique pdpcols summed??? 20750???

#DO PCP FOR 3MERS, F=1, have deno, need numerator/error number
pcp3 <- subset(flanks, select = c(Prev.Seq, Consensus, Post.Seq, Number))
#get prev and post to 1 letter
#loop through rows, few/4 mins?
for (i in 1 :nrow(pcp3)){
  pcp3$Prev.Seq[i] <- substr(pcp3$Prev.Seq[i], nchar(pcp3$Prev.Seq[i]), nchar(pcp3$Prev.Seq[i]))
  pcp3$Post.Seq[i] <- substr(pcp3$Post.Seq[i], 1, 1)
}
#paste, unique
pcp3$pcpcol <- paste(pcp3$Prev.Seq, pcp3$Consensus, pcp3$Post.Seq, sep='')
pcp3$Prev.Seq <- NULL
pcp3$Consensus <- NULL
pcp3$Post.Seq <- NULL
pcp3 <- ddply(pcp3,.(pcpcol),summarize,Number=sum(Number))

#take out 2mers
pcp3no2 <- pcp3 #need to keep difference col for new pcp plot
pcp3no2$pcpcol <- gsub('\\s+', '', pcp3no2$pcpcol) #replacing one or more spaces with no space
pcp3no2$nchar <- nchar(pcp3no2$pcpcol)
pcp3no2 <- subset(pcp3no2, nchar == 3, select=c(pcpcol, Number))

#calc pct, merge
pcp3P <- merge(pcp3no2, combos3, by.x = "pcpcol", by.y = "combos") #merge, these columns have the common row names 
pcp3P$pct <- NULL
pcp3P$count2 <- NULL #remove unneccessary cols that came from combos3
pcp3P$pct2 <- NULL
pcp3P$pct <- pcp3P$Number/pcp3P$count *100 #error rate
pcp3P <- pcp3P[rev(order(pcp3P$pct)),]

#plot
pcp3plot <- ggplot(pcp3P, aes(x=reorder(pcpcol, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(pct,1), y=(pct)+0.1)) + labs(title="Proportion of sequence fragments that are the Consensus flanking sequence", x= "Fragment (pre-cons-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pct)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymax = pct, ymin=pct)) # + geom_errorbar(data=pcpP, aes(ymax = pcpP$pct + pcpP$se, ymin= pcpP$pct - pcpP$se), stat = "identity", position = "dodge", width=0.25)


#NEW PCP DF FOR STACKED PLOT
#set up lists
pcp <- pcp[order(pcp$pcpcol),] #either alpha or rev
npcpcol <- vector(mode='list', nrow(upcp))
anum <- vector(mode='list', nrow(upcp))
cnum <- vector(mode='list', nrow(upcp))
gnum <- vector(mode='list', nrow(upcp))
tnum <- vector(mode='list', nrow(upcp))
pcpcol <- "BBBBBBB"
n <- 0
for (i in 1:nrow(pcp)){ #50 000 rows
  num <- pcp$Number[i] #number for that diff
  diff <- pcp$Difference[i]
  #printf("diff %s", diff)
  if (pcpcol == pcp$pcpcol[i]){
    #just add diff nums
    if (diff=="A"){
      anum[[n]] <- num
    }
    if (diff=="C"){
      cnum[[n]] <- num
    }
    if (diff=="G"){
      gnum[[n]] <- num
    }
    if (diff=="T"){
      tnum[[n]] <- num
    }
  } else { #if null checking here what happens for first row
    #new row 
    #increment n
    n <- n + 1
    #printf("n %f", n)
    #new pcp
    pcpcol <- pcp$pcpcol[i]
    #printf("pcp %s", pcpcol)
    npcpcol[[n]] <- pcpcol
    #add nums
    if (diff=="A"){
      anum[[n]] <- num
    }
    if (diff=="C"){
      cnum[[n]] <- num
    }
    if (diff=="G"){
      gnum[[n]] <- num
    }
    if (diff=="T"){
      tnum[[n]] <- num
    }
    #make prev row empty cells 0
    if (is.null(anum[[n]])){ #not the best spot in loop, goes through every row
      anum[[n]] <- 0
    }
    if (is.null(cnum[[n]])){
      cnum[[n]] <- 0
    }
    if (is.null(gnum[[n]])){
      gnum[[n]] <- 0
    }
    if (is.null(tnum[[n]])){
      tnum[[n]] <- 0
    }
  }
  #printf("anum %f", anum[[n]])
  #printf("gnum %f", gnum[[n]])
  #printf("cnum %f", cnum[[n]])
  #printf("tnum %f", tnum[[n]])
  printf("for this row n: %f pcp: %s diff: %s num: %f anum: %f cnum: %f gnum: %f tnum: %f", n, pcpcol, diff, num, anum[[n]], cnum[[n]], gnum[[n]], tnum[[n]])
}
#make df or add lists to df as cols
Npcp <- cbind(npcpcol, anum, cnum, gnum, tnum) #get matrix
Npcp <- as.data.frame(Npcp)
Npcp$anum <- as.numeric(Npcp$anum) #to unlist cell entries, str npcpcol?
Npcp$cnum <- as.numeric(Npcp$cnum)
Npcp$gnum <- as.numeric(Npcp$gnum)
Npcp$tnum <- as.numeric(Npcp$tnum)
Npcp$npcpcol <- as.character(Npcp$npcpcol)
total <- vector(mode='list', nrow(Npcp))
for (i in 1:nrow(Npcp)){
  total[[i]] <- sum(Npcp$anum[i], Npcp$cnum[i], Npcp$gnum[i], Npcp$tnum[i])
}
Npcp$total <- total
Npcp$total <- as.numeric(Npcp$total)

#order df by total then remove?
Npcp <- Npcp[rev(order(Npcp$total)),]#descending, still intact when long for plots?

Npcp$total <- NULL

#plot this df instead for pcp!!! w error bars
Npcp <- plyr::rename(Npcp, c("anum"="A", "cnum"="C", "gnum"="G", "tnum"="T")) #so legend has right labels

Npcplong <- melt(Npcp, id.var="npcpcol") #"melt" into long table form for plotting
#error - Can't melt data.frames with non-atomic 'measure' columns
#just as.numeric the columns

NpcpP <- ggplot(Npcplong, aes(x = reorder(npcpcol, -value), y = value, fill=variable, colour=variable)) + geom_bar(stat = "identity") + labs(title="Frequency of F-mers, broken up by the difference", x="F-mer", y="Frequency", fill="Difference") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0)) + scale_colour_discrete(guide = FALSE)

#merge df
pcpP <- merge(upcp, p_pcombos, by.x = "pcpcol", by.y = "combos") #merge, these columns have the common row names 
pcpP$pct <- pcpP$Number/pcpP$count *100
pcpP <- pcpP[rev(order(pcpP$pct)),]
#divide Number by overall count
npcpplot <- ggplot(pcpP, aes(x=reorder(pcpcol, -pct), y=Number/count*100)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(Number/count*100,1), y=(Number/count*100)+0.1)) + labs(title="Proportion of sequence fragments that are the Consensus flanking sequence", x= "Fragment (pre-cons-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max((pcpP$Number/pcpP$count)*100)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymax = pct + se, ymin=pct - se)) # + geom_errorbar(data=pcpP, aes(ymax = pcpP$pct + pcpP$se, ymin= pcpP$pct - pcpP$se), stat = "identity", position = "dodge", width=0.25)
#change all num/count to pct

#ERROR Error in pmin(y, 0) : object 'y' not found
#add standard error bars
#se needs to be column in df? = sd/ sqrt(n), add col to pcpP, n will be count, sd?
pcpP$se <- (sqrt(((pcpP$pct/100)*(1-(pcpP$pct/100)))/pcpP$count))*100

#compare package didn't work
#use dplyr 
#what type of join?
#join p_pcombos and pcpP, it will find combos and give it it's pcp "count"
#left join, p_pcombos x/first arg, just leaves stuff from x
#anti_join, no match
#full_join, all cols from both, NA when no match
p_pApcp <- dplyr::full_join(p_pcombos, upcp, by = c("combos"="pcpcol")) #RUN WITH UPCP
p_pApcp[is.na(p_pApcp)] <- 0 #Number col from pcp, so error freq 
#join
#then make pct col
p_pApcp$pct <- p_pApcp$Number/p_pApcp$count *100

#replace pcpP? fix  whole PCP plot so zeros show at end!!!!!

p_pcombosSort <- ggplot(p_pApcp, aes(x=reorder(combos, -pct), y=count)) + geom_bar(stat="identity", width=0.8, fill="#FF6666", colour="#FF6666")+ geom_text(data=subset(p_pApcp, count>=50000), aes(label=combos, y=count*1.05, angle=90)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(p_pApcp$count)*1.2)) + theme(axis.text.x = element_blank())
#sort by pcp order/"Number"

#first 200
pcpPless <- pcpP[1:200,] #this could be changed, no biggie
npcpPfront <- ggplot(pcpPless, aes(x=reorder(pcpcol, -(pct)), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#CC99CC") + labs(title="Proportion of sequence fragments that are the Consensus flanking sequence", x= "Fragment (pre-cons-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pcpPless$pct)+ max(pcpPless$se)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymax = pcpPless$pct + se, ymin= pcpPless$pct - se), position="dodge", width=0.5) + geom_text(aes(label=round(pcpPless$pct,1), y=pcpPless$pct + se + 0.25, angle=90))

#debug pcp double ups
pcpd <- subset(flanks, select = c(Prev.Seq, Difference, Consensus, Post.Seq, Number)) #removes cons
pcpd$pcpcol <- paste(pcpd$Prev.Seq, pcpd$Consensus, pcpd$Post.Seq, sep='')

#sort by pct
p_pApcp <- p_pApcp[with(p_pApcp, order(-pct, combos)),]
p_pApcpfront <- p_pApcp[1:200,] #comma then empty means all cols
#plot freq of those 7mers
p_pcombosSortCfront <- ggplot(p_pApcpfront, aes(x=reorder(combos, -pct), y=count)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(count,0), y=(p_pApcpfront$count)+0.1, angle=90)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(p_pApcpfront$count)*1.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot last 200
pcpPend <- p_pApcp[16185:16384,] 
npcpPback <- ggplot(pcpPend, aes(x=reorder(pcpcol, -(Number/count)), y=Number/count*100)) + geom_bar(stat="identity", width=0.8, fill="#CCCCFF")+ geom_text(aes(label=round(Number/count*100,2), y=(Number/count*100)+0.025, angle=90)) + labs(title="Proportion of sequence fragments that are the Difference flanking sequence", x= "Fragment (pre-diff-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max((pcpPend$Number/pcpPend$count)*100)+0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot freq of those 7mers
p_pcombosSortCback <- ggplot(pcpPend, aes(x=reorder(combos, -pct), y=count)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(count,0), y=(pcpPend$count)*1.05, angle=90)) + labs(title="Frequency of f-mers in sample", x= "F-mer", y="Frequency")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pcpPend$count)*1.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#4 FRAGS OF LESS plots 
firstfew <- prevSeqs[grep(" ", prevSeqs$Prev.Seq),] #rows where prev seq has a space
pos2 <- subset(prevSeqs, grepl("  [[:upper:]]", Prev.Seq)) #rid white space?
pos3 <- prevSeqs[grep(" [[:upper:]][[:upper:]]", prevSeqs$Prev.Seq),] #full stop for any letter, 'metacharacter', need to use regex/grep
#no subset for 3 spaces cause no prev 

#multiple subsets, 2 in this case

#pos2 plot
#after subset, sum number col of each subset to get num errors at each pos, normalising by
sumP2 <- sum(pos2$Number)

#put number from each row of subset on top of total number for subset
pos2$pct <- pos2$Number/sumP2*100

#plot Prev.Seq against pct from subset pos2
pos2P <- ggplot(pos2, aes(x=reorder(Prev.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#990099")+ geom_text(aes(label=round(pos2$pct,0), y=pos2$pct+2)) + labs(title="2nd nt errors", x= "Previous nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pos2$pct)+5))

#pos3 plot
sumP3 <- sum(pos3$Number)
pos3$pct <- pos3$Number/sumP3*100
pos3P <- ggplot(pos3, aes(x=reorder(Prev.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#990099")+ geom_text(aes(label=round(pos3$pct,1), y=pos3$pct+1)) + labs(title="3rd nt errors", x= "Previous nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pos3$pct)+2))

#replicate, no subset for all spaces cause no post seq
lastfew <- postSeqs[grep(" ", postSeqs$Post.Seq),]
#pos=len no post
#len-1
len1 <- postSeqs[grep("[[:upper:]]  ", postSeqs$Post.Seq),] #not working?? UP TO
#len-2
len2 <- postSeqs[grep("[[:upper:]][[:upper:]] ", postSeqs$Post.Seq),]

#len-1 plot
sumL1 <- sum(len1$Number)
len1$pct <- len1$Number/sumL1*100
len1P <- ggplot(len1, aes(x=reorder(Post.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#660066")+ geom_text(aes(label=round(len1$pct,1), y=len1$pct+2)) + labs(title="2nd last nt errors", x= "Following nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(len1$pct)+4))

#len-2 plot
sumL2 <- sum(len2$Number)
len2$pct <- len2$Number/sumL2*100
len2P <- ggplot(len2, aes(x=reorder(Post.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#660066")+ geom_text(aes(label=round(len2$pct,1), y=len2$pct+0.5)) + labs(title="3rd last nt errors", x= "Following nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(len2$pct)+2))

#AUTO FOR FLANKNO

#LOOP flankno-1 until that equals 1, everything inc plot in loop?

#for first few, pos=flankno is innermost 
#rid white space from prevSeqs
#prevSeqs$Prev.Seq <- gsub('\\s+', '', prevSeqs$Prev.Seq)
i  <- 1
P <- list() #loop will make a list of plots
while (flankno-i > 0){
  pos <- subset(prevSeqs, nchar(Prev.Seq)==(flankno-i)) #rename?
  sum <- sum(pos$Number)
  pos$pct <- pos$Number/sum*100
  P[[i]] <- ggplot(pos, aes(x=reorder(Prev.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#990099")+ geom_text(aes(label=round(pct,2), y=pct+2)) + labs(title="First few errors", x= "Previous nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pos$pct)+3)) 
  i = i + 1
}
#creates 2 plots; pos[1], pos[2]

#for last few, fix above first
#postSeqs$Post.Seq <- gsub('\\s+', '', postSeqs$Post.Seq) #done earlier 
i  <- 1
G <- list() #loop will make a list of plots
while (flankno-i > 0){
  pos <- subset(postSeqs, nchar(Post.Seq)==(flankno-i)) #rename?
  sum <- sum(pos$Number)
  pos$pct <- pos$Number/sum*100
  G[[i]] <- ggplot(pos, aes(x=reorder(Post.Seq, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#990099")+ geom_text(aes(label=round(pct,2), y=pct+2)) + labs(title="First few errors", x= "Following nucleotide", y="Percentage (%)") + theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pos$pct)+3)) 
  i = i + 1
}


printf <- function(...) invisible(print(sprintf(...)))

info <- matrix(NA, nrow=3555963, ncol=6) #3 mill errors, storing all might be slow

n <- 1 
cstart <- 0
gend <- 0

#half hour loop
for (i in 1:length(lines)){ #for each line of fastq
  #print(i)
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    seq = lines[i]
    seqno <- (i/4 + 0.5)
    printf("seqno %f", seqno) #i for integer
    len <- nchar(line)
    #if C first
    cstartB <- 0 #means no
    gendB <- 0 
    if (substring(seq,1,1) == 'C'){ #works?? 
      cstart <- cstart + 1
      cstartB <- 1
    }
    #if G last
    if (substring(seq,len,1)=='G'){
      gend <- gend + 1
      gendB <- 1
    }
  }
  if (i %% 4 == 3){
    diffs = lines[i]
    #print(diffs)
    diffs = strsplit(diffs, " ") #creates a list of vector elements 
    print(diffs)
    #error no.
    errsperseq <- length(diffs[[1]])
    printf("errs %e", errsperseq) 
    #iterate through diffs, to get posnos and store data
    for (i in 1:errsperseq){
      diff = diffs[[1]][i]
      pos <- as.numeric(gsub("([0-9]+).*$", "\\1", diff)) #THIS gets me first number, understand
      printf("posno %g", pos) 
      info[n,] <- c(seqno,errsperseq,diff,pos, cstartB, gendB) #print i (error no) or errsperseq
      n <- n + 1
    }
  }
  if (i %% 4000 == 2){ #every 1000seqs
    print(i)
    print(Sys.time()-time) #to process a single seq
  }
}
#don't need to print
#matrix to df
info <- data.frame(info) 
info <- plyr::rename(info, c("X1"="Seq no.", "X2"="Errors in seq", "X3"="Diff", "X4"="Pos no.", "X5"="C start", "X6"="G end"))


#INFO2 1 row per seq, vectors in cells
info2 <- matrix(NA, nrow=608856, ncol=13) #3 mill errors, storing all might be slow

n <- 1 
cstart <- 0
astart <- 0
gstart <- 0
tstart <- 0
gend <- 0
aend <- 0
cend <- 0
tend <- 0

#set size of pos and poslist
poslist <- vector(mode='list', 608856) #ave errs per seq 5.84 X 608856 or just 608856
#posstartList <- vector(mode='list', 608856)
#posrightList <- vector(mode='list', 608856)
#posendList <- vector(mode='list', 608856)
#posleftList <- vector(mode='list', 608856)

for (i in 1:length(lines)){ #for each line of fastq
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    seq = lines[i]
    seqno <- (i/4 + 0.5)
    printf("seqno %f", seqno) #i for integer
    len <- nchar(seq)
    #if C first
    cstartB <- 0 #means no
    astartB <- 0
    gstartB <- 0
    tstartB <- 0
    gendB <- 0 
    aendB <- 0
    cendB <- 0
    tendB <- 0
    if (substring(seq,1,1) == 'C'){ #works?? 
      cstart <- cstart + 1
      cstartB <- 1
    }
    if (substring(seq,1,1) == 'A'){ 
      astart <- astart + 1
      astartB <- 1
    }
    if (substring(seq,1,1) == 'G'){ 
      gstart <- gstart + 1
      gstartB <- 1
    }
    if (substring(seq,1,1) == 'T'){  
      gstart <- gstart + 1
      gstartB <- 1
    }
    #if G last
    if (substring(seq,len,len)=='G'){ #substring args are data,start,stop
      gend <- gend + 1
      gendB <- 1
    }
    if (substring(seq,len,len)=='A'){ #slow?? with all these loops NAH 15MINS
      aend <- aend + 1
      aendB <- 1
    }
    if (substring(seq,len,len)=='C'){ 
      cgend <- cend + 1
      cendB <- 1
    }
    if (substring(seq,len,len)=='T'){ 
      tend <- tend + 1
      tendB <- 1
    }
  }
  if (i %% 4 == 3){
    diffs = lines[i]
    #remove plus start
    diffs <- substring(diffs, 2)
    diffs = strsplit(diffs, " ") #creates a list of vector elements 
    print(diffs)
    #error no.
    errsperseq <- length(diffs[[1]])
    printf("errs %e", errsperseq) 
    #iterate through diffs, to get posnos
    m <- 1 #reset
    #c <- 1
    #g <- 1
    pos <- logical(errsperseq) 
    #posstart <- logical(errsperseq)
    #posrightmaj <- logical(errsperseq)
    #posend <- logical(errsperseq)
    #posleftmaj <- logical(errsperseq)
    #print(pos)
    if (errsperseq > 0){
      for (i in 1:errsperseq){ #TWO DIFF Is WORKS???? yes, cause used just under their own loops
        diff = diffs[[1]][i]
        pos[m] <- as.numeric(gsub("([0-9]+).*$", "\\1", diff))#THIS gets me first number, understand
        #if (cstartB == 1){
          #if (pos[m] <= (0.1*len)){
            #posstart[c] <- pos[m]
          #} else {
            #posrightmaj[c] <- pos[m]
          #}
          #c <- c+1
        #} else {
          #posstart <- NULL
          #posrightmaj <- NULL
        #}
        #if (gendB == 1){
          #if (pos[m]>= (0.9*len)){
            #posend[g] <- pos[m]
          #} else {
            #posleftmaj[g] <- pos[m]
          #}
          #g <- g + 1
        #} else {
          #posend <- NULL
          #posleftmaj <- NULL
        #}
        #print(pos)
        m <- m + 1
      }
      print(pos)
      poslist[[n]] <- pos #add row/element, double bracket to access single element?
      #remove 0s or falses??
      #posstart <- posstart[posstart != "0"] #takes out zeros????
      #posstartList[[n]] <- posstart
      #posrightmaj <- posrightmaj[posrightmaj != "0"]
      #posrightList[[n]] <- posrightmaj
      #posend <- posend[posend != "0"]
      #posendList[[n]] <- posend
      #posleftmaj <- posleftmaj[posleftmaj != "0"]
      #posleftList[[n]] <- posleftmaj
    }
    info2[n,] <- c(seqno, len, errsperseq, diffs, cstartB, gstartB, tstartB, astartB, gendB, cendB, tendB, aendB) #order of error occurance
    n <- n + 1  
  }
  if (i %% 4000 == 2){ #every 1000seqs
    print(i)
    print(Sys.time()-time) #to process a single seq
  }
}
#don't need to print
#matrix to df
info2 <- data.frame(info2) 
info2<- plyr::rename(info2, c("X1"="SeqNo.", "X2"="Length", "X3"="Errors", #ADD DIFFS and num reads!! "X4"="Cstart", "X5"="Gstart", "X6"="Tstart", "X7"="Astart", "X8"="Gend","X9"="Cend","X10"="Tend","X11"="Aend"))
#make pos list then add list as column
info2$pos <- poslist

info2$posStart <- posstartList
info2$posRight <- posrightList
info2$posEnd <- posendList
info2$posLeft <- posleftList 


#if C first, % of error after
#subset flanks - sum, subset of info2, with c start- nrows
#AUTO FLANKNO, spaces before
CstartER2 <- flanks[grep("  C", flanks$Prev.Seq),] #C then anything, full stop gets spaces, flanks should still have spaces, or extract first letter 
numCstartER2 <-sum(CstartER2$Number)

Cprob <- (numCstartER2/cstart)*100 #ANSWER

#start vs throughout? pos, count how many before and after a boundary number, first 10%??
#info2, count len of each row of column? iterate through rows or just add len column?

#CAN I SPLIT POS OUTSIDE LOOP, need to iterate through pos elements?
#subset info2 to get all cstart seqs
Cstart <- subset(info2, Cstart==1)
Cstart$posStart <- Cstart$pos[Cstart$pos <= (0.1*Cstart$Length)] #split 10/90 of len
Cstart$posRightMaj <- Cstart$pos[Cstart$pos > (0.1*Cstart$Length)] 
#(list) object cannot be coerced to type 'double'
#works for single vector, iterate thru rows *********** then run pcp/call, then C end, hists, other qs

CstartEr <- subset(info2, Cstart==1 & Errors > 0)

Cstartrow <- CstartEr[1:2303,] #extract rows 

#iterate, backup
#NEW SHORTER VERSION
posStartList <- vector(mode='list', nrow(CstartEr))
posRightList <- vector(mode='list', nrow(CstartEr))
m <- 1
for (i in 1:nrow(CstartEr)){ #ABORTING!, from adding to list, from brackets
  print(i)
  poss <- CstartEr$pos[i][[1]] #all pos from that seq row
  posStart <- poss[poss <= (0.1*CstartEr$Length[i])]
  posRightMaj <- poss[poss > (0.1*CstartEr$Length[i])]
  if (length(posStart)>0){ #avoid? what happens? null/empty?
    posStartList[[i]] <- posStart
  } else {
    posStartList[[i]] <- "none" #slower 
  }
  if (length(posRightMaj)>0) {
    posRightList[[i]] <- posRightMaj #to get to a list element you have to double bracket!!!
  } else {
    posRightList[[i]] <- "none"
  }
}
CstartEr$posStart <- posStartList
CstartEr$posRightMajority <- posRightList

#COMPARE NT AT START TO NUM ERRORS!!!!!!!!

#CstartEr$posStart <- CstartEr$pos[[1]][CstartEr$pos[[1]] <= (0.1*CstartEr$Length)] #doesn't get me whole pos column, so still need to iterate through rows but not pos vector

#ERRORS START VS THROUGHOUT

CstartposStart <- subset(CstartEr, posStart != "none")
numposStart <- vector(mode='list', nrow(CstartposStart))
for (i in 1:nrow(CstartposStart)){ 
  print(i)
  poss <- CstartposStart$pos[i][[1]] #all pos from that seq row
  numposStart[i] <- length(poss)
  numposStar <- length(poss)
  numCstartposStart <- numCstartposStart + numposStar 
}
CstartposStart$numposStart <- numposStart
#error rate, need total number of nts in fron 10% of all cstart seqs...

#don't need separate df, CstartEr has the posStart column, REMOVE OTHER LOOPS, RUNNNNNN
#192,881 seqs have Cstart and zero errors
#3 lists - overall seq error rate, start error rate, rightmaj error rate
seqErRates <- vector(mode='list', nrow(CstartEr))
startErRates <- vector(mode='list', nrow(CstartEr))
rightErRates <- vector(mode='list', nrow(CstartEr))
for (i in 1:nrow(CstartEr)){ #have to do [[1]] for posStart too
  print(i)
  seqErRate <- length(CstartEr$pos[i][[1]])/CstartEr$Length[i] *100
  startErRate <- 0 #reset for each seq
  rightErRate <- 0
  if (CstartEr$posStart[i][[1]] != "none") {
    startErRate <- (length(CstartEr$posStart[i][[1]])/(CstartEr$Length[i]*0.1))*100
  }
  if (CstartEr$posRightMajority[i][[1]] != "none") {
    rightErRate <- (length(CstartEr$posRightMajority[i][[1]])/(CstartEr$Length[i]*0.9))*100
  }
  seqErRates[[i]] <- seqErRate
  startErRates[[i]] <- startErRate
  rightErRates[[i]] <- rightErRate
}
#add cols to df
CstartEr$ErrorRate <- seqErRates
CstartEr$StartErrorRate <- startErRates
CstartEr$RightErrorRate <- rightErRates
#determine average?

#do right
CstartposRight <- subset(CstartEr, posRightMajority != "none")
numposRight <- vector(mode='list', nrow(CstartposRight))
numCstartposRight <- 0
for (i in 1:nrow(CstartposRight)){ 
  print(i)
  poss <- CstartposRight$pos[i][[1]] #all pos from that seq row
  numposRight[i] <- length(poss)
  numposRi <- length(poss)
  numCstartposRight <- numCstartposRight + numposRi #ANSWERRR?
}
CstartposRight$numposRight <- numposRight

#if G last 
GendERb4 <- flanks[grep("G  ", flanks$Post.Seq),] 
numGendERb4 <-sum(GendERb4$Number)

Gprob <- (numGendERb4/gend)*100 #ANSWER 

#end vs throughout? pos, count how many before and after a boundary number, need error rate, to compare to other nts 
COPY ABOVE

#BOTH/OVERLAPS make new df where both = 1


#C end, triplets
Cend <- subset(info2, Cend==1)

CendEr <- subset(info2, Cend==1 & Errors>0)
Cendtriplist <- logical(nrow(CendEr))
for (i in 1:nrow(CendEr)){
  seq=lines[CendEr$SeqNo.[i]*4-2]
  len <- nchar(seq) #NEW, TRY AGAIN
  Cendtrip=substr(seq,len-2,len)
  Cendtriplist[i] <- Cendtrip
}
CendEr$CendTriplet <- Cendtriplist

#get count, w unique or table or plot...
#make output var, do for start and all nts
CendErTrips <- sort(table(CendEr$CendTriplet)/nrow(CendEr)*100) #percents
#plot

#CendTrips

#CendER2trips
  
CendER2 <- subset(flanks, Post.Seq=="C..")

#C start, triplets
#Cstart dfs made above
Cstarttriplist <- logical(nrow(Cstart)) #CHANGED FROM CstartEr
for (i in 1:nrow(Cstart)){
  seq=lines[Cstart$SeqNo.[i]*4-2]
  len <- nchar(seq) 
  Cstarttrip=substr(seq,1,3)
  Cstarttriplist[i] <- Cstarttrip
}
Cstart$CstartTriplet <- Cstarttriplist

#get count, w unique or table or plot...
#make output var, do for start and all nts
CstartTrips <- sort(table(Cstart$CstartTriplet)/nrow(Cstart)*100) #percents
#plot


#how many pcp frags above or below 3.2%
#subset pcp$pct, pcp RUN
belowEER <- subset(pcpP, pct<3.20189657) #good
belowEER <- nrow(belowEER) #ANSWER

aboveEER <- subset(pcpP, pct>=3.20189657) #sep at exact error rate
aboveEER <- nrow(aboveEER)


#how many seqs have errors at pos 2 and 3
#subset info
#WRITE
Er2and3 <- subset(info2, pos==2 & 3)# ***
#2 or 3 total?

Er2 <- info2[grep("2", info2$pos),] #works??
Er3 <- info2[grep("3", info2$pos),]
#put together and unique by seq no?
Er2and3 <- info2[grep("2", info2$pos),] & info2[grep("3", info2$pos),]

Er2and3 <- info2[grep("2" & "3", info2$pos[[1]]),] #just gets me row1?

#pos column is lists so iterate through rows...
info2er <- subset(info2, Errors>0) #RUN

#Er2and3list <- vector(mode='list', nrow(info2)*0.25) #quarter of how many seqs
numEr2and3 <- 0
for (i in 1:nrow(info2er)){ 
  print(i)
  poss <- info2er$pos[i][[1]] #all pos from that seq row
  match <- match(c(2,3),poss, nomatch=0) #tells you pos of matches
  #remove 0s
  match <- match[match != "0"]
  len <- length(match)
  #iterate through pos or pull out from vector?
  if (len == 2){ #=true? what's output of grep WORKS? with spaces surrounding
    numEr2and3 <- numEr2and3 + 1 #ANSWERRRRR 
  }
}
#info2er$Er2and3 <- Er2and3list


#per base error rates loop, considering all reads 45MINS
#set things up
nts <- c("A", "C", "G", "T") 
#vector of nts, could be useful
seqnolist <- vector(mode='list', 608856)
numreadslist <- vector(mode='list', 608856)
errsperseqlist <- vector(mode='list', 608856)
Anolist <- vector(mode='list', 608856)
Cnolist <- vector(mode='list', 608856)
Gnolist <- vector(mode='list', 608856)
Tnolist <- vector(mode='list', 608856)
s <- 0
Ano <- 0
Cno <- 0
Gno <- 0
Tno <- 0
ACcount <- 0 #replacements in R for acons, ccons etc, number of times it is the cons char in errors
CCcount <- 0
GCcount <- 0
TCcount <- 0
#now a total number/var, doesn't have to be a column/per seq
Adcount <- 0 #when A is the cons, how many diffs in that error
Cdcount <- 0
Gdcount <- 0
Tdcount <- 0
for (i in 1:length(lines)){ #for each line of fastq
  #print(i)
  #header line has num reads
  if (i %% 4 == 1) { 
    header <- lines[i]
    #numreads <- substring(header, 51, 53) #strip colons later
    numreads <- strsplit(header, split=":")[[1]] #need num form cause going to multiply later
    length <- length(numreads)
    numreads <- as.numeric(numreads[length-3]) #get 4th last thing in header
  }
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    s <- s + 1
    #time <- Sys.time()
    seq = lines[i]
    seqno <- (i/4 + 0.5)
    printf("seqno %f", seqno) #i for integer
    #len <- nchar(seq)
    #seq specific nt counts, replaced each loop 
    Ano <- str_count(seq, "A")
    Cno <- str_count(seq, "C")
    Tno <- str_count(seq, "T")
    Gno <- str_count(seq, "G")
  }
  if (i %% 4 == 3){
    diffs = lines[i]
    #remove plus start
    diffs <- substring(diffs, 2)
    diffs = strsplit(diffs, " ")[[1]] #creates a list of vector elements 
    #print(diffs)
    #error no.
    errsperseq <- length(diffs)
    #printf("errs %e", errsperseq) 
    if (errsperseq > 0){
      for (k in 1:errsperseq){ #going through errors FASTER WAY???
        diff = diffs[k]
        #pull out nums and letters
        posno <- as.numeric(strsplit(diff, split="[[:upper:]]")[[1]][1]) #first num, strsplit makes list
        conschar <- substr(seq, posno, posno) #start and stop
        #consnum, could take it out and be left with diffnum
        cons <- str_extract(diff, paste0(conschar, "[[:digit:]]+"))
        if (conschar == "A"){
          ACcount <- ACcount + 1
        }
        if (conschar == "C"){
          CCcount <- CCcount + 1
        }
        if (conschar == "G"){
          GCcount <- GCcount + 1
        }
        if (conschar == "T"){
          TCcount <- TCcount + 1
        }
        #consnum <- str_extract(cons, "[[:digit:]]+") #number after cons letter NEED???
        #take out cons from error
        diffwocons <- gsub(cons, "", diff)
        diffnum <-strsplit(diffwocons, split="[[:upper:]]")[[1]][-1] #takes out first num/posno, length 2 OR MORE/5
        diffchar <-strsplit(diffwocons, split="[[:digit:]]+")[[1]][-1] #plus means any size digit, first one empty for some reason?
        #could be 4 diff chars in one err!
        #so need loop through diffnums
        for (j in 1:length(diffnum)){ #goes through numbers of each error 
          #diffnum not in quotes?
          adiffnum <- as.numeric(diffnum[j]) #the current diff num
          #add the diff letter to the right Ndcount
          if (conschar == "A"){ #or loop through nts ???
            Adcount <- Adcount + adiffnum #usually 1 but sometimes more
          }
          if (conschar == "C"){
            Cdcount <- Cdcount + adiffnum
          }
          if (conschar == "G"){
            Gdcount <- Gdcount + adiffnum
          }
          if (conschar == "T"){
            Tdcount <- Tdcount + adiffnum
          }
          #if (conschar == "N"){ #cons will never be N
            #Ndcount <- Ndcount + adiffnum
          #}
        }
        #keep the one that is not consnum
        #if (diffnum[1]!=consnum){ #both characters/in quotes
          #diffnum <- diffnum[1]
        #} else if (diffnum[2]!= consnum) { 
          #diffnum <- diffnum[2]
        #}
      }
    }
    #for each seq
    seqnolist[[s]] <- seqno
    numreadslist[[s]] <- numreads
    errsperseqlist[[s]] <- errsperseq
    Anolist[[s]] <- Ano
    Cnolist[[s]] <- Cno
    Gnolist[[s]] <- Gno
    Tnolist[[s]] <- Tno
  }
  #if (i %% 4000 == 2){ #every 1000seqs
    #print(i)
    #print(Sys.time()-time) #to process a single seq
  #}
}
#list to columns/df
pberr <- data.frame(cbind(seqnolist, numreadslist, errsperseqlist, Anolist, Cnolist, Gnolist, Tnolist))
#rename?

#unlist, can as.num the whole df? just means I can use operators, not unlisting
pberr$seqnolist <- as.numeric(pberr$seqnolist)
pberr$errsperseqlist <- as.numeric(pberr$errsperseqlist)

pberr$numreadslist <- as.numeric(pberr$numreadslist)
pberr$Anolist <- as.numeric(pberr$Anolist)
pberr$Cnolist <- as.numeric(pberr$Cnolist)
pberr$Gnolist <- as.numeric(pberr$Gnolist)
pberr$Tnolist <- as.numeric(pberr$Tnolist)

#calc, add cols and sum cols
#(Ano x numreads) + Ad, col for each nt DON'T NEED ANYMORE
#pberr$Apa <- (pberr$Anolist * pberr$numreadslist) + pberr$Adlist #UPDATE DF
#pberr$Cpa <- (pberr$Cnolist * pberr$numreadslist) + pberr$Cdlist #can use num operators on lists? or as.num?
#pberr$Gpa <- (pberr$Gnolist * pberr$numreadslist) + pberr$Gdlist
#pberr$Tpa <- (pberr$Tnolist * pberr$numreadslist) + pberr$Tdlist
#then sum cols

#put adiffs var over sum of col to see pct
#AerrDR6 <- sum(pberr$Adlist)/sum(pberr$Apa) * 100 #adiffs based on 600 000, using all for now
#CerrDR6 <- sum(pberr$Cdlist)/sum(pberr$Cpa) * 100
#GerrDR6 <- sum(pberr$Gdlist)/sum(pberr$Gpa) * 100
#TerrDR6 <- sum(pberr$Tdlist)/sum(pberr$Tpa) * 100 #NUMERATOR DOESN'T USE DIFF NUM SO WRONG, change to sum(Tdlist)
#ANSWERS
#SHOULD BE CALCULATING ERROR RATES FOR CONS

rm(CerrCC) #for var or df

$prob <- NULL #for column

#take out all the seqs where less than 3 reads contrib
pberr2 <- subset(pberr, numreadslist > 2) #LESS SEQS consseq with 3 or more reads contrib

#USING SEQS WITH MORE THAN 2 READS CONTRIBUTING, 400 000
#subset fastq/lines file

#2 types of per base error rates, 3mers, 7mers
#per base error rates considering consensus seqs, 4th graph DIFF OR CONS, ASK WHICH ERROR!! 
#num times seen as error /num times seen

#adiffs/numerrors
#should be??? AerrDC, over 100%/higher num that real?

#acons/numerrors 
#should be acons/totalacount RIGHT
AerrCC4 <- ACcount/sum(pberr2$Anolist) *100 
CerrCC4 <- CCcount/sum(pberr2$Cnolist) *100
GerrCC4 <- GCcount/sum(pberr2$Gnolist) *100
TerrCC4 <- TCcount/sum(pberr2$Tnolist) *100

#per base error rates considering all reads SO HAVE TO MINUS THE CONSENSUS SEQ?? IT IS, pberr2

#diffs, do as b4 with pberr2 df
#AerrDR4 <- sum(pberr2$Adlist)/sum(pberr2$Apa) * 100 
#CerrDR4 <- sum(pberr2$Cdlist)/sum(pberr2$Cpa) * 100
#GerrDR4 <- sum(pberr2$Gdlist)/sum(pberr2$Gpa) * 100
#TerrDR4 <- sum(pberr2$Tdlist)/sum(pberr2$Tpa) * 100 

#THIS IS WHAT I WANT
#cons AerrCR, 
pberr2$AnoR <- pberr2$Anolist*pberr2$numreadslist #DENOMINATOR FOR ERROR RATE AisCONS
pberr2$CnoR <- pberr2$Cnolist*pberr2$numreadslist
pberr2$GnoR <- pberr2$Gnolist*pberr2$numreadslist
pberr2$TnoR <- pberr2$Tnolist*pberr2$numreadslist 

AerrCR4 <- Adcount/sum(pberr2$AnoR) * 100 #sum the above cols for deno
CerrCR4 <- Cdcount/sum(pberr2$CnoR) * 100
GerrCR4 <- Gdcount/sum(pberr2$GnoR) * 100
TerrCR4 <- Tdcount/sum(pberr2$TnoR) * 100

#plot correct error rates? 
  

#3mers, like fig2 in other art, add diff to pcpcol somehow to be x axis, do same for others?
*****
#CHECK ERROR RATE CALC, is numero/deno numbers for error rate for centre cons char wrong?
#meant to consider all reads?
#deno number of times the 3mer seen in cons seqs, it is overlapping atm
#numero, number of times that 3mer is seen where centre char is the conschar of error
#IT IS RIGHT

pcp3D <- subset(flanks, select = c(Prev.Seq, Consensus, Post.Seq, Number, Difference))
#get prev and post to 1 letter
#loop through rows, few/4 mins?
for (i in 1 :nrow(pcp3D)){
  pcp3D$Prev.Seq[i] <- substr(pcp3D$Prev.Seq[i], nchar(pcp3D$Prev.Seq[i]), nchar(pcp3D$Prev.Seq[i]))
  pcp3D$Post.Seq[i] <- substr(pcp3D$Post.Seq[i], 1, 1)
}
#paste, unique
pcp3D$pcpdcol <- paste(pcp3D$Prev.Seq, pcp3D$Consensus, pcp3D$Post.Seq, ">", pcp3D$Difference, sep='')
pcp3D$pcpcol <- paste(pcp3D$Prev.Seq, pcp3D$Consensus, pcp3D$Post.Seq, sep='')
pcp3D$Prev.Seq <- NULL
pcp3D$Consensus <- NULL
pcp3D$Post.Seq <- NULL
pcp3D$Difference <- NULL
pcp3D <- ddply(pcp3,.(pcpdcol),summarize,Number=sum(Number))

#take out 2mers
pcp3Dno2 <- pcp3D #need to keep difference col for new pcp plot
pcp3Dno2$pcpdcol <- gsub('\\s+', '', pcp3Dno2$pcpdcol) #replacing one or more spaces with no space
pcp3Dno2$nchar <- nchar(pcp3Dno2$pcpdcol)
pcp3Dno2 <- subset(pcp3Dno2, nchar == 5, select=c(pcpcol, pcpdcol, Number))

#calc pct, merge
pcp3DP <- merge(pcp3Dno2, combos3, by.x = "pcpcol", by.y = "combos") #merge, these columns have the common row names 
pcp3DP$pct <- NULL
pcp3DP$count2 <- NULL #remove unneccessary cols that came from combos3
pcp3DP$pct2 <- NULL
pcp3DP$pct <- pcp3DP$Number/pcp3DP$count *100 #error rate
pcp3DP <- pcp3DP[rev(order(pcp3DP$pct)),] 

#plot
pcp3Dplot <- ggplot(pcp3DP, aes(x=reorder(pcpdcol, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(pct,1), y=(pct)+0.1)) + labs(title="Proportion of sequence fragments that are the Consensus flanking sequence", x= "Fragment (pre-cons-post)", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(pct)+1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymax = pct, ymin=pct)) # + geom_errorbar(data=pcpP, aes(ymax = pcpP$pct + pcpP$se, ymin= pcpP$pct - pcpP$se), stat = "identity", position = "dodge", width=0.25)

#above is wrong for just considering cons seqs e.g. CTC, T gets 2 counts in py if error is 33T48C1A1   

#consider all reads    
#RUN ANOTHER LOOP
#enviro for number of times see 3mer in cons seqs x num reads for deno DENOMINATOR
combos3R <- do.call(CJ, replicate(3, nts, FALSE))
combos3R <- apply(combos3R, 1, paste, collapse="") #char array/vector
combos3R <- data.frame(combos=combos3R, stringsAsFactors = FALSE) #keys/index for counter

count3mersR <- new.env(parent=emptyenv(), size=nrow(combos3R)*1.2) 

#loop through seqs 50 MINS
for (i in 1:length(lines)){ #for each seq ALL 600,000
  if (i %% 4 == 1) { 
    header <- lines[i]
    #numreads <- substring(header, 51, 53) #strip colons later
    numreads <- strsplit(header, split=":")[[1]] #need num form cause going to multiply later
    length <- length(numreads)
    numreads <- as.numeric(numreads[length-3]) #get 4th last thing in header
  }
  if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
    time <- Sys.time()
    line = lines[i]
    seqno <- (i/4 + 0.5)
    print(seqno)
    len <- nchar(line)
    begin = 1
    endd = 3 #reset for each seq
    if (numreads > 2){
      while (endd <= len){ #200ish loops #for faster than while, if only once
        frag <- substr(line, begin, endd)
        count3mersR[[frag]] <- mget(frag, count3mersR, ifnotfound = 0)[[1]] + numreads
        #print(counters[[frag]]) 
        begin = begin + 1
        endd = endd + 1
        #print(begin)
    }
   } 
  }
}
combos3R$count <- unlist(mget(combos3R$combos, envir = count3mersR, ifnotfound = 0)) 

combos3R$numero <- numeric(192) #creates numeric vector of 64 zeros

#loop through seqs and errors like for pberr
for (i in 1:length(lines)){ #for each line of fastq
  if (i %% 4 == 1) { 
    header <- lines[i]
    #numreads <- substring(header, 51, 53) #strip colons later
    numreads <- strsplit(header, split=":")[[1]] #need num form cause going to multiply later
    length <- length(numreads)
    numreads <- as.numeric(numreads[length-3]) #get 4th last thing in header
  }
  if (numreads > 2 ) {
    if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
      print(i/4 + 0.5) #i for integer, seqno, only prints the seq it goes to
      seq = lines[i]
    }
    if (i %% 4 == 3){
      diffs = lines[i]
      #remove plus start
      diffs <- substring(diffs, 2)
      diffs = strsplit(diffs, " ")[[1]] #creates a list of vector elements 
      #print(diffs)
      errsperseq <- length(diffs)
      #print(errsperseq) 
      if (errsperseq > 0){
        for (k in 1:errsperseq){ #going through errors 
          diff = diffs[k]
          #print(diff)
          #pull out nums and letters
          posno <- as.numeric(strsplit(diff, split="[[:upper:]]")[[1]][1]) #first num, strsplit makes list
          conschar <- substr(seq, posno, posno) #start and stop
          #print(conschar)
          conschar3mer <- substr(seq, posno-1, posno+1) #ADD VAR OR COL WITH DIFF
          #print(conschar3mer)
          #consnum, could take it out and be left with diffnum
          cons <- str_extract(diff, paste0(conschar, "[[:digit:]]+"))
          #consnum <- str_extract(cons, "[[:digit:]]+") #number after cons letter NEED???
          #take out cons from error
          diffwocons <- gsub(cons, "", diff)
          diffnum <-strsplit(diffwocons, split="[[:upper:]]")[[1]][-1] #takes out first num/posno, length 2 OR MORE/5
          #print(diffnum)
          diffchar <-strsplit(diffwocons, split="[[:digit:]]+")[[1]][-1] #plus means any size digit, first one empty for some reason?
          #could be 4 diff chars in one err!
          #so need loop through diffnums
          for (l in 1:nrow(combos3R)) { #don't have to loop?
            if (conschar3mer == combos3R$combos[l]){ #USE COMBOS AS ROW NAMES FOR FASTER? make a list entry, letters, number
              for (j in 1:length(diffnum)){
                combos3R$numero[l] <- combos3R$numero[l] + as.numeric(diffnum[j])
              }
            break  
            }
          }
        }
      }
    }
  }
}
#calc, pct col
combos3R$pct <- combos3R$numero/combos3R$count * 100

#above didn't have diffs included
combos <- vector(mode='list', 4000000) #4 mill, too slow? but preseting, otherwise use enviro DID I OVERWRITE ANOTHER VAR
comboswD <- vector(mode='list', 4000000)
numero <- vector(mode='list', 4000000)
w <- 1
for (i in 1:length(lines)){ #for each line of fastq
  if (i %% 4 == 1) { 
    header <- lines[i]
    #numreads <- substring(header, 51, 53) #strip colons later
    numreads <- strsplit(header, split=":")[[1]] #need num form cause going to multiply later
    length <- length(numreads)
    numreads <- as.numeric(numreads[length-3]) #get 4th last thing in header
  }
  if (numreads > 2 ) {
    if (i %% 4 == 2) { #if line index divided by 4, gets remainder = 2
      print(i/4 + 0.5) #i for integer, seqno, only prints the seq it goes to
      seq = lines[i]
    }
    if (i %% 4 == 3){
      diffs = lines[i]
      #remove plus start
      diffs <- substring(diffs, 2)
      diffs = strsplit(diffs, " ")[[1]] #creates a list of vector elements 
      #print(diffs)
      errsperseq <- length(diffs)
      #print(errsperseq) 
      if (errsperseq > 0){
        for (k in 1:errsperseq){ #going through errors 
          diff = diffs[k]
          #print(diff)
          #pull out nums and letters
          posno <- as.numeric(strsplit(diff, split="[[:upper:]]")[[1]][1]) #first num, strsplit makes list
          conschar <- substr(seq, posno, posno) #start and stop
          #print(conschar)
          conschar3mer <- substr(seq, posno-1, posno+1) #ADD VAR OR COL WITH DIFF
          #print(conschar3mer)
          #consnum, could take it out and be left with diffnum
          cons <- str_extract(diff, paste0(conschar, "[[:digit:]]+"))
          #consnum <- str_extract(cons, "[[:digit:]]+") #number after cons letter NEED???
          #take out cons from error
          diffwocons <- gsub(cons, "", diff)
          diffnum <-strsplit(diffwocons, split="[[:upper:]]")[[1]][-1] #takes out first num/posno, length 2 OR MORE/5
          #print(diffnum)
          diffchar <-strsplit(diffwocons, split="[[:digit:]]+")[[1]][-1] #plus means any size digit, first one empty for some reason?
          #could be 4 diff chars in one err!
          #so need loop through diffnums
          for (j in 1:length(diffnum)){
            cons3merwDiff <- paste0(conschar3mer, ">", diffchar[j])
            #add to the 3 cols
            comboswD[[w]] <- cons3merwDiff
            combos[[w]] <- conschar3mer
            numero[[w]] <- diffnum[j]
            w <- w + 1
          }
        }
      }
    }
  }
}
#dataframe, rid zeros from lists??
combos3errR4wD <- data.frame(cbind(combos, comboswD, numero))

combos <- unlist(combos) #RERUN WITH PLAIN NAMES 
comboswD <- unlist(comboswD)
numero <- as.numeric(unlist(numero))
combos3errR4wDv2 <- data.frame(cbind(combos, comboswD, numero)) #HAS LEVELS
combos3errR4wDv2$numero <- as.numeric(combos3errR4wDv2$numero) #get rid of levels
combos3errR4wDv2$comboswD <- as.character(combos3errR4wDv2$comboswD)
combos3errR4wDv2$combos <- as.character(combos3errR4wDv2$combos)

#combos3errR4wD2 <- data.frame(cbind(unlist(combos), unlist(comboswD), unlist(numero)))#not a list but has levels? UP TO try unique and need to rename

#then unlist, unlist before make columns? unlist will unlist lists within lists to vectors
#combos3errR4wD$combos <- unlist(combos3errR4wD$combos)
#combos3errR4wD$comboswD <- unlist(combos3errR4wD$comboswD)
#combos3errR4wD$numero <- unlist(combos3errR4wD$numero) #error - unlisting removes unfilled/zeros of list

#unique, based on comboswD col
combos3errR4wD <- ddply(combos3errR4wDv2,.(comboswD),summarize,numero=sum(numero)) #works?
#lose combos col 

#char to get rid of 2mers
combos3errR4wD <- subset(combos3errR4wD, nchar(combos3errR4wD$comboswD) == 5)

#add combos col yourself, 3 first letters
combos3errR4wD$combos <- substr(combos3errR4wD$comboswD, 1, 3)

#merge combos 3R and new df
Mcombos3errR4wD <- merge(combos3errR4wD, combos3R, by.x = "combos", by.y = "combos")

#calc pct
Mcombos3errR4wD$pct <- Mcombos3errR4wD$numero/Mcombos3errR4wD$count *100

#plot, order right
combos3R4plot <- ggplot(Mcombos3errR4wD, aes(x=reorder(comboswD, -pct), y=pct)) + geom_bar(stat="identity", width=0.8, fill="#CC3366")+ geom_text(aes(label=round(pct,4), y=(pct)+0.02), angle=90) + labs(title="Error rate of 3-mers", x= "3-mers", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(Mcombos3errR4wD$pct)+0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
# + geom_errorbar(aes(ymax = pct, ymin=pct)) 

#add diff col, order by that
Mcombos3errR4wD$diff <- substr(Mcombos3errR4wD$comboswD, 5, 5)
  
#another version without Ns
#data frame without Ns
Mcombos3errR4wDnoNs <- subset(Mcombos3errR4wD, diff != "N")

#within diff order, order by pct or alphabet? do both
#Mcombos3errR4wDnoNs <- Mcombos3errR4wDnoNs[order(Mcombos3errR4wDnoNs$diff),] #this does diff then alpha
#ggplot won't stick to this

Mcombos3errR4wDnoNs$comboswD <- factor(Mcombos3errR4wDnoNs$comboswD, levels=Mcombos3errR4wDnoNs$comboswD[order(Mcombos3errR4wDnoNs$diff)])

#plot
combos3R4plotnoNs <- ggplot(Mcombos3errR4wDnoNs, aes(x=comboswD, y=pct, fill=diff)) + geom_bar(stat="identity", width=0.8)+ geom_text(aes(label=round(pct,4), y=(pct)+0.06), angle=90) + labs(title="Error rate of 3-mers", x= "3-mers", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(Mcombos3errR4wDnoNs$pct)+0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#ordered by alpha/3mers (after diff)

#diffs be diff colours, should be ordered by cons/central char like in fig2!!!!!!!!!
#add cons col
Mcombos3errR4wDnoNs$cons <- substr(Mcombos3errR4wDnoNs$comboswD, 2, 2)
#order
Mcombos3errR4wDnoNs$comboswD <- factor(Mcombos3errR4wDnoNs$comboswD, levels=Mcombos3errR4wDnoNs$comboswD[order(Mcombos3errR4wDnoNs$cons)])
#plot
combos3R4plotnoNs <- ggplot(Mcombos3errR4wDnoNs, aes(x=comboswD, y=pct, fill=cons)) + geom_bar(stat="identity", width=0.8)+ geom_text(aes(label=round(pct,4), y=(pct)+0.06), angle=90) + labs(title="Error rate of 3-mers", x= "3-mers", y="Percentage (%)")+ theme_bw() + theme(panel.border = element_blank()) + scale_y_continuous(expand=c(0,0), limits = c(0, max(Mcombos3errR4wDnoNs$pct)+0.2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


#after diffs, order by pct


#0.1, low, what were the pcp rates before? denominator was count for just cons seqs
#article has mutation rates (known cancer mutations)
  
#7mers for 400 000, stack, SE 
#all


#front


#back
  


#single change plot


#put R on github


#pink things/move MD




