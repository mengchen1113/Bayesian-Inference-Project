Ec1127 Final Project Code

#Section 1: Treatment assigned

DS<-read.csv(file="Documents/Ec1127/project5data.csv",header=TRUE)
head(DS)

#just for ease of manipulation:

female<-DS$female
age<-DS$age
haschld<-DS$haschld
educ<-DS$educ
currjob<-DS$currjob
mosinjob<-DS$mosinjob
stdearn<-DS$std.earn.yr
white<-DS$white
partnered<-DS$partnered
TreatAssigned<-DS$TreatAssigned
TreatReceived<-DS$TreatReceived
yobs<-DS$WklyEarningsAt4Yrs


#Now condition on gender. 
	#create new vectors:

#dataframe w/ females only:	
DSF<-as.data.frame(DS[female==1,]) #Subgroup after control units with propensity scores outside bounds are removed

DSM<-as.data.frame(DS[female==0,])
		
ageF<-DSF$age
ageM<-DSM$age
haschldF<-haschld[female==1]
haschldM<-haschld[female==0]
educF<-educ[female==1]
educM<-educ[female==0]
currjobF<-currjob[female==1]
currjobM<-currjob[female==0]
mosinjobF<-mosinjob[female==1]
mosinjobM<-mosinjob[female==0]
stdearnF<-stdearn[female==1]
stdearnM<-stdearn[female==0]
whiteF<-white[female==1]
whiteM<-white[female==0]
partneredF<-partnered[female==1]
partneredM<-partnered[female==0]
TreatAssignedF<-TreatAssigned[female==1]
TreatAssignedM<-TreatAssigned[female==0]
TreatReceivedF<-TreatReceived[female==1]
TreatReceivedM<-TreatReceived[female==0]
yobsF<-yobs[female==1]
yobsM<-yobs[female==0]

#add 2nd order & 3rd order interactions. Define as new variables:

intjobF<-mosinjobF*currjobF
int2F<-partneredF*educF*whiteF
intjobM<-mosinjobM*currjobM
int2M<-partneredM*educM*whiteM

#log terms
sqrtageM<-sqrt(ageM)
sqrtmosjobM<-sqrt(mosinjobM)
sqrtageF<-sqrt(ageF)
sqrtmosjobF<-sqrt(mosinjobF)


#get pvalues for male group and female group separately:
#



#for final balance, only add added interaction terms!!!!

mylogit=glm(TreatAssignedM~ageM+haschldM+educM+
currjobM+mosinjobM+stdearnM+whiteM+partneredM+sqrtageM+sqrtmosjobM,family=binomial("logit"))
pscoresM<-fitted.values(mylogit) #pscores corresponding to males

mylogit=glm(TreatAssignedF~ageF+haschldF+educF+
currjobF+mosinjobF+stdearnF+whiteF+partneredF+sqrtageF+sqrtmosjobF,family=binomial("logit"))
pscoresF<-fitted.values(mylogit) #pscores corresponding to females

#1st time design w/ original 8 terms only
hist(pscoresM,xlim=c(0.5,0.75), ylim=c(0,800), main="Histogram of Pscores Conditioning on Gender", density=NULL,xlab="Pscore",col="blue")
par(new=TRUE)
hist(pscoresF,xlim=c(0.5,0.75), ylim=c(0,800),main="",xlab="",col="orange")
legend("topright", title="Legend", cex=0.75, pch=16, 
  col=c("blue", "orange"), density=TRUE, legend=c("Male", "Female"), ncol=2)
  
#Males: pscores by treatment group 

pscoreMT<-rep(NA,0)
pscoreMC<-rep(NA,0)

for(ii in 1:3863) #population size
{
	if(TreatAssignedM[ii]==1)
	{pscoreMT<-rbind(pscoreMT,pscoresM[ii])}
		else
		{pscoreMC<-rbind(pscoreMC,pscoresM[ii])}	
}


#Females: pscores by treatment group

pscoreFT<-rep(NA,0)
pscoreFC<-rep(NA,0)

for(ii in 1:2616) #population size
{
	if(TreatAssignedF[ii]==1)
	{pscoreFT<-rbind(pscoreFT,pscoresF[ii])}
		else
		{pscoreFC<-rbind(pscoreFC,pscoresF[ii])}
}

par(mfrow=c(2,1)) #create a better looking histogram this time

hist(pscoreMT, xlim=c(0.5,0.8), ylim=c(0,800), main="Pscores with Transformations by Treatment Group for Males", density=FALSE,xlab="Pscore",col="blue")
par(new=TRUE)
hist(pscoreMC,xlim=c(0.5,0.8), ylim=c(0,800),main="",xlab="",col="red")

hist(pscoreFT,xlim=c(0.5,0.8), ylim=c(0,600), main="Pscores with Transformations by Treatment Group for Females", density=FALSE,xlab="Pscore",col="blue")
par(new=TRUE)
hist(pscoreFC,xlim=c(0.5,0.8), ylim=c(0,600), main="",xlab="",col="red")
legend("topright", title="Legend", cex=0.75, pch=16, 
col=c("blue", "red"), legend=c("Treated", "Control"), ncol=2)


#Throw out control group units whose pscores are outside the range of pscores for treated group units:


maxPMT<-max(pscoreMT) #range of propensity scores for males
minPMT<-min(pscoreMT)

maxPFT<-max(pscoreFT) 
minPFT<-min(pscoreFT)


DSM2<-as.data.frame(DSM[pscoresM>=minPMT & pscoresM <= maxPMT,]) #Subgroup after control units with propensity scores outside bounds are removed

length(DSM$female)
[1] 3863
length(DSM2$female)
[1] 3861 #note that only 2 units were removed after this process
#instead of redoing results, assume these two units do not significantly diminish robustness of findings

DSF2<-as.data.frame(DSF[pscoresF>=minPFT & pscoresF <= maxPFT,])


#subclassify again using quintiles:
#values shown are from rebalancing w/ transformations

quantile(pscoresF,c(.2,.4,.6,.8))#for females
      20%       40%       60%       80% 
0.6403928 0.6591047 0.6730441 0.6929311 

Q1F<-0.6403928
Q2F<-0.6591047
Q3F<-0.6730441
Q4F<-0.6929311

quantile(pscoresM,c(.2,.4,.6,.8))#for males
      20%       40%       60%       80% 
0.5649519 0.5757407 0.5890158 0.6031824 

Q1M<-0.5649519
Q2M<-0.5757407
Q3M<-0.5890158
Q4M<-0.6031824

subIndM<-rep(NA,3863)
subIndF<-rep(NA,2616)
 
#to make subclasses for units:
 for(ii in 1:3863)
 {
 	if(pscoresM[ii]<=Q1M)
 	{subIndM[ii]<-1}else
 	
 	if(pscoresM[ii]<=Q2M & pscoresM[ii]>Q1M)
 	{subIndM[ii]<-2}else
 	
 	if(pscoresM[ii]<=Q3M & pscoresM[ii]>Q2M)
 	{subIndM[ii]<-3}else
 	
 	if(pscoresM[ii]<=Q4M & pscoresM[ii]>Q3M)
 	{subIndM[ii]<-4}else
 	
 	{subIndM[ii]<-5}
 	
}

 for(ii in 1:2616)
 {
 	if(pscoresF[ii]<=Q1F)
 	{subIndF[ii]<-1}else
 	
 	if(pscoresF[ii]<=Q2F & pscoresF[ii]>Q1F)
 	{subIndF[ii]<-2}else
 	
 	if(pscoresF[ii]<=Q3F & pscoresF[ii]>Q2F)
 	{subIndF[ii]<-3}else
 	
 	if(pscoresF[ii]<=Q4F & pscoresF[ii]>Q3F)
 	{subIndF[ii]<-4}else
 	
 	{subIndF[ii]<-5}
 	
}

#check length in subclasses:

SCsizeFT<-rep(NA,5)
SCsizeFC<-rep(NA,5)
SCsizeMT<-rep(NA,5)
SCsizeMC<-rep(NA,5)

for(ii in 1:5)
{
	SCsizeFT[ii]<-length(subIndF[subIndF==ii & TreatAssignedF==1])
	SCsizeFC[ii]<-length(subIndF[subIndF==ii & TreatAssignedF==0])
	SCsizeMT[ii]<-length(subIndM[subIndM==ii & TreatAssignedM==1])
	SCsizeMC[ii]<-length(subIndM[subIndM==ii & TreatAssignedM==0])
	}
	

	
	
#find means for subclasses for each gender stratum:

meancovFT<-matrix(data=NA,nrow=10,ncol=5) #want empty 12x5 matrices
meancovFC<-matrix(data=NA,nrow=10,ncol=5) 
meancovMT<-matrix(data=NA,nrow=10,ncol=5) 
meancovMC<-matrix(data=NA,nrow=10,ncol=5) 



for(ii in 1:5)
{
	#create a column of mean covariates for each subclass
	
	#female treatment assigned group
		meancovFT[,ii]<-rbind(mean(ageF[subIndF==ii & TreatAssignedF==1]),mean(haschldF[subIndF==ii & TreatAssignedF==1]),mean(educF[subIndF==ii & TreatAssignedF==1]),mean(currjobF[subIndF==ii & TreatAssignedF==1]),mean(mosinjobF[subIndF==ii & TreatAssignedF==1]),mean(stdearnF[subIndF==ii & TreatAssignedF==1]),mean(whiteF[subIndF==ii & TreatAssignedF==1]),mean(partneredF[subIndF==ii & TreatAssignedF==1]),mean(sqrtageF[subIndF==ii & TreatAssignedF==1]),mean(sqrtmosjobF[subIndF==ii & TreatAssignedF==1]))
	
	#female control assigned group
	meancovFC[,ii]<-rbind(mean(ageF[subIndF==ii & TreatAssignedF==0]),mean(haschldF[subIndF==ii & TreatAssignedF==0]),mean(educF[subIndF==ii & TreatAssignedF==0]),mean(currjobF[subIndF==ii & TreatAssignedF==0]),mean(mosinjobF[subIndF==ii & TreatAssignedF==0]),mean(stdearnF[subIndF==ii & TreatAssignedF==0]),mean(whiteF[subIndF==ii & TreatAssignedF==0]),mean(partneredF[subIndF==ii & TreatAssignedF==0]),mean(sqrtageF[subIndF==ii & TreatAssignedF==0]),mean(sqrtmosjobF[subIndF==ii & TreatAssignedF==0]))
	
	
	#male treatment assigned group
	meancovMT[,ii]<-rbind(mean(ageM[subIndM==ii & TreatAssignedM==1]),mean(haschldM[subIndM==ii & TreatAssignedM==1]),mean(educM[subIndM==ii & TreatAssignedM==1]),mean(currjobM[subIndM==ii & TreatAssignedM==1]),mean(mosinjobM[subIndM==ii & TreatAssignedM==1]),mean(stdearnM[subIndM==ii & TreatAssignedM==1]),mean(whiteM[subIndM==ii & TreatAssignedM==1]),mean(partneredM[subIndM==ii & TreatAssignedM==1]),mean(sqrtageM[subIndM==ii & TreatAssignedM==1]),mean(sqrtmosjobM[subIndM==ii & TreatAssignedM==1]))
	
	#male control assigned group
	meancovMC[,ii]<-rbind(mean(ageM[subIndM==ii & TreatAssignedM==0]),mean(haschldM[subIndM==ii & TreatAssignedM==0]),mean(educM[subIndM==ii & TreatAssignedM==0]),mean(currjobM[subIndM==ii & TreatAssignedM==0]),mean(mosinjobM[subIndM==ii & TreatAssignedM==0]),mean(stdearnM[subIndM==ii & TreatAssignedM==0]),mean(whiteM[subIndM==ii & TreatAssignedM==0]),mean(partneredM[subIndM==ii & TreatAssignedM==0]),mean(sqrtageM[subIndM==ii & TreatAssignedM==0]),mean(sqrtmosjobM[subIndM==ii & TreatAssignedM==0]))
	
	}
	

	
		#Some GRAFXXXXXX: Histograms for each covariate by treatment assigned 
	#For females (repeat for males):
	
 #create a better looking histogram this time

#Subclass n (changing n each time to get histograms of all 8 covariates for each subclass)

#Subclass 2: 
n=5
par(mfrow=c(2,4))
#X1
hist(ageM[subIndM==n & TreatAssignedM==1], ylim=c(0,150), xlim=c(15,25),main="Age Distribution", density=FALSE,xlab="Age",col="blue") #treated units
par(new=TRUE)
hist(ageM[subIndM==n & TreatAssignedM==0], ylim=c(0,150), xlim=c(15,25),main="", xlab="", density=FALSE,col="red") #control units
 
#X2
hist(haschldM[subIndM==n & TreatAssignedM==1], ylim=c(0,450), main="Distribution of Having a Child", density=FALSE,xlab="Indicator",col="blue") #treated units
par(new=TRUE)
hist(haschldM[subIndM==n & TreatAssignedM==0], ylim=c(0,450), main="", xlab="",density=FALSE,col="red") #control units

#X3
hist(educM[subIndM==n & TreatAssignedM==1], ylim=c(0,400), main="Education", density=FALSE,xlab="Indicator",col="blue") #treated units
par(new=TRUE)
hist(educM[subIndM==n & TreatAssignedM==0], ylim=c(0,400), main="", xlab="",density=FALSE,col="red") #control units

#X4
hist(currjobM[subIndM==n & TreatAssignedM==1], ylim=c(0,500), main="Whether currently holding job", density=FALSE,xlab="Indicator",col="blue") #treated units
par(new=TRUE)
hist(currjobM[subIndM==n & TreatAssignedM==0], ylim=c(0,500), main="", xlab="",density=FALSE,col="red") #control units
 
#5
hist(mosinjobM[subIndM==n & TreatAssignedM==1], ylim=c(0,150), main="Months in Job", density=FALSE,xlab="Months",col="blue") #treated units
par(new=TRUE)
hist(mosinjobM[subIndM==n & TreatAssignedM==0], ylim=c(0,150), main="", xlab="",density=FALSE,col="red") #control units
 
#6
hist(stdearnM[subIndM==n & TreatAssignedM==1], ylim=c(0,350), xlim=c(-1,7), main="Standardized yearly earnings", density=FALSE,xlab="Standardized Earnings",col="blue") #treated units
par(new=TRUE)
hist(stdearnM[subIndF==n & TreatAssignedM==0], ylim=c(0,350), xlim=c(-1,7), main="", xlab="",density=FALSE,col="red") #control units
 
#7
hist(whiteM[subIndM==n & TreatAssignedM==1], ylim=c(0,500), main="White", density=FALSE,xlab="Indicator",col="blue") #treated units
par(new=TRUE)
hist(whiteM[subIndM==n & TreatAssignedM==0], ylim=c(0,500), main="", xlab="",density=FALSE,col="red") #control units
 
#8
hist(partneredM[subIndM==n & TreatAssignedM==1], ylim=c(0,500), main="Partnered", density=FALSE,xlab="Indicator",col="blue") #treated units
par(new=TRUE)
hist(partneredM[subIndM==n & TreatAssignedM==0], ylim=c(0,500), main="", xlab="",density=FALSE,col="red") #control units
legend("topright", cex=0.75, pch=16, 
+ col=c("blue", "red"), density=FALSE, legend=c("Treated", "Control"), ncol=2)


#get 4 weight vectors: proportion of control units in each stratum subclass, "" treated (each will be 5 values):

#Subclass size in final design w/ 2 log terms:
SCsizeFT
[1] 321 352 331 375 367
> SCsizeFC
[1] 203 170 193 148 156
> SCsizeMT
[1] 417 438 452 468 482
> SCsizeMC
[1] 356 335 320 304 291



weightMT<-SCsizeMT/2257
weightMC<-SCsizeMC/1606
weightFT<-SCsizeFT/1746
weightFC<-SCsizeFC/870


#now can calculate estimated weighted mean covariate values:

tempMT<-weightMT*meancovMT
tempMC<-weightMC*meancovMC
tempFT<-weightFT*meancovFT
tempFC<-weightFC*meancovFC





wtmeancovMT<-rep(NA,10) #each row is a covariate- we'll want the weighted sum for each covariate over all subclasses (columns)
wtmeancovMC<-rep(NA,10) 
wtmeancovFT<-rep(NA,10)
wtmeancovFC<-rep(NA,10)

for(n in 1:10)
{
	wtmeancovMT[n]<-sum(tempMT[n,])
	wtmeancovMC[n]<-sum(tempMC[n,])
	wtmeancovFT[n]<-sum(tempFT[n,])
	wtmeancovFC[n]<-sum(tempFC[n,])
	
	}


#create a loveplot:

#let y be subclass 1-10 and x be difference in covariate means


covMdiffs<-abs(tempMT-tempMC) #10x5 matrices, note that columns are subclasses!
covFdiffs<-abs(tempFT-tempFC)

covdiffs<-cbind(covMdiffs,covFdiffs) #now 10x10 matrix, still columns are subclasses



for(ii in 1:10)
{
	y<-rep(ii,10)
	plot(covdiffs[,ii],y,xlim=c(0,0.3),ylim=c(0,11),col=c("red","pink","blue","brown","black","purple","yellow","green","orange","gray"),main="Assessment of Covariate Balance",xlab="Mean difference in covariate values in treatment groups",ylab="subclass")
	par(new=TRUE)
	}
	
	abline(h=c(1:10))
	
	legend("topright", cex=0.5, pch=10, col=c("red","pink","blue","brown","black","purple","yellow","green","orange","gray"), density=FALSE, legend=c("age", "child","education","current job","months in job","yrly earn","white","partnered","sqrt mos in job", "sqrt age"), ncol=5)

#Weight vectors for proportion of population in each of the 10 subclasses 

wtsSC<-c(SCsizeMT+SCsizeMC,SCsizeFT+SCsizeFC)/6479 #create weight vector-proportion of total population in each subclass 1-10:

wtsSC #w/ log terms: 
 [1] 0.11930854 0.11930854 0.11915419
 [4] 0.11915419 0.11930854 0.08087668
 [7] 0.08056799 0.08087668 0.08072233
[10] 0.08072233

#to estimate 

#to estimate average causal effect:
meandiffsM<-rep(NA,5)
meandiffsF<-rep(NA,5)

for(ii in 1:5)
{
	meandiffsM[ii]<-mean(yobsM[subIndM==ii & TreatAssignedM==1])-mean(yobsM[subIndM==ii & TreatAssignedM==0])
	meandiffsF[ii]<-mean(yobsF[subIndF==ii & TreatAssignedF==1])-mean(yobsF[subIndF==ii & TreatAssignedF==0])
	}
	

> meandiffsM #w/o sqrt terms:
[1] 22.399528 21.793578 35.399669  3.527369 13.859094
> meandiffsF
[1] 14.837343  3.336152 18.669302 16.419690 28.701413


#w/ transformations:

> meandiffsM
[1] 34.77062 20.89613 17.64614 10.51337 11.75131
> meandiffsF
[1] 26.622029 24.965784 12.618011  5.378497 21.323408
> 

colMdiffs<-c(meandiffsM,meandiffsF)
colMdiffs
 [1] 34.770621 20.896128 17.646137 10.513369 11.751306
 [6] 26.622029 24.965784 12.618011  5.378497 21.323408

 
 #Now apply weights:

meandiffsY<-sum(colMdiffs*wtsSC)
meandiffsY
[1] 18.73936

#calculate weighted variance:

wtdvar<-rep(NA,10) #finding variance w/in each subclass


for(ii in 1:5)
{
	#ii=1-5 is for male subclasses
	wtdvar[ii]<-(wtsSC[ii])^2 *((1/length(yobsM[subIndM==ii & TreatAssignedM==1]))*var(yobsM[subIndM==ii & TreatAssignedM==1])+(1/length(yobsM[subIndM==ii & TreatAssignedM==0]))*var(yobsM[subIndM==ii & TreatAssignedM==0]))
	
	#for females, add 5 to index:
	wtdvar[ii+5]<-(wtsSC[ii+5])^2 *((1/length(yobsF[subIndF==ii & TreatAssignedF==1]))*var(yobsF[subIndF==ii & TreatAssignedF==1])+(1/length(yobsF[subIndF==ii & TreatAssignedF==0]))*var(yobsF[subIndF==ii & TreatAssignedF==0])) 



if(ii==5)
	{print(c("weighted variance=",sum(wtdvar)))}
	
}

[1] "weighted variance=" "29.5134904496964"  



wSD<-sqrt(sum(wtdvar))



#95% Neyman conf int

CInterval<-c(meandiffsY-1.96*wSD, meandiffsY+1.96*wSD)

CInterval
[1]  8.091401 29.387319 #virtually the same the second time.


#Section 2: Treatment received

#we have TreatReceivedF and TreatReceivedM already
#need CACE estimate within each subclass, then averaged over subclasses:

#Have: wtsSC, colMdiffs- which is the subclass ITT estimates
#Want: proportion compliers=1-proportion noncompliers per subclass

#For each male subclass n:

n=1
A<-length(TreatAssignedM[TreatAssignedM==1 & TreatReceivedM==0 & subIndM==n])
B<-length(TreatAssignedM[subIndM==n])

#A/B = proportion noncompliers per subclass

propCM<-rep(NA,5) #prop noncompliers for males
propCF<-rep(NA,5) #prop noncompliers for females

totCACE<-rep(NA,10) 

#have colMdiffs already for subclasses 1-10, where 1-5 males, 6-10 females (both increasing in pscores)

for(n in 1:5)
{

A<-length(TreatAssignedM[TreatAssignedM==1 & TreatReceivedM==0 & subIndM==n]) # total # noncompliers in subclass n for males

B<-length(TreatAssignedM[subIndM==n]) # total units in subclass n for males

C<-length(TreatAssignedF[TreatAssignedF==1 & TreatReceivedF==0 & subIndF==n]) #total # noncompliers in subclass n for females

D<-length(TreatAssignedM[subIndF==n]) #total units in subclass n

propCM[n]<-1-A/B #prop compliers for males
propCF[n]<-1-C/D #"" for females

totCACE[n]<-colMdiffs[n]/propCM[n]
totCACE[n+5]<-colMdiffs[n+5]/propCF[n]

}

propC<-c(propCM,propCF)

[1] 0.8783959 0.8706339 0.8626943 0.8665803 0.8382924 0.8894668 0.8740360 0.8613990 0.8359788 0.8324873

totCACE
[1] 39.584226 24.001050 20.454682 12.132019 14.018148 29.930322 28.563795 14.648277  6.433772 25.614094


#compare with corresponding ITT:

[1] 34.770621 20.896128 17.646137 10.513369 11.751306 26.622029 24.965784 12.618011  5.378497 21.323408 #slightly lower per, which makes sense

CACE wtd estimate: 

wtdCACE<-sum(totCACE*wtsSC)
[1] 21.63527

ITTy<-colMdiffs

#Now calculate the variance using the delta method:

#for each subclass, we have CACE, ITT, and P_complier estimates:
#CACE=totCACE, ITT=colMdiffs=ITTy, Proportion complier=propC



vCACE<-rep(NA,10) #empty matrix for subclass var(CACE) 

#calculate for men first (indexing issue...)
for(n in 1:5)
{
	A<-ITTy[n]
	B<-propC[n]
	#cov(x,y)<-1/
	vCACE[n]<-getvarCACEM(a=A,b=B,m=n)
	
	}

#then for females

for(n in 1:5)
{
	A<-ITTy[n+5]
	B<-propC[n+5]
	#cov(x,y)<-1/
	
	vCACE[n+5]<-getvarCACEF(a=A,b=B,m=n)
	
	}

#first for males (since I indexed subsets of data using 1-5, easier this way to create separate functions per stratum)

getvarCACEM<-function(a,b,m)
{
	
	
	Nt<-length(TreatAssignedM[TreatAssignedM==1 & subIndM==m]) #treated group subclass population
	Nc<-length(TreatAssignedM[TreatAssignedM==0 & subIndM==m]) #control group subclass population
	
	#to get cov(ITTy,Pc), call ITTy=X and Pc=Y
	Ytdiff<-yobsM[TreatAssignedM==1 & subIndM==m]-mean(yobsM[TreatAssignedM==1 & subIndM==m]) #diff btwn unit's earnings and mean earnings- generally quite negative
	Wtdiff<-TreatReceivedM[TreatAssignedM==1 & subIndM==m]-mean(TreatReceivedM[TreatAssignedM==1 & subIndM==m]) #diffs btwn unit indicator of whether treatment was received and the average W in the subclass 
	
	Ycdiff<-yobsM[TreatAssignedM==0 & subIndM==m]-mean(yobsM[TreatAssignedM==0 & subIndM==m])
	Wcdiff<-TreatReceivedM[TreatAssignedM==0 & subIndM==m]-mean(TreatReceivedM[TreatAssignedM==0 & subIndM==m])
	
	varYWt<-sum(Ytdiff*Wtdiff)
	varYWc<-sum(Ycdiff*Wcdiff)
	
	#now can plug values into formula for cov(ITTy,Pc)
	covXY<-(1/(Nt*(Nt-1)))*varYWt+(1/(Nc*(Nc-1)))*varYWc
	
	#To get sample variances within treatment groups for treatment assigned:
	svarTX<-var(yobsM[TreatAssignedM==1 & subIndM==m])
	svarCX<-var(yobsM[TreatAssignedM==0 & subIndM==m])
	
	#To get sample variances within treatment groups for treatment received:
	svarTY<-var(TreatReceivedM[TreatAssignedM ==1 & subIndM==m])
	svarCY<-var(TreatReceivedM[TreatAssignedM ==0 & subIndM==m])
	
	#Can now plug in above values for Neyman variances- overestimates
	
	varX<-svarTX/Nt+svarCX/Nc 
	varY<-svarTY/Nt+svarCY/Nc	
	
	#overall subclass variance, using equation derived from the variance of the Taylor approx of ITTy/Pc evaluated at the point estimates a, b in each subclass 
	varSC<-(1/b^2)*varX+(a^2/b^4)*varY-(2*a/b^3)*covXY
	
	
	return(varSC)
	
	} 


#for females, separate due to indexing issues

getvarCACEF<-function(a,b,m)
{
	
	Nt<-length(TreatAssignedF[TreatAssignedF==1 & subIndF==m]) #treated group subclass population
	Nc<-length(TreatAssignedF[TreatAssignedF==0 & subIndF==m]) #control group subclass population
	
	#to get cov(ITTy,Pc), call ITTy=X and Pc=Y
	Ytdiff<-yobsF[TreatAssignedF==1 & subIndF==m]-mean(yobsF[TreatAssignedF==1 & subIndF==m])
	Wtdiff<-TreatReceivedF[TreatAssignedF==1 & subIndF==m]-mean(TreatReceivedF[TreatAssignedF==1 & subIndF==m])
	
	Ycdiff<-yobsF[TreatAssignedF==0 & subIndF==m]-mean(yobsF[TreatAssignedF==0 & subIndF==m])
	Wcdiff<-TreatReceivedF[TreatAssignedF==0 & subIndF==m]-mean(TreatReceivedF[TreatAssignedF==0 & subIndF==m])
	
	varYWt<-sum(Ytdiff*Wtdiff)
	varYWc<-sum(Ycdiff*Wcdiff)
	
	#now can plug values into formula for cov(ITTy,Pc)
	covXY<-(1/(Nt*(Nt-1)))*varYWt+(1/(Nc*(Nc-1)))*varYWc
	
	#To get sample variances within treatment groups for treatment assigned:
	svarTX<-var(yobsF[TreatAssignedF==1 & subIndF==m])
	svarCX<-var(yobsF[TreatAssignedF==0 & subIndF==m])
	
	#To get sample variances within treatment groups for treatment received:
	svarTY<-var(TreatReceivedF[TreatAssignedF ==1 & subIndF==m])
	svarCY<-var(TreatReceivedF[TreatAssignedF ==0 & subIndF==m])
	
	#Can now plug in above values for Neyman variances- overestimates
	
	varX<-svarTX/Nt+svarCX/Nc 
	varY<-svarTY/Nt+svarCY/Nc	
	
	#overall subclass variance, using equation derived from the variance of the Taylor approx of ITTy/Pc evaluated at the point estimates a, b in each subclass 
	varSC<-(1/b^2)*varX+(a^2/b^4)*varY-(2*a/b^3)*covXY
	
	
	return(varSC)
	
	} 
	
	#var(CACE) for subclasses 1-10
	
	vCACE
 [1] 360.8623 336.3588 389.2671 308.9306 426.3738 288.4228 463.9767
 [8] 364.9020 687.1064 354.9446
	
	
	#get weighted var(CACE):
	
	wtdvCACE<-sum(vCACE*wtsSC)
	
wtdvCACE
[1] 391.5848

sdCACE<-sqrt(wtdvCACE)

cInt<-c(wtdCACE-1.96*sdCACE,wtdCACE+1.96*sdCACE)

cInt
[1] -17.15020  60.42073





#Complier Types:

>length(TreatAssigned[TreatAssigned==0 & TreatReceived==1])
[1] 0

0 Always takers-> 1-sided non-compliance

> length(TreatAssigned[TreatAssigned==0 & TreatReceived==0])
[1] 2476 N,C

> length(TreatAssigned[TreatAssigned==1]) & TreatReceived==0])
[1] 1074 N

> length(TreatAssigned[TreatAssigned==1 & TreatReceived==1])
[1] 2929 C 










#	Covariate means w/ transformations

	> meancovMT
             [,1]         [,2]        [,3]        [,4]         [,5]
 [1,] 18.37486254 17.940082374 18.40528569 19.17090494 20.867327718
 [2,]  0.13429257  0.054794521  0.07079646  0.09188034  0.186721992
 [3,]  0.12470024  0.146118721  0.44469027  0.53418803  0.784232365
 [4,]  0.04316547  0.006849315  0.11061947  0.42521368  0.688796680
 [5,]  2.20413043  3.164600027  4.22022636  6.32113746  7.275319908
 [6,] -0.25948895 -0.143715089  0.02780518  0.41723112  0.909428937
 [7,]  0.52997602  0.296803653  0.38716814  0.27991453  0.265560166
 [8,]  0.22062350  0.020547945  0.01769912  0.01495726  0.008298755
 [9,]  4.28351574  4.231815616  4.28474411  4.37288867  4.559375514
[10,]  0.85681015  1.302832196  1.68690526  2.31570984  2.504200977
> meancovMC
             [,1]        [,2]       [,3]        [,4]        [,5]
 [1,] 18.28774048 17.99996737 18.4079450 19.13879806 20.89564570
 [2,]  0.13483146  0.06865672  0.0750000  0.08552632  0.17525773
 [3,]  0.10955056  0.16716418  0.4312500  0.52631579  0.82817869
 [4,]  0.04213483  0.01791045  0.1468750  0.45065789  0.63917526
 [5,]  2.23284606  2.91574159  4.8099404  6.35604579  7.10770115
 [6,] -0.26261004 -0.16177789  0.1202729  0.46683130  0.81997708
 [7,]  0.52808989  0.25970149  0.4250000  0.27302632  0.26804124
 [8,]  0.19943820  0.03283582  0.0281250  0.01315789  0.01718213
 [9,]  4.27349073  4.23873513  4.2852617  4.36973514  4.56235238
[10,]  0.86266377  1.21872190  1.8622067  2.31612007  2.49461398
> meancovFT
             [,1]        [,2]        [,3]        [,4]        [,5]
 [1,] 17.71424072 18.54575028 18.65370227 19.50380277 21.68367736
 [2,]  0.39875389  0.32102273  0.26586103  0.28266667  0.36239782
 [3,]  0.49844237  0.51420455  0.43202417  0.49066667  0.60762943
 [4,]  0.36760125  0.27272727  0.20241692  0.20000000  0.19073569
 [5,]  4.59570360  4.97056290  4.19061897  3.90021317  4.34431090
 [6,]  0.11139133  0.03540664 -0.03452061 -0.10761316  0.03891576
 [7,]  0.02492212  0.16761364  0.20845921  0.30666667  0.41689373
 [8,]  0.00623053  0.01704545  0.03625378  0.07466667  0.29700272
 [9,]  4.20609035  4.30283460  4.31293685  4.41056562  4.65143683
[10,]  1.95881146  1.93484465  1.49938117  1.34350024  1.41662304
> meancovFC
             [,1]          [,2]        [,3]        [,4]       [,5]
 [1,] 17.91043931 18.4563416471 18.47821767 19.63461669 21.5449241
 [2,]  0.45320197  0.2764705882  0.25906736  0.26351351  0.3717949
 [3,]  0.55665025  0.5058823529  0.44041451  0.48648649  0.5512821
 [4,]  0.28571429  0.2705882353  0.24870466  0.19594595  0.2371795
 [5,]  4.56110927  4.9162031485  4.66263638  4.18652483  3.6238745
 [6,]  0.15666316 -0.0005485901  0.05361593 -0.04415497 -0.1428519
 [7,]  0.03940887  0.1705882353  0.25388601  0.28378378  0.3653846
 [8,]  0.01477833  0.0470588235  0.05181347  0.08783784  0.2051282
 [9,]  4.22880802  4.2928217257  4.29259754  4.42546759  4.6364778
[10,]  1.96359386  1.9549595266  1.58908551  1.45407686  1.1908422
	
	