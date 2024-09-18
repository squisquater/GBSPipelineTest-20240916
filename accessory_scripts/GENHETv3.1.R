#!/usr/bin/env Rscript

# Load necessary libraries
library(gtools)

# Read in the parameters from Snakemake
genhet_input_file <- snakemake@input$genhet
locus_input_file <- snakemake@input$loci
output_file <- snakemake@output$genhet_output
estimfreq <- snakemake@params$estimfreq

# Load the input data
dat <- read.table(genhet_input_file, sep="\t", header = FALSE)
locusname = scan(locus_input_file,what="character",sep="\t")

#####################################################################################
############################ GENHET FUNCTION #############################
#####################################################################################
"GENHET"<-
function(dat,estimfreq,locname,alfuser){

nbloc=(ncol(dat)-1)/2
nbind=nrow(dat)

#estimation of allele frequencies (only if estimfreq=T)

if(estimfreq=="T")

{

#creation of the list of alleles
datv=vector(length=nbind*nbloc*2)
for (i in 2:ncol(dat)) datv[(nrow(dat)*(i-2)+1):(nrow(dat)*(i-1))]=dat[,i]
al=sort(na.omit(unique(datv)))

#count of the number of times each allele appears + nb of missing data
alcount=matrix(nrow=(length(al)+1),ncol=(nbloc+1))
alcount[,1]=c(al,NA)
for(j in 1:(nrow(alcount)-1))
	for(k in 1:(ncol(alcount)-1))
		alcount[j,(k+1)]=sum(dat[,(k*2):(k*2+1)]==alcount[j,1],na.rm=T)
for(l in 2:ncol(alcount))
		alcount[nrow(alcount),l]=(2*nbind-sum(alcount[1:(nrow(alcount)-1),l]))


#creation of the table of allele frequencies
alfreq=matrix(nrow=length(al),ncol=(nbloc+1))
colnames(alfreq)=c("Allele",locname)
alfreq[,1]=al
for(m in (1:nrow(alfreq)))
	for (n in 2:ncol(alfreq)) alfreq[m,n]=alcount[m,n]/(nbind*2-alcount[nrow(alcount),n])

}

else alfreq=alfuser


dat=as.data.frame(dat)
library(gtools)
res=matrix(nrow=nrow(dat),ncol=6)
colnames(res)=c("sampleid","PHt","Hs_obs","Hs_exp","IR","HL")
res[,1]=as.character(dat[,1])



#estimation of E per locus (for HL and Hs_exp)
E=vector(length=nbloc)
alfreq2=alfreq[,2:ncol(alfreq)]*alfreq[,2:ncol(alfreq)]
for(k in 1:ncol(alfreq2)) E[k]=1-sum(alfreq2[,k],na.rm=T)


#estimation of the mean heterozygosity per locus
mHtl=vector(length=nbloc)
ctNAl=0
ctHtl=0
for(l in 1:ncol(dat))
	{if (even(l)==T)
	 	{
		for (m in 1:nrow(dat))
			{ if (is.na(dat[m,l])==T) ctNAl=(ctNAl+1)
				  else if (is.na(dat[m,(l+1)])==T) ctNAl=(ctNAl+1)
			      	   else if (dat[m,l]!=dat[m,(l+1)]) ctHtl=(ctHtl+1)
			}
		mHtl[l/2]=ctHtl/(nrow(dat)-ctNAl)
		ctNAl=0
		ctHtl=0
		}
	}

#the program in itself

ctHt=0
ctNA=0

ctHm=0
smHtl=0
mmHtl=0
sE=0
mE=0

sfl=0

sEh=0
sEj=0

for(i in 1:nrow(dat))
	{ for (j in 2:(nbloc*2))
		{ if (even(j)==T)
			  {
			  if (is.na(dat[i,j])==T) ctNA=(ctNA+1)
			  else if (is.na(dat[i,(j+1)])==T) ctNA=(ctNA+1)
			       else {
			       		 if (dat[i,j]!=dat[i,(j+1)])
			       			{
			       		 	 ctHt=(ctHt+1)
			       		 	 sEj=sEj+E[j/2]
			       			}
			       		else sEh=sEh+E[j/2]
			       		smHtl=smHtl+mHtl[j/2]
			       		sE=sE+E[j/2]
			       		sfl=sfl+alfreq[alfreq[,1]==as.numeric(dat[i,j]),(j/2+1)]+alfreq[alfreq[,1]==as.numeric(dat[i,j+1]),(j/2+1)]
			       		}
			  }
		}
	  res[i,2]=ctHt/(nbloc-ctNA)
	  mmHtl=smHtl/(nbloc-ctNA)
	  res[i,3]=(ctHt/(nbloc-ctNA))/mmHtl
	  mE=sE/(nbloc-ctNA)
	  res[i,4]=(ctHt/(nbloc-ctNA))/mE
	  ctHm=nbloc-ctHt-ctNA
	  res[i,5]=(2*ctHm-sfl)/(2*(nbloc-ctNA)-sfl)
	  res[i,6]=sEh/(sEh+sEj)
	  ctHt=0
	  ctNA=0
	  ctHm=0
	  smHtl=0
	  mmHtl=0
	  sE=0
	  mE=0
	  sfl=0
	  sEh=0
	  sEj=0
	  }
return(res)
}

#####################################################################################
#####################################################################################

# Run the simplified GENHET function 
cat("Running GENHET analysis...\n")
res <- GENHET(dat=dat,estimfreq=estimfreq,locname=locusname)

# Write the output to a file
write.table(res, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)

# Final message
cat("GENHET analysis completed and output saved to:", output_file, "\n")
