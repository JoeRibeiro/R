an <- as.numeric; af <- as.factor; ac <- as.character
h <- head; u <- unique
count<- function(DATA){length(unique(DATA))}
neg<-function(x) -x 

#tapply.ID is rather like aggregate, but I wrote it before I discovered "aggregate" and actually it does a few more things!
#tapply.ID takes a data frame, a vector of the names of the data column, a vector of the factors to apply across, the name of the function to apply (i.e. sum) and the name of the new variable
#we then create an ID field in the dataframe and use that to run the tapply on.  this way we don't try to expand out to include missing combinations of strata.

#returns a dataframe with the new value
tapply.ID <- function(df,datacols, factorcols, func, newnames, pad, padval, na.stuff, func.args)
  {

    #if pad omitted, defalut is false
    #if padval omitted, default is "NA"
    if (missing(pad)) pad=FALSE
    if(missing(padval))padval=NA
    if(missing(na.stuff)) na.stuff=TRUE
    if(missing(func.args)) func.args=NA
    
    len <- length(factorcols)
    colnum <- 0

good<-TRUE
#check that all of the data columns and all of the factors are contained in the DF
temp<-match( datacols, names(df))
if(is.na(sum(temp)))
{
#print(datacols[is.na(temp)])
warning(paste("returned nothing, requested data columns", datacols[is.na(temp)], "do not exist\n", sep=""))
good<-FALSE
}

temp<-match( factorcols, names(df))
if(is.na(sum(temp)))
{
#print(factorcols[is.na(temp)])
warning(paste("returned nothing, requested factor columns", factorcols[is.na(temp)], "do not exist\n",sep=""))
good<-FALSE
}

if(good)
	{
	for(d in c(1:length(datacols)))
          {	


                                        #find the column number whose name matches datacol
            for(i in 1:length(df[1,]))
              { 
                if(names(df)[i]==datacols[d]) {
                  colnum <- i}
              }
            
                                        #only proceed if colnum>0
            if(colnum>0)
              {
                                        #first up create the ID999 field
                
                df$ID999 <-  df[,factorcols[1]]
                if(len>1)
                  {
                    for(i in 2:len)
                      {
                        df$ID999 <- paste(df$ID999, "@", df[,factorcols[i]], sep="")
                      }
                    
                  }
                                        #now run the tapply
                if(is.na(func.args))new.df <- vectorise.names(with(df, tapply(df[,colnum], ID999, func, na.rm=na.stuff)), c(newnames[d], "ID999"))
                if(!is.na(func.args))new.df <- vectorise.names(with(df, tapply(df[,colnum], ID999, func,func.args, na.rm=na.stuff)), c(newnames[d], "ID999"))
                
                if(pad==T)
                  {
                   #now work out all potential permutations of the factors - remove those that are already done and
                   #put padval into those remaining ones.
                    
                                        #get all the possible values for the first factorcols
                    t <- as.character((unique(df[, factorcols[1]])))
                    t <- as.data.frame(t)
                    names(t)[1] <- "ID999"
                    if(length(factorcols)>1)
                      {
                        for(i in c(2:length(factorcols)))
                          {
                                        #how many values for the next factor
                            v <- sort(unique(df[, factorcols[i]]))
                            v.len <- length(v)
                                        #print(v)
                            
                                        #create a new vector long enough to hold the interaction, fill with the existing vector,
                                        #repeated as many times as there are new unique factors
                            t.len <- length(t$ID999)
                            x <-rep(t$ID999,  each=v.len)
                            t2 <- as.data.frame(x)
                            names(t2)[1] <- "ID999"
                                        #print(t)
                                        #print(rep(t,  each=v.len))
                                        #print(v)
                            t2$temp <- rep(v, times=t.len)
                            t2$ID999 <- paste(t2$ID999, t2$temp, sep="@")
                            t <- data.frame.drop(t2, "temp")
                          }
                      }
                    
                    t$temp <- rep(padval, length(t$ID999))
                    done <- unique(new.df$ID999)
                    names(t)[2] <- newnames[d]
                    blank <- t[!t$ID999 %in% done,]
                                        #print(summary(blank))
                    
                                        #now stick onto the bottom of new.df
                    new.df <- rbind(new.df, blank[c(2,1)])
                                        #print(summary(new.df))
                  }
                
                
                
                
                                        #and now unstitch the ID999 field
                s <- strsplit(new.df$ID999, "@")
                sdf <- as.data.frame(unlist(s))
                names(sdf)[1] <- "val"
                sdf$num <- rep(c(1:len), length(sdf[,1])/len)
                                        #work through each of the parts of the ID999 field
                                        #for some reason, integers need special handling to coerce them back into the same format they started in
                for(i in 1:len)
                  {
                    if(class(df[,factorcols[i]])=="integer" || class(df[,factorcols[i]])=="numeric")
                      {
                        new.df[,(i+2)] <-as.character(sdf$val[sdf$num==i])
                      }
                    else if (class(df[,factorcols[i]])=="character" )
                      {
                        new.df[,(i+2)] <-as.character(sdf$val[sdf$num==i])
                      }
                    else
                      {
                        new.df[,(i+2)] <-as.factor(as.character(sdf$val[sdf$num==i]))
                      }
                    
                    
                    class(new.df[,i+2]) <- class(df[,factorcols[i]])
                    names(new.df)[i+2] <- factorcols[i]
                  }
                                        #print(new.df)
              }
            
            if(d==1) {final.df<-new.df}
            if(d>1) 
              {
                final.df<-cbind(final.df, new.df[,1])
                names(final.df)[length(names(final.df))]<-names(new.df)[1]
              }
          }
        
        data.frame.drop(final.df,"ID999")
      }
}
##########################################################################
data.frame.drop<-function(df, cols){
#drops selected columns out of the data frame

    for(i in 1:length(cols)){
	a<-names(df)
	b<-a[a!=cols[i]]
	df<-df[b]
	}
df
}
##########################################################################
NA.to.val<-function( col, val){

pos<-c(1:length(col))
reppos<-pos[is.na(col)]

res<-replace(col, reppos, val)
res
}
##########################################################################
NA.to.0<-function( col){

pos<-c(1:length(col))
reppos<-pos[is.na(col)]

res<-replace(col, reppos, 0)
res
}
##########################################################################
##########################################################################################

vectorise.names<-function (tab,name){

n<-length(attributes(tab)[[1]])
dims<-attributes(tab)[[1]]
len<-prod(attributes(tab)[[1]])
d2<-c(dims, 0)

n1 <- name[1]
n2 <- name[2:length(name)]

#set up the data frame to be the correct length
df<-data.frame(as.vector(tab))
names(df)<-"value"
j<-2

for(i in 1 : n){
ech<- max(1,prod(dims[0:(i-1)]))  # this is the number of sets
reps<-max(1,prod(d2[(i+1):n]))  # this is the number of repeats of each number within a set
df[j]<-rep(dimnames(tab)[[i]],reps,each=ech)
j<-j+1
}
names(df)<-c("value", n2)
names(df)[1] <- n1
df
}


##############################################################

sort.data.frame <- function(form,dat){ 
# Author: Kevin Wright 
# Some ideas from Andy Liaw 
# http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html 


# Use + for ascending, - for decending. 
# Sorting is left to right in the formula 
   

# Useage is either of the following: 
# sort.data.frame(~Block-Variety,Oats) 
# sort.data.frame(Oats,~-Variety+Block) 
   

# If dat is the formula, then switch form and dat 
  if(inherits(dat,"formula")){ 
    f=dat 
    dat=form 
    form=f 
  } 
  if(form[[1]] != "~") 
    stop("Formula must be one-sided.") 

# Make the formula into character and remove spaces 
  formc <- as.character(form[2]) 
  formc <- gsub(" ","",formc) 
# If the first character is not + or -, add + 
  if(!is.element(substring(formc,1,1),c("+","-")))     formc <- paste("+",formc,sep="") 
# Extract the variables from the formula 
# Remove spurious "" terms 
  vars <- unlist(strsplit(formc, "[\\+\\-]"))
   vars <- vars[vars!=""] 


# Build a list of arguments to pass to "order" function 
  calllist <- list() 
  pos=1 # Position of + or - 
  for(i in 1:length(vars)){ 
    varsign <- substring(formc,pos,pos) 
    pos <- pos+1+nchar(vars[i]) 
    if(is.factor(dat[,vars[i]])){ 

      if(varsign=="-")
        calllist[[i]] <- -rank(dat[,vars[i]])
      else
        calllist[[i]] <- rank(dat[,vars[i]])

    } 
    else { 
      if(varsign=="-")
        calllist[[i]] <- -dat[,vars[i]]
      else
        calllist[[i]] <- dat[,vars[i]]


    } 
  } 
  dat[do.call("order",calllist),] 
} 


