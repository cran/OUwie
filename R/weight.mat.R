#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters

#written by Jeremy M. Beaulieu

weight.mat<-function(phy, edges, Rate.mat, root.state, assume.station=TRUE){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	k=max(as.numeric(phy$node.label))
	edges=edges
	oldregime=root.state
	oldtime=0
	nodevar=rep(0,max(edges[,3]))
	edges[,4]<-1-edges[,4]
	alpha=Rate.mat[1,]
	
	if(assume.station==TRUE){
		
		W<-matrix(0,ntips,k)
		
		for(j in 1:k){
			
			n.cov=matrix(0, n, n)
			nodecode=matrix(c(ntips+1,0,root.state),1,3)
			#Weight calculated for the root
			#W[,1]<-exp(-alpha[1]*1)
			#Weights for each species per regime
			for(i in 1:length(edges[,1])){
				anc = edges[i, 2]
				desc = edges[i, 3]
				newregime=which(edges[i,5:(k+4)]==1)
				current=edges[i,4]
				
				if(anc%in%nodecode[,1]){
					start=which(nodecode[,1]==anc)
					oldtime=nodecode[start,2]
					oldregime=nodecode[start,3]
				}
				else{
					newrow=c(anc,newtime,oldregime)
					nodecode=rbind(nodecode,newrow)
					oldtime=newtime
				}
				
				if(oldregime==newregime){
					newtime=current
					
					if(newregime==j){
						nodevar[i]=exp(-alpha[oldregime])*(exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime))
									   }
					else{
						nodevar[i]=0
					}
				}
				else{
					newtime=current-((current-oldtime)/2)
					epoch1=exp(-alpha[oldregime])*(exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime))
					oldtime=newtime
					newtime=current
					epoch2=exp(-alpha[newregime])*(exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime))
					if(oldregime==j){
						nodevar[i]=epoch1
					}
					if(newregime==j){
						nodevar[i]=epoch2
					}
				}
				oldregime=newregime
				n.cov[edges[i,2],edges[i,3]]=nodevar[i]
			}
			
			#Generates the mystical S matrix
			S<-matrix(0,n,n)
			for(i in 1:n){
				nn<-Ancestors(phy,i)
				S[nn,i]<-1
			}
			
			#Convert to a sparse matrix
			S<-as.matrix.csr(S)
			#Create a matrix of sums for each node
			
			temp<-n.cov%*%S+n.cov
			#Remove S, nodecode, and n.cov matrices from memory
			rm(nodecode)
			rm(n.cov)
			rm(S)
			
			n.covsums=apply(as.matrix(temp), 2, sum)
			
			rm(temp)
			
			W[1:(ntips),j]<-t(n.covsums[1:ntips])	
			rm(n.covsums)
			
		}
		
	}
	
	if(assume.station==FALSE){
		
		W<-matrix(0,ntips,k+1)
		
		for(j in 1:k){
			
			n.cov=matrix(0, n, n)
			nodecode=matrix(c(ntips+1,0,root.state),1,3)
			#Weight calculated for the root
			W[,1]<-exp(-alpha[1])
			#Weights for each species per regime
			for(i in 1:length(edges[,1])){
				anc = edges[i, 2]
				desc = edges[i, 3]
				newregime=which(edges[i,5:(k+4)]==1)
				current=edges[i,4]
				
				if(anc%in%nodecode[,1]){
					start=which(nodecode[,1]==anc)
					oldtime=nodecode[start,2]
					oldregime=nodecode[start,3]
				}
				else{
					newrow=c(anc,newtime,oldregime)
					nodecode=rbind(nodecode,newrow)
					oldtime=newtime
				}
				
				if(oldregime==newregime){
					newtime=current
					if(newregime==j){
						nodevar[i]=exp(-alpha[root.state])*(exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime))
					}
					else{
						nodevar[i]=0
					}
				}
				else{
					newtime=current-((current-oldtime)/2)
					epoch1=exp(-alpha[root.state])*(exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime))
					oldtime=newtime
					newtime=current
					epoch2=exp(-alpha[root.state])*(exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime))
					if(oldregime==j){
						nodevar[i]=epoch1
					}
					if(newregime==j){
						nodevar[i]=epoch2
					}
				}
				
				oldregime=newregime
				n.cov[edges[i,2],edges[i,3]]=nodevar[i]
			}
			
			#Generates the mystical S matrix
			S<-matrix(0,n,n)
			for(i in 1:n){
				nn<-Ancestors(phy,i)
				S[nn,i]<-1
			}
			
			#Convert to a sparse matrix
			S<-as.matrix.csr(S)
			#Create a matrix of sums for each node
			
			temp<-n.cov%*%S+n.cov
			#Remove S, nodecode, and n.cov matrices from memory
			rm(nodecode)
			rm(n.cov)
			rm(S)
			
			n.covsums=apply(as.matrix(temp), 2, sum)
			
			rm(temp)
			
			W[1:(ntips),j+1]<-t(n.covsums[1:ntips])	
			rm(n.covsums)
			
		}
		
	}

	#Restandardizes W so that the rows sum to 1 -- Generalized. Will reduce to the simpler model if assuming 1 alpha parameter
	W<-W/rowSums(W)
	
	W
}

#Utility function for obtaining mrcas for each species pair
Ancestors<-function (x, node, type = c("all", "parent")) 
{
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    pvector <- numeric(max(x$edge)) # parents
    pvector[child] <- parents    
    type <- match.arg(type)
    if (type == "parent") 
	return(pvector[node])
    res <- numeric(0)
    repeat {
        anc <- pvector[node]
        if (anc == 0) break
        res <- c(res, anc)
        node <- anc
    }
    res
}

