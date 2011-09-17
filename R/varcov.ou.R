#OU variance-covariance matrix generator

#written by Jeremy M. Beaulieu

varcov.ou<-function(phy, edges, Rate.mat, root.state){

	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	k=max(as.numeric(phy$node.label))
	edges[,4]<-1-edges[,4]
	oldregime=root.state
	oldtime=0
	nodevar1=rep(0,max(edges[,3]))
	nodevar2=rep(0,max(edges[,3]))
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]
	n.cov1=matrix(rep(0,n*n), n, n)
	n.cov2=matrix(rep(0,n*n), n, n)
	nodecode=matrix(c(ntips+1,0,root.state),1,3)

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
			nodevar1[i]=alpha[oldregime]*(newtime-oldtime)
			nodevar2[i]=sigma[oldregime]*exp(2*alpha[oldregime]*(newtime+oldtime)/2)*(newtime-oldtime)
		}
		else{
			newtime=current-((current-oldtime)/2)
			epoch1a=(alpha[oldregime])*(newtime-oldtime)
			epoch1b=sigma[oldregime]*exp(2*alpha[oldregime]*(newtime+oldtime)/2)*(newtime-oldtime)
			oldtime=newtime
			newtime=current
			epoch2a=alpha[newregime]*(newtime-oldtime)
			epoch2b=sigma[newregime]*exp(2*alpha[newregime]*(newtime+oldtime)/2)*(newtime-oldtime)
			nodevar1[i]<-epoch1a+epoch2a
			nodevar2[i]<-epoch1b+epoch2b
		}
		oldregime=newregime
		n.cov1[edges[i,2],edges[i,3]]=nodevar1[i]
		n.cov2[edges[i,2],edges[i,3]]=nodevar2[i]
	}
		
	#Remove nodecode matrix from memory
	rm(nodecode)
	vcv1<-mat.gen(n.cov1,phy)
	vcv2<-mat.gen(n.cov2,phy)
	vcv<-exp(-2*diag(vcv1))*vcv2
	rm(vcv1)
	rm(vcv2)

	vcv
	
}

#Utility for building a summary matrix -- slow for now, need to work on speeding up

mat.gen<-function(mat,phy){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	
	A=matrix(rep(0,n*n), n, n)
	A[phy$edge]=1
	mm<-c(1:n)
	Nt<-t(t(A)*mm)
	#Convert to a sparse matrix
	Nt<-as.matrix.csr(Nt)
	rm(A)
	
	#Generates the mystical S matrix
	S<-matrix(0,n,n)
	for(i in 1:n){
		nn<-Ancestors(phy,i)
		S[nn,i]<-1
	}
	#Convert to a sparse matrix
	S<-as.matrix.csr(S)
	
	#Create a matrix of sums for each node
	temp<-mat%*%S+mat
	n.covsums=apply(as.matrix(temp), 2, sum)
	rm(temp)
	rm(mat)
	#Generates a matrix that lists the descendants for each ancestral node
	H.temp=Nt%*%S+Nt
	rm(Nt)
	rm(S)
	H.mat=as.matrix(H.temp[(ntips+1):n,])
	rm(H.temp)
	H.mat
	
	new.mat<-matrix(0,ntips,ntips)
	diag(new.mat)=n.covsums[1:ntips]
	#Enters covariances
	for(i in 1:length(H.mat[,1])){
		temp=unique(H.mat[i,])
		temp=temp[temp!=0]
		tempR=which(H.mat[i,]==temp[1])
		tempR=subset(tempR,tempR%in%c(1:ntips))
		tempL=which(H.mat[i,]==temp[2])
		tempL=subset(tempL,tempL%in%c(1:ntips))
		
		new.mat[tempL,tempR]=n.covsums[i+ntips]	
		new.mat[tempR,tempL]=n.covsums[i+ntips]	
	}
	
	rm(H.mat)
	
	new.mat
	
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

