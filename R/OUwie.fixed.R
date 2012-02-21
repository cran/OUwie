#OUwie likelihood calculator

#written by Jeremy M. Beaulieu

#Allows the user to calculate the likelihood given a specified set of parameter values. 

OUwie.fixed<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),root.station=TRUE, alpha=NULL, sigma.sq=NULL, theta=NULL, clade=NULL){
	
	#Makes sure the data is in the same order as the tip labels
	data<-data.frame(data[,2], data[,3], row.names=data[,1])
	data<-data[phy$tip.label,]
	
	#Values to be used throughout
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	#Will label the clade of interest
	if(is.null(clade)){
		phy=phy
	}
	else{
		node<-mrca(phy)[clade[1],clade[2]]
		int<-c(node,Descendants(phy,node, "all"))
		tips<-int[int<ntips]
		data[,1]<-1
		data[tips,1]<-2
		int<-int[int>ntips]
		phy$node.label<-rep(1,phy$Nnode)
		pp<-int-length(phy$tip.label)
		phy$node.label[pp]<-2
	}
	int.states<-factor(phy$node.label)
	phy$node.label=as.numeric(int.states)
	tip.states<-factor(data[,1])
	data[,1]<-as.numeric(tip.states)
	k<-length(levels(int.states))
	
	#A boolean for whether the root theta should be estimated -- default is that it should be.
	root.station=root.station
	
	if (is.character(model)) {
		
		if (model == "BM1"| model == "OU1"){
			##Begins the construction of the edges matrix -- similar to the ouch format##
			#Makes a vector of absolute times in proportion of the total length of the tree
			branch.lengths=rep(0,(n-1))
			branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
			
			#Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
			root.state<-1
			#New tree matrix to be used for subsetting regimes
			edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
			edges=edges[sort.list(edges[,3]),]
			
			edges[,4]=branch.lengths
			regime <- matrix(0,nrow=length(edges[,1]),ncol=2)
			regime[,1]<-1
			regime[,2]<-0
			
			edges=cbind(edges,regime)
		}
		else{
			##Begins the construction of the edges matrix -- similar to the ouch format##
			#Makes a vector of absolute times in proportion of the total length of the tree
			branch.lengths=rep(0,(n-1))
			branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
			
			#Obtain root state and internal node labels
			root.state<-phy$node.label[1]
			int.state<-phy$node.label[-1]
			#New tree matrix to be used for subsetting regimes
			edges=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
			edges=edges[sort.list(edges[,3]),]
			
			edges[,4]=branch.lengths
			
			mm<-c(data[,1],int.state)
			
			regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
			#Generates an indicator matrix from the regime vector
			for (i in 1:length(mm)) {
				regime[i,mm[i]] <- 1 
			}
			#Finishes the edges matrix
			edges=cbind(edges,regime)
			
		}
	}
	#Resort the edge matrix so that it looks like the original matrix order
	edges=edges[sort.list(edges[,1]),]
	x<-as.matrix(data[,2])
	
	#Matches the model with the appropriate parameter matrix structure
	if (is.character(model)) {
		Rate.mat <- matrix(1, 2, k)
		
		if (model == "BM1"){
			np=1
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-0.0000000001
			Rate.mat[2,1:k]<-sigma.sq
			bool=TRUE
		}
		if (model == "BMS"){
			np=k
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-0.0000000001
			Rate.mat[2,1:k]<-sigma.sq
			bool=FALSE
		}
		if (model == "OU1"){
			np=2
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-alpha
			Rate.mat[2,1:k]<-sigma.sq
			bool=root.station
		}
		if (model == "OUM"){
			np=2
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-alpha
			Rate.mat[2,1:k]<-sigma.sq
			bool=root.station
		}
		if (model == "OUMV") {
			np=k+1
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-alpha
			Rate.mat[2,1:k]<-sigma.sq
			bool=root.station
		}
		if (model == "OUMA") {
			np=k+1
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-alpha
			Rate.mat[2,1:k]<-sigma.sq
			bool=root.station
		}
		if (model == "OUMVA") {
			np=k*2
			index<-matrix(TRUE,2,k)
			Rate.mat[1,1:k]<-alpha
			Rate.mat[2,1:k]<-sigma.sq
			bool=root.station
		}
	}
	obj<-NULL
	
	#Likelihood function for estimating model parameters
	dev<-function(){
		
		N<-length(x[,1])
		V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state)
		W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, assume.station=bool)
		
		if(is.null(theta)){
			theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
			se<-sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
		}
		else{
			theta=theta
			se=rep(NA,length(theta))
		}

		theta.est<-cbind(theta,se)
		
		DET<-determinant(V, logarithm=TRUE)
		
		res<-W%*%theta-x		
		q<-t(res)%*%solve(V,res)
		logl <- -.5*(N*log(2*pi)+as.numeric(DET$modulus)+q[1,1])
		
		list(-logl,theta.est)
	}
	#Informs the user that the optimization routine has started and starting value is being used (default=1)
	cat("Calculating likelihood using fixed parameter values:",c(alpha,sigma.sq,theta), "\n")
	
	loglik <- dev()
	obj$loglik<- -loglik[[1]]
	
	if (is.character(model)) {
		if (model == "BM1"){
			obj$AIC <- -2*obj$loglik+2*(np+1)
			obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
			obj$Param.est <- Rate.mat
			obj$Param.est[1,1:k] <- NA
			rownames(obj$Param.est)<-c("alpha","sigma.sq")
			colnames(obj$Param.est) <- levels(int.states)
			theta <- loglik[[2]]
			obj$ahat <- matrix(theta[1,], 1,2)
			colnames(obj$ahat) <- c("Estimate", "SE")
		}
		if (model == "BMS"){
			obj$AIC <- -2*obj$loglik+2*(np+1)
			obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
			obj$Param.est <- Rate.mat
			obj$Param.est[1,1:k] <- NA
			rownames(obj$Param.est)<-c("alpha","sigma.sq")
			colnames(obj$Param.est) <- levels(int.states)
			theta <- loglik[[2]]
			obj$ahat <- matrix(theta[1,], 1,2)
			colnames(obj$ahat) <- c("Estimate", "SE")
		}
		if (root.station == TRUE){
			if (model == "OU1"){
				obj$AIC <- -2*obj$loglik+2*(np+1)
				obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
				obj$Param.est<- Rate.mat
				rownames(obj$Param.est)<-c("alpha","sigma.sq")
				colnames(obj$Param.est)<-levels(int.states)
				theta <- loglik[[2]]
				obj$theta <- matrix(theta[1,], 1,2)
				colnames(obj$theta) <- c("Estimate", "SE")			}
		}
		if (root.station == FALSE){
			if (model == "OU1"){
				obj$AIC <- -2*obj$loglik+2*(np+1)
				obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
				obj$Param.est<- Rate.mat
				rownames(obj$Param.est)<-c("alpha","sigma.sq")
				colnames(obj$Param.est)<-levels(int.states)
				theta <- loglik[[2]]
				obj$theta<-theta[1:2,1:2]
				rownames(obj$theta)<-c("Root", "Primary")
				colnames(obj$theta)<-c("Estimate", "SE")
			}
		}
		if (root.station == FALSE){
			if (model == "OUM"| model == "OUMV"| model == "OUMA" | model == "OUMVA"){ 
				obj$AIC <- -2*obj$loglik+2*(np+k)
				obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
				obj$Param.est <- Rate.mat
				rownames(obj$Param.est)<-c("alpha","sigma.sq")
				colnames(obj$Param.est)<-levels(int.states)
				obj$theta<-loglik[[2]]
				rownames(obj$theta)<-c("Root",levels(int.states))
				colnames(obj$theta)<-c("Estimate", "SE")
			}
		}
		if (root.station == TRUE){
			if (model == "OUM"| model == "OUMV"| model == "OUMA" | model == "OUMVA"){ 
				obj$AIC <- -2*obj$loglik+2*(np+k)
				obj$AICc <- -2*obj$loglik+(2*(np+k)*(ntips/(ntips-(np+k)-1)))
				obj$Param.est<- Rate.mat
				rownames(obj$Param.est)<-c("alpha","sigma.sq")
				colnames(obj$Param.est)<-levels(int.states)				
				obj$theta<-loglik[[2]]
				rownames(obj$theta)<-c(levels(int.states))
				colnames(obj$theta)<-c("Estimate", "SE")
			}
		}
	}
	
	obj
}
