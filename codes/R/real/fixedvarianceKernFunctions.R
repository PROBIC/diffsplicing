# Copyright (c) 2014, Hande TOPA
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

fixedvarianceKernCompute <-
function (kern, x, x2) {
	if ( nargs()<3 ) {
		k <- diag(kern$fixedvariance[,1])
	} else {
		x1dim <- dim(as.array(x))[1]
		x2dim <- dim(as.array(x2))[1]
		k <- matrix(0, nrow=x1dim, ncol=x2dim)
	}
	return (k)
}

fixedvarianceKernExtractParam <-
function (kern, only.values=TRUE,
untransformed.values=TRUE) {
	params <- c()
	
	return (params)
}


fixedvarianceKernDiagCompute <-
function (kern, x) {
	k <- kern$fixedvariance
	return (k)
}


fixedvarianceKernGradient <-
function(kern,x,x2,covGrad) {
	
	g=0

	return(g)
	
}

fixedvarianceKernExpandParam <-
function(kern,params) {
	kern=kern
	return(kern)
	
}

fixedvarianceKernParamInit <-
function(kern) {
	kern$fixedvariance=kern$options$variance
	kern$fixedvariance_times=kern$options$input
	kern$use_fixedvariance=1
	kern$nParams=0
	kern$isStationary=TRUE
	
	return(kern)
	
}
