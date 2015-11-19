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

getInitParams <-
function (model) {

	no_of_params=model$kern$nParams
	grid_size=5

	no_of_kernels=length(model$kern$comp)
	paramNames=list()
	for (i in 1:no_of_kernels) {
		paramNames=append(paramNames,model$kern$comp[[i]]$paramNames)
	}

	paramLims=list()

	if ("inverseWidth" %in% paramNames) {
		iw_ind=which(paramNames=="inverseWidth")		
		iw_bounds=model$kern$comp[[iw_ind]]$options$inverseWidthBounds
		paramLims[[iw_ind]]=seq((log(iw_bounds[1])),(log(iw_bounds[2])),length=grid_size)
	}

	var_ind=which(paramNames=="variance")
	for (i in 1:length(var_ind)) {
		paramLims[[var_ind[i]]]=seq(-10,1,length=grid_size)
	}

	gridPoints=expand.grid(paramLims)
	
	lenGrid=nrow(gridPoints)

	if (no_of_params!=ncol(gridPoints)) {
		print("Not all of the parameters are included in the grid search.")
		initial_params=c()
	} else {
		LogLik=matrix(0,lenGrid,1)
		for (i in 1:lenGrid) {
			model1 = gpExpandParam(model, as.vector(as.matrix(gridPoints[i,])))
			LogLik[i] = gpLogLikelihood(model1)
		}	
		ind_maxLogLik=which.max(LogLik)
		initial_params=as.vector(as.matrix(gridPoints[ind_maxLogLik,]))
	}

	return(initial_params)


}
