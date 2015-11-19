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

bbgp_test <-
function(x,y,v,indModelCovTypes,depModelCovTypes) {

	model=constructModel(x,y,v,indModelCovTypes)	
	params_init=getInitParams(model)
	model = modelExpandParam(model,params_init)		
	model0 = gpOptimise(model, 0)	
	
        LogLik0=gpLogLikelihood(model0)
	white_var=model0$kern$comp[[1]]$variance

	model=constructModel(x,y,v,depModelCovTypes)	
	params_init=getInitParams(model)
	model = modelExpandParam(model,params_init)		
	model1 = gpOptimise(model, 0)	
	
        LogLik1=gpLogLikelihood(model1)

	BF=LogLik1-LogLik0

	result=list("independentModel"=model0,"dependentModel"=model1,"BF"=BF)

	return(result)

}
