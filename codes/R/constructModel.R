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

constructModel <-
function (x,y,v,covarianceTypes) {

        list_covTypes=list()

        if ("rbf" %in% covarianceTypes) {
		l_bound = calculateLbound(x)
		iw_upperbound=1/(l_bound^2)
		iw_lowerbound=1/(tail(x,1)^2)
		list_rbf=list(type="rbf",options=list(inverseWidthBounds=c(iw_lowerbound,iw_upperbound)))
		list_covTypes=append(list_covTypes,list(list_rbf))
	}
        if ("bias" %in% covarianceTypes) {
                list_bias=list(type="bias")
                list_covTypes=append(list_covTypes,list(list_bias))
        }
        if ("white" %in% covarianceTypes) {
		list_white=list(type="white")
		list_covTypes=append(list_covTypes,list(list_white))
	}
        if ("fixedvariance" %in% covarianceTypes) {
		list_fixedvar=list(type="parametric", realType="fixedvariance",options=list(variance=v,input=x))
		list_covTypes=append(list_covTypes,list(list_fixedvar))
	}


	options=gpOptions(approx="ftc")
	options$kern=list(type="cmpnd",comp=list_covTypes)
	model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)

	return(model)

}




