# Copyright (c) 2015, Hande TOPA
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

negloglik <-
function(theta,y,mu) {

	# theta = [alpha; beta]
	# y = the data vector, M*R matrix for M transcripts R replicates
 
	M=dim(y)[1]
	R=dim(y)[2]
	alpha=theta[1]
	beta=theta[2]

	if (alpha<0 | beta<0) {
   		negloglik=-log(0)
	} else {
  		y_2=y^2
  		A=as.matrix(rowSums(y_2))
  		B=2*mu*as.matrix(rowSums(y))-R*(mu^2)
  		C=A-B
  		negloglik=-(M*alpha*log(beta)-M*lgamma(alpha)+M*lgamma(alpha+R/2)-(alpha+R/2)*sum(log(beta+0.5*C)))
	}
	
	return(negloglik)
	
}
