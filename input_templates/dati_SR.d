&dati_SR
SR_num_max_WO_MIN=20          !number of maximum SR steps without a new minimum
SR_beta=0.0001                !initial SR_beta
SR_beta_Rp=0.01               !initial SR_beta_Rp (for protonic derivatives)
fSR=T                         !use the fast Stochastic Reconfiguration (fSR) ?
SR_max_SVD_MIN=0.000001       !maximum Singula Value Decomposition minimum value to consider when inverting the SR matrix
SR_change_bound=T             !set some bounds on the change of the variational parameters ?
SR_min_change=0.              !minimum change (in percentage) of the variational parameters in a SR step
SR_max_change=20.             !maximum change (in percentage) of the variational parameters in a SR step
SR_adaptative_beta=F          !use the adaptative beta scheme ? (increase beta when SR does not found a new minimum for too long)
SR_lambda=F                   !use the lambda adaptative scheme ? (explained in reference ...)
max_lambda=2.                 !maximum value for lambda
SR_lambda_Rp=F                !use the lambda adaptative scheme for the protonic coordinates ?
/
