&dati_funzione_onda
SDe_kind='bat'        !Electronic Slater Determinant: 'pw_'=simple plane waves, 'lda'=DFT orbitals (to be specifies), 'prf'=diagonalized nuclear field Hamiltonian (dnfH), 'fre'=dnfH without interactions, 'atm'/'atp'=exp(-r*C_atm) with or without PC, 'bat'=bi-atomic, 'no_'=no SD
Jee_kind='yup'        !Electron-Electron Jastrow: 'yuk'/'yup'=Yukawa with or without PC, 'no_'=no Jee
Jep_kind='yup'        !Electronic-Proton Jastrow: 'yuk'/'yup'=Yukawa, 'atm'/'atp'=exp(-F*r), 'no_'=no Jep
Jpp_kind='no_'        !Proton-Proton Jastrow: 'no_'=no Jpp
SDse_kind='no_'       !eShadow Slater Determinant: 'pw_'=simple plane waves, 'pw2'=plane waves squared(for sign problem), 'lda'=DFT orbitals, 'no_'=no SD, 'gem'=geminal of gaussians, 'gss'=gaussians centered on the protons, 'atm'/'atp'=exp(-r)
Jse_kind='no_'        !eShadow-eShadow Jastro: 'yuk'/'yup'=Yukawa, 'pot'=effective potential (DEPRECATED), 'no_'=no Jse
Kse_kind='no_'        !eShadow Kernel: 'gss'/'gsp'=gaussian, 'gsc'=gaussian with cut-off 'gsd'/'gdp'=gaussian determinant, 'gdc'=gsd with cut-off, 'no_'=no Kernel, 'atm'/'atp'=exp(-C_ATM*(r-se)), 'atc'=atm with ctf
Jsesp_kind='no_'      !eShadow-pShadow (at the moment =Proton) Jastrow: 'yuk'/'yup'=Yukawa, 'gss'=gaussian, 'gsd'=gaussian determinant (for sign problem), 'no_'=no Jsesp
split_Aee=T           !Spin split for the Aee Yukawa parameters
split_Aep=T           !Spin split for Aep
split_Asese=T         !Spin split for Asese
split_Asesp=T         !Spin split for Asesp
split_Fee=T           !Spin split for Fee
split_Fep=T           !Spin split for Fep
split_Fsese=T         !Spin split for Fsese
split_Fsesp=T         !Spin split for Fsesp
AEE_YUK=  1.2819801758793221     ,
AEE_UD_YUK=  2.2098611377210831     ,
FEE_YUK=  1.1178362429152959     ,
FEE_UD_YUK=  1.0559440836004022     ,
AEP_YUK=-0.45378310533588501     ,
AEP_UD_YUK=-0.45833583582447029     ,
FEP_YUK=  3.8627648124100249     ,
FEP_UD_YUK=  3.8228039597443910     ,
C_KERN_E= 0.4     ,
ASESE_YUK=  1.2819801758793221    ,
ASESE_UD_YUK= 2.2098611377210831     ,
FSESE_YUK=  1.1178362429152959     ,
FSESE_UD_YUK=  1.0559440836004022     ,
ASESP_YUK= -0.45378310533588501     ,
ASESP_UD_YUK= -0.45833583582447029     ,
FSESP_YUK=  3.8627648124100249     ,
FSESP_UD_YUK= 3.8228039597443910     ,
GSWF=1.                           !For SDse='gss'/'gsd'
C_ATM=1.2d0                       !For SDse='atm'
N_ritraccia_coppie=1000           !(DEPRECATED) if < 0 the molecular pairs are not tracked
N_mc_relax_traccia_coppie=10      !(DEPRECATED)
A_POT_se=1.                       !(DEPRECATED) For Jse='pot'
D_POT_se=1.                       !(DEPRECATED)
Gsesp=1.                          !For Jsesp='gss'/'gsd'
c_se=1.    !1.1d0                 !(DEPRECATED) for the Jse
B_se=-1.4d0                       !(DEPRECATED)
D_se=-0.7d0                       !(DEPRECATED) For the bounding in the Jse
C_sesp=1.                         !(DEPRECATED)
lda_path='orbitals'               !Path to the file (or files, with TABC) that specifies the DFT orbitals generated with Quantum Espresso
kf_coeff_dnfH=1.d0             !(DEPRECATED)
flag_usa_coeff_dnfH=.FALSE.    !(DEPRECATED)
/

