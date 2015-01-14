&dati_funzione_onda
SDe_kind='bat'        !Electronic Slater Determinant: 'pw_'=simple plane waves, 'lda'=DFT orbitals (to be specifies), 'prf'=diagonalized nuclear field Hamiltonian (dnfH), 'fre'=dnfH without interactions, 'atm'/'atp'=exp(-r*C_atm) with or without PC, 'bat'=bi-atomic, '1sb'=dynamic e-p backflow with 1s orbitals, 'spb'=dynamic e-p backflow with spline orbital, 'no_'=no SD
Jee_kind='spp'        !Electron-Electron Jastrow: 'yuk'/'yup'=Yukawa with or without PC, 'spl'/'spp'=spline with or without PC, inizialized to fit the Yukawa Jastrow, 'no_'=no Jee
Jep_kind='spp'        !Electronic-Proton Jastrow: 'yuk'/'yup'=Yukawa, 'spl'/'spp'=spline with or without PC, 'atm'/'atp'=exp(-F*r), 'no_'=no Jep
Jpp_kind='no_'        !Proton-Proton Jastrow: 'no_'=no Jpp
SDse_kind='no_'       !eShadow Slater Determinant: 'pw_'=simple plane waves, 'pw2'=plane waves squared(for sign problem), 'lda'=DFT orbitals, 'no_'=no SD, 'gem'=geminal of gaussians, 'gss'=gaussians centered on the protons, 'atm'/'atp'=exp(-r)
Jse_kind='no_'        !eShadow-eShadow Jastro: 'yuk'/'yup'=Yukawa, 'pot'=effective potential (DEPRECATED), 'no_'=no Jse
Kse_kind='no_'        !eShadow Kernel: 'gss'/'gsp'=gaussian, 'gsc'=gaussian with cut-off 'gsd'/'gdp'=gaussian determinant, 'gdc'=gsd with cut-off, 'no_'=no Kernel, 'atm'/'atp'=exp(-C_ATM*(r-se)), 'atc'=atm with ctf
Jsesp_kind='no_'      !eShadow-pShadow (at the moment =Proton) Jastrow: 'yuk'/'yup'=Yukawa, 'gss'=gaussian, 'gsd'=gaussian determinant (for sign problem), 'no_'=no Jsesp
split_Aee=F           !Spin split for the Aee Yukawa parameters
split_Aep=F           !Spin split for Aep
split_Asese=F         !Spin split for Asese
split_Asesp=F         !Spin split for Asesp
split_Fee=F           !Spin split for Fee
split_Fep=F           !Spin split for Fep
split_Fsese=F         !Spin split for Fsese
split_Fsesp=F         !Spin split for Fsesp
m_Bsplep=1
nknots_Bsplep=5
cutoff_Bsplep=F
m_Jsplee=1
nknots_Jsplee=5
cutoff_Jsplee=F
m_Jsplep=1
nknots_Jsplep=5
cutoff_Jsplep=F
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
C_ATM=1.                          !For SDse='atm'
N_ritraccia_coppie=1000           !(DEPRECATED) if < 0 the molecular pairs are not tracked
N_mc_relax_traccia_coppie=10      !(DEPRECATED)
A_POT_se=3.                       !For SDe_kind='1sb'
D_POT_se=1.                      
Gsesp=1.                          !For Jsesp='gss'/'gsd'
c_se=1.                           !(DEPRECATED) for the Jse
B_se=-1.4                         !(DEPRECATED)
D_se=-0.7                         !(DEPRECATED) For the bounding in the Jse
C_sesp=1.                         !(DEPRECATED)
lda_path='orbitals'               !Path to the file (or files, with TABC) that specifies the DFT orbitals generated with Quantum Espresso
kf_coeff_dnfH=1.                  !(DEPRECATED)
flag_usa_coeff_dnfH=.FALSE.       !(DEPRECATED)
/

