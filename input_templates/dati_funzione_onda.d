&dati_funzione_onda
SDe_kind='bat'        !'pw_'=plane wave   'lda'=orbitali da LDA   'prf'=proton field (hartree)  'fre'=free con hartree  'atm/atp'=exp(-r*C_atm)   'bat'=bi-atomic   'no_'=no SD
Jee_kind='yup'        !'yuk'=Yukawa  'yup'=Yukawa con Periodic Coordinate     'no_'=no Jee
Jep_kind='yup'        !'yuk'=Yukawa  'yup'=Yukawa con Periodic Coordinate  'atm'=exp(-F*r)   'atp'=atm with PC  'no_'=no Jep
Jpp_kind='no_'        !'no_'=no Jpp
SDse_kind='no_'       !'pw_'=plane wave   'pw2'=plane wave squared    'lda'=orbitali da LDA   'no_'=no SD    'gem'=geminal of gaussians    'gss'=gaussians centered on the protons    'atm'=exp(-r)   'atp'=atm con PC
Jse_kind='no_'        !'pot'=potenziale effettivo riscalato (richiede flag_traccia_coppie in dati_fisici)   'bou'=bounding B*(r-D)^2  'ppb'=potenziale riscalato piú bounding 'yuk'=Yukawa  'no_'=no Jse
Kse_kind='no_'        !'gsd'=gaussian determinant   'gdc'=gsd con ctf  'gss'=gaussiana   'gsc'=gaussiana con ctf  'no_'=no Kernel   'gsp'=gaussiana con PC    'atm'=exp    'atc'=exp con ctf
Jsesp_kind='no_'      !'pot'=potenziale riscalato    'yuk'=Yukawa   'gss'=gaussian    'gsd'=gaussian determinant   'no_'=no Jsesp
split_Aee=T
split_Aep=T
split_Asese=T
split_Asesp=T
split_Fee=T
split_Fep=T
split_Fsese=T
split_Fsesp=T
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
GSWF=1.
C_ATM=1.2d0
N_ritraccia_coppie=1000           !con un numero <0 non vengono ricalcolate 
N_mc_relax_traccia_coppie=10
A_POT_se=1.                      !per il jastrow pot se-se
D_POT_se=1.
Gsesp=1.                         !per il Jsesp gaussiano
c_se=1.    !1.1d0                !per il jastrow fra le shadow-e, Jse
B_se=-1.4d0
D_se=-0.7d0                          !per il bounding nel Jse
C_sesp=1.                      !per il jastrow per shadow-e e shadow-p
lda_path='orbitals'       !'../lda_orbitals/r_s=1.31/BCC-54/ctf=10/'       !path per trovare i file necessari per usare orbitali LDA o Hartree. 'genera_on_the_fly' li fa generare
kf_coeff_hartree=1.d0   !se =0.d0 si usa il file vecchio
flag_usa_coeff_hartree=.FALSE.
/

