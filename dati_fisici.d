&dati_fisici
r_s=1.31d0                  !solido 1.31, liquido2 2.61, liquidH 1.251, per avere L=1 con N=16 0.246186255460674
crystal_cell='bcc__'          !sono permessi 'bcc' o 'fcc' o 'hcp' o 'hcp_w' o 'mhcpo' o 'sc_' o 'mol' o 'dat' o 'datex' o 'grp__'
file_reticolo='reticolo/guglielmo.pos'  !file che specifica le posizioni degli atomi
flag_molecular=F            !implementata solo per hcp e mhcpo
strecthing_cov_bond=1.
N_cell_side=2           !se crystal_cell='dat' allora bisogna scrivere qua il numero di particelle
/

! Le distanze sono espresse in unitá di a_0, il raggio di Bohr
! quindi per esprimerle in Angstrom bisogna moltiplicare per 0.53
! Le energie sono espresse in Rydberg per atomo (per averle in Hartree bisogna dividere per 2)
!LiquidH2 T=300K  P~1GPa
!LiquidH T=1000K
!
!'datex' sono dati con posizioni non normalizzate, quindi devono avere un L giusto (la prima riga riporta L)
!'dat' sono posizioni normalizzate
