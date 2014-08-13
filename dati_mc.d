&dati_simulazione_mc
N_mc=-10000              !numero di passi MC, se negativo non viene diviso fra processori
N_blank=1000            !numero di passi MC da fare a vuoto (per ogni processore)
N_1ppt=-150                 !con un numero negativo: N_1ppt=-(N_1ppt/100.)*N_part ###(ottimale -150)
flag_TABC=F               !Twist Averaged Boundary Conditions
N_TABC=2000              !numero di passi dopo cui viene eseguito un twist. Se N_TABC<0, ci saranno -N_TABC twist per ogni processo
N_mc_relax_TABC=100
step_e=0.5
step_se=0.1
step_p=0.075
acceptance_rate=50     !in percentuale
flag_continua=F        !continuare da dove si era rimasti con la simulazione precedente?
howtomove='1ppt'        !'allp'=tutti i walkers vengono mossi insieme    '1ppt'=ogni walker viene mosso singolarmente
propmove='gaus'         !'flat'=mossa viene proposta con distribuzione piatta   'gaus'=distribuzione gaussiana e^(-8*(x-x_0)^2/step^2)
trimer_steps=F     !ancora non compatibile con il kernel della forma GD
flag_elettroni=T
flag_protoni=F
flag_shadow=F
flag_E_tot=T           !calcolo l'energia totale?
flag_E_kin=T       !calcolo l'energia cinetica?
flag_E_pot=T       !calcolo l'energia potenziale?
flag_somme_ewald=T
alpha_ewald=-1.0d0
num_k_ewald=515
flag_gr=F             !calcolo le varie g(r)?
N_hist=250
flag_posizioni=F
flag_normalizza_pos=T
N_AV=-6400              !-10000 
flag_MPI=T
what_to_do='simpcal'
stampa_dati_funzione_onda=T
path_dati_funzione_onda='ottimizzazione/OPT_wf.d'           !'dati_funzione_onda.d'      'ottimizzazione/OPT_wf.d'
accuracy_energy_opt=0.001d0     !accuratezza per l'energia quando ottimizzo la funzione d'onda
flag_disk=T       !flag che segnala se i dati vanno salvati sul disco per calcolare la media o no
flag_output=T      !per salvare tutto l'output nel file output.d
quick_error=0		!numero di blocchi per valutare l'errore. Se =0 l'errore viene valutato accuratamente, con la blocking technique
flag_random_file=T          !inizializzare il random generator usando un file specifico?
random_seed_path='random_seed/randomseed1.d'     !file da usare come input per il random seed
/
