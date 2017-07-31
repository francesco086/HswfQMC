#!/bin/bash

HswfQMC_NAME="HswfQMC"

HswfQMC_PATH=$(which pilot-HswfQMC.sh  | sed -e "s/\/pilot-HswfQMC.sh//")
#echo "Path to the executable: "${HswfQMC_PATH}

CURRENT_PATH=$(pwd)
#echo "Current path: "${CURRENT_PATH}

pilot_PATH=$(which pilot-HswfQMC.sh | sed -e "s/\/pilot-HswfQMC.sh//")
#echo "Path to the pilot executable: "${pilot_PATH}

OS_NAME=$(uname)
#echo "The Operating System is: "${OS_NAME}

case ${OS_NAME} in
    "Darwin")
        NUM_CPU=$(sysctl hw.ncpu | sed -e "s/hw.ncpu://" | sed -e "s/ *//" )
        ;;
    "Linux")
        NUM_CPU=$(grep -c ^processor /proc/cpuinfo)
        ;;
esac
#echo "Number of CPU: ${NUM_CPU}, -j${NUM_CPU}"

LAPACK_FOLDER="lapack_lib"

while :
do
    case $1 in
        help)
            echo "These are the options you have:
                                "
            echo "      --- Build and set the HswfQMC code ---"
            echo "set_path - set the PATH shell variable in order to be able to use the pilot script, the HswfQMC_exe executalbe, and the helpers"
            echo "install_markuspline - Download and compile the markuspline library (https://github.com/francesco086/markuspline)"
            echo "install_lapack - Download and compile the lapack library (recommended) "
            echo "set_makefile - Set the Makefile automatically"
            echo "build - Compile (make) the HswfQMC source code in order to have the HswfQMC_exe executable"
            echo "rebuild - Recompile the code from scratch"
            echo "make_ipic - Compile (make) the i-pi client, which connects HswfQMC as a force generator to i-pi."
            echo "git_pull - Update the code to the latest version on GitHub (it might require a rebuild)"

            echo ""
            echo "      --- Use HswfQMC ---"
            echo "set_dir - Make the current folder a working folder (with all necessary input files and folders)"
            echo "generate_orbitals - Compute the DFT orbitals corresponding to the given input file dati_fisici.d and dati_DFT.d using Quantum Espresso (it requires pw.x and iotk from Quantum Espresso)"
            echo "find_SR_minimum - Find the minimum of a SR minimization, by looking at the file ottimizzazione/SR_energies.d"
            echo "clean - Clean all old data from previous simulations"
            echo "wash - Clean all old data from previous simulations but the optimized wf and lattice positions"

            echo ""
            echo "      --- For developers only --- "
            echo "commit - commit and push on github"

            exit
            ;;
        set_dir)
            echo "Make the current folder a working folder (with all necessary input files and folders)"
            mkdir estimatori estimatori/gradiente ottimizzazione ottimizzazione/splines posizioni reticolo
            cp ${pilot_PATH}/input_templates/* .
            cp ${pilot_PATH}/random_seed/randomseed1.d .
            exit
            ;;

        generate_orbitals)
            echo "### Calculate LDA orbitals via Quantum Espresso and place them in the folder orbitals ###"
            #generate the crystal structure with HswfQMC
            FLAG=$(ls | grep -w dati_mc.d)
            if [ "${FLAG}" == "" ]
            then
                echo "ERROR: Impossible to find a file dati_mc.d . ABORT!"
                exit
            fi
            cp dati_mc.d dati_mc.original
            cp dati_funzione_onda.d dati_funzione_onda.original
            cp ${pilot_PATH}/helpers/qespresso/dati_mc.generate_crystal dati_mc.d
            cp ${pilot_PATH}/helpers/qespresso/dati_funzione_onda.generate_crystal dati_funzione_onda.d
            FLAG=$(ls | grep -w posizioni)
            if [ "${FLAG}" == "" ]
            then
                echo "ERROR: Impossible to find the folder posizioni/ . ABORT!"
                exit
            fi
            \rm -f output.d
            HswfQMC_exe > /dev/null 2>&1
            #Extract from the output file the size of the simulation box and save it in the file L.d
            cat output.d | grep "L=" | sed "s/L=   //" | sed "s/   \[bohr\]//" > L.d
            cat output.d | grep "L=" | sed "s/L=   //" | sed "s/   \[bohr\]//"
            #Get the file which contains the crystal structure and put it in the file crystal.d
            mv posizioni/inizio-QEgen_p_0000.pos crystal.d
            #Delete the useless temporary files
            \rm -r -f posizioni/*
            \rm -f output.d
            echo "Generated the crystal structure with HswfQMC"
            #go back to the original dati_mc.d
            mv dati_mc.original dati_mc.d
            mv dati_funzione_onda.original dati_funzione_onda.d
            #fetch the ORBITALS_FOLDER variable from dati_DFT.d
            ORBITALS_FOLDER=$( cat dati_DFT.d | grep -w "ORBITALS_FOLDER" | sed "s/ORBITALS_FOLDER\=//" | sed "s/\"//g" | sed "s/'//g" )
            if [ "${ORBITALS_FOLDER}" == ""  ]
            then
                echo "WARNING: the ORBITALS_FOLDER provided in dati_DFT.d was an empty string, therefore the default has been adopted (qespresso)"
                ORBITALS_FOLDER=qespresso
            fi

            #Build the folders where the orbitals will be created
            mkdir -p orbitals
            cd orbitals/
            #fetch the OVERWRITE variable from dati_DFT.d
            OVERWRITE=$( cat ../dati_DFT.d | grep -w "OVERWRITE" | sed "s/OVERWRITE\=//" )
            #if OVERWRITE is equal to T, then the following check is skipped
            if [ "${OVERWRITE}" != "T" ] 
            then
                #check if there is already a folder with the same name, abort if OVERWRITE=T
                FLAG=$(ls | grep -w ${ORBITALS_FOLDER})
                if [ "${FLAG}" != "" ]
                then
                    echo "ERROR: A folder with the provided name already exist. ABORT!"
                    exit
                fi
            fi
            #create the folder to be used to store the DFT orbitals
            \rm -r -f ${ORBITALS_FOLDER}
            mkdir ${ORBITALS_FOLDER}
            echo "Created the folder orbitals/"${ORBITALS_FOLDER}
            #enter in the folder where the DFT orbitals will be stored
            cd ${ORBITALS_FOLDER}
            #generate the scf.in file, i.e. the input file for Quantum Espresso
            mv ../../L.d .
            mv ../../crystal.d .
            generate_scf.in_from_dati_DFT.d.py
            #run Quantum Espresso
            QE_run
            echo "Quantum Espresso (pw.x) has been executed"
            #convert output .dat files into .xml
            QE_convert_dat_into_xml
            echo "Converted the *.dat files into *.xml (readable for HswfQMC)"
            cd ..
            cd ..
            #fetch the WF variable from dati_DFT.d
            WF=$( cat dati_DFT.d | grep -w "WF" | sed "s/WF\=//" | sed "s/\"//g" | sed "s/'//g" )
            sed -i.sedbak "s/lda_path=.*/lda_path='orbitals\/${ORBITALS_FOLDER}'/" dati_funzione_onda.d
            echo "Set "${WF}" for using the generated orbitals"
            exit
            ;;
        find_SR_minimum)
            find_SR_minimum_energy.py ottimizzazione/SR_energies.dat
            exit
            ;;
        build)
            cd ${pilot_PATH}
            echo "Build the executable file HswfQMC_exe"
            cd source
            make
            mv HswfQMC* ../
            cd $CURRENT_PATH
            exit
            ;;
        rebuild)
            echo "Clean compiled files and compile again?"
            echo "Are you sure? [y/n] "
            read ANSW
            if [ "$ANSW" = "y" ]
            then
                cd $HswfQMC_PATH
                cd source/
                make clean
                make
                mv HswfQMC* ../
                cd $CURRENT_PATH
            else
                echo "Aborted"
            fi
            exit
            ;;
        make_ipic)
            cd ${pilot_PATH}
            echo "Build the i-pi client executable HswfQMC_ipic"
            cd ipi_client
            make clean
            make
            mv ipi_client.x ../HswfQMC_ipic
            cd $CURRENT_PATH
            exit
            ;;
        commit)
            echo "Commit and push on github."
            echo "Are you sure? [y/n] "
            read ANSW
            if [ "$ANSW" = "y" ]
            then
                cd $pilot_PATH
                git add -A
                echo "Provide a comment for this commit"
                read GIT_COMMENT
                git commit -m "$GIT_COMMENT"
                git push -u origin master
                cd $CURRENT_PATH
            else
                echo "Aborted"
            fi
            exit
            ;;
        git_pull)
            echo "Update to the latest version on GitHub by executing 'git pull origin master'."
            echo "Are you sure? [y/n]"
            read ANSW
            if [ "$ANSW" = "y" ]
            then
                git pull origin master
            else
                echo "Aborted"
            fi
            cd
            exit
            ;;
        clean)
            echo "Clean all old datas from previous simulations!"
            echo "Are you sure? [y/n] "
            read ANSW
            if [ "$ANSW" = "y" ]
            then
                rm -f -r estimatori posizioni ottimizzazione
                rm -f output.d
                rm -f nohup.out
                rm -f reticolo/*.d
                mkdir estimatori
                mkdir posizioni
                mkdir estimatori/gradiente
                mkdir ottimizzazione
                mkdir ottimizzazione/splines
            else
                echo "Aborted"
            fi
            exit
            ;;
        wash)
            echo "Clean everything but the optimized wf and lattice positions!"
            echo "Are you sure? [y/n] "
            read ANSW
            if [ "$ANSW" = "y" ]
            then
                rm -f -r estimatori posizioni
                rm -f ottimizzazione/*.dat
                rm -f output.d
                rm -f nohup.out
                rm -f reticolo/SR_Rp-*.d
                mkdir estimatori
                mkdir posizioni
                mkdir estimatori/gradiente
            else
                echo "Aborted"
            fi
            exit
            ;;
        set_path)
            echo "Set the PATH variable in order to be able to use pilot-HswfQMC"
            case ${OS_NAME} in
                "Darwin")
                    echo "-Recognized an OS X system. Setting .bash_profile."
                    FILE_TO_SET="bash_profile"
                    ;;
                "Linux")
                    echo "-Recognized a Linux operating system. Setting .bashrc"
                    FILE_TO_SET="bashrc"
                    ;;
                *)
                    echo "Operative system not recognized."
                    ;;
            esac
            chmod u+x pilot-HswfQMC.sh
            COUNT_SETPATH=$(grep ~/.${FILE_TO_SET} -e "#add path for HswfQMC" | wc -l)
            if [ $COUNT_SETPATH \> "0" ]
            then
                echo "The path was already set, overwriting."
                #delete the previous text
                case ${OS_NAME} in
                    "Darwin")
                        sed -i -e "/#add path for HswfQMC/{N;N;N;N;N;N;N;N;d;}" ~/.${FILE_TO_SET}
                        ;;
                    "Linux")
                        sed -i -e "/#add path for HswfQMC/,+8d" ~/.${FILE_TO_SET}
                        ;;
                esac
            fi
            PATH=${CURRENT_PATH}:${CURRENT_PATH}/helpers/qespresso:${CURRENT_PATH}/helpers/SR:\$PATH
            echo "
#add path for HswfQMC
export PATH=$PATH
#introduce autocompletation feature to pilot-HswfQMC
_pilot-HswfQMC.sh()
{
    local cur=\${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( \$(compgen -W \"git_pull set_path install_markuspline install_lapack set_makefile build rebuild make_ipic set_dir generate_orbitals find_SR_minimum clean wash commit\" -- \$cur) )
}
complete -F _pilot-HswfQMC.sh pilot-HswfQMC.sh" >> ~/.${FILE_TO_SET}
            exit
            ;;
        install_lapack)
            echo "Download and compile the lapack library (For IBM XL please do it manually / use system's libraries!)"
            cd ${pilot_PATH}
            echo "Which fortran compiler do you use? "
            read FF
            if hash svn 2>/dev/null; then
                svn co https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk
            else
                echo "svn command is missing! Please install subversion and run the lapack installation again."
                exit
            fi
            mv trunk ${LAPACK_FOLDER}
            cd ${LAPACK_FOLDER}
            sed -i.bak "s/FORTRAN  = gfortran/FORTRAN = ${FF}/" make.inc.example
            sed -i.bak "s/OPTS     = -O2 -frecursiv/OPTS     = -O3 -frecursiv/" make.inc.example
            sed -i.bak "s/LOADER   = gfortran/LOADER   = ${FF}/" make.inc.example
            \rm make.inc.example.bak
            mv make.inc.example make.inc
            make -j${NUM_CPU} blaslib
            make -j${NUM_CPU} lapacklib
            mv librefblas.a libblas${HswfQMC_NAME}.a
            mv liblapack.a liblapack${HswfQMC_NAME}.a
            cd $CURRENT_PATH
            echo ""
            echo "Lapack library compiled! In order to use it, insert in the Makefile:"
            echo "LIBS=-llapack${HswfQMC_NAME} -lblas${HswfQMC_NAME}"
            echo "LDFLAGS=-L${pilot_PATH}/${LAPACK_FOLDER}"
            exit
            ;;
        install_markuspline)
            echo "Download and compile the markuspline library (https://github.com/francesco086/markuspline)"
            cd ${pilot_PATH}
            echo "Which fortran compiler do you use? [press enter for default: gfortran]"
            read FF
            if [ "${FF}" = "" ]
            then
                FF="gfortran"
                ANSWBG="n"
            else
                echo "Are you using an IBM XL compiler on Blue Gene? [y/n]"
                read ANSWBG
            fi
            \rm -r -f markuspline/
            git clone https://github.com/francesco086/markuspline
            cd markuspline/
            if [ "$ANSWBG" = "y" ]
            then
                patch <../source/markuspline_XL.patch module_markuspline.f90
                ${FF} -c -O3 -qstrict -qarch=qp -qtune=qp -qsimd=auto -qessl -qmaxmem=-1 -qfree=f90 -qport=mod module_markuspline.f90 -L${LAPACK_LIB} -L/bgsys/local/lib -llapack -lesslbg
            else
                ${FF} -c -O3 module_markuspline.f90 -L${pilot_PATH}/${LAPACK_FOLDER} -llapack${HswfQMC_NAME} -lblas${HswfQMC_NAME}
            fi
            ar rcv libmarkuspline.a *.o
            ranlib libmarkuspline.a
            mv libmarkuspline.a libmarkuspline${HswfQMC_NAME}.a
            cd $CURRENT_PATH
            echo ""
            echo "markuspline library compiled! In order to use it, insert in the Makefile:"
            echo "LIBS=-lmarkuspline${HswfQMC_NAME}"
            echo "LDFLAGS=-I${pilot_PATH}/markuspline -L${pilot_PATH}/markuspline"
            exit
            ;;
        set_makefile)
            echo "Set automatically the Makefile in source/"
            COMPUTER=$(hostname)              #hostname
            USERNAME=$(whoami)                #username
            IDENTIFIER="${USERNAME}@${COMPUTER}"
            cd ${pilot_PATH}
            ALREADYSET=$(grep source/makefile.users_settings -e "ifeq (\$(IDENTIFIER),${IDENTIFIER})" | wc -l)
            if [ ${ALREADYSET} \> "0" ]
            then
                echo "Your profile was already saved in 'source/makefile.users_settings'. Do you want to delete it, and set it again? [y/n]"
                read ANSW
                if [ "${ANSW}" = 'y' ]
                then
                    case ${OS_NAME} in
                        "Darwin")
                            sed -i -e "/ifeq (\$(IDENTIFIER),${IDENTIFIER})/{N;N;N;N;N;d;}" source/makefile.users_settings
                            ;;
                        "Linux")
                            sed -i -e "/ifeq (\$(IDENTIFIER),${IDENTIFIER})/,+5d" source/makefile.users_settings
                            ;;
                    esac
                    FLAG=true
                else
                    FLAG=false
                fi
            else
                FLAG=true
            fi
            if [ ${FLAG} == true ]
            then
                echo "Your profile is going to be set."
                echo "ifeq (\$(IDENTIFIER),${IDENTIFIER})" >> source/makefile.users_settings
                echo "        EXEC=\$(EXEC1)" >> source/makefile.users_settings
                echo "Did you use 'pilotHswfQMC install_lapack' and 'pilotHswfQMC install_markuspline'? [y/n]"
                read ANSW
                if [ "${ANSW}" = "y" ]
                then
                    echo "        LIBS=-llapack${HswfQMC_NAME} -lblas${HswfQMC_NAME} -lmarkuspline${HswfQMC_NAME}" >> source/makefile.users_settings
                    echo "        LDFLAGS=-L${pilot_PATH}/${LAPACK_FOLDER} -I${pilot_PATH}/markuspline -L${pilot_PATH}/markuspline" >> source/makefile.users_settings
                else
                    echo "Provide the followings:"
                    echo "LIBS=[press enter for default: -llapack -lblas -lmarkuspline]"
                    read ANSWLIBS
                    if [ "${ANSWLIBS}" = "" ]
                    then
                        echo "        LIBS=-llapack -lblas -lmarkuspline" >> source/makefile.users_settings
                    else
                        echo "        LIBS=${ANSWLIBS}" >> source/makefile.users_settings
                    fi
                    echo "LDFLAGS=[press enter for default: -empty-]"
                    read ANSWLDFLAGS
                    echo "        LDFLAGS=${ANSWLDFLAGS}" >> source/makefile.users_settings
                fi
                echo "Which MPI-FORTRAN compiler are you using?[press enter for default: mpif90]"
                read ANSWFC
                if [ "$ANSWFC" = "" ]
                then
                    echo "        FC=mpif90" >> source/makefile.users_settings
                else
                    echo "        FC=${ANSWFC}" >> source/makefile.users_settings
                    echo "Are you using an IBM XL compiler? [y/n]"
                    read ANSWXL
                    if [ "$ANSWXL" = "y" ]
                    then
                        echo "        XL=true" >> source/makefile.users_settings
                    fi
                fi
                echo "endif" >> source/makefile.users_settings
                echo "" >> source/makefile.users_settings
            fi
            cd ${CURRENT_PATH}
            exit
            ;;
        *)
            echo "The command you typed is not supported. Run \"pilot-HswfQMC.sh help\" to read the supported ones"
            exit
            ;;
    esac
done
