#!/bin/bash

HswfQMC_NAME="HswfQMC"

HswfQMC_PATH=$(which HswfQMC_exe  | sed -e "s/\/HswfQMC_exe//")
#echo "Path to the executable: "${HswfQMC_PATH}

CURRENT_PATH=$(pwd)
#echo "Current path: "${CURRENT_PATH}

pilot_PATH=$(which pilot-HswfQMC.sh | sed -e "s/\/pilot-HswfQMC.sh//")
#echo "Path to the pilot executable: "${pilot_PATH}

OS_NAME=$(uname)
#echo "The Operating System is : "${OS_NAME}

LAPACK_FOLDER="lapack_lib"

while :
do
	case $1 in
		help)
			echo "These are the options you have:
								"
			echo "      --- Build and set the HswfQMC code ---"
			echo "setpath - set the PATH shell variable in order to be able to use the pilot script and the HswfQMC_exe executalbe"
			echo "install_lapack - Download and compile the lapack library (recommended) "
			echo "setmakefile - Set the Makefile automatically"
			echo "build - Compile (make) the HswfQMC source code in order to have the HswfQMC_exe executable"
			echo "recompile - Recompile the code from scratch"
			
			echo ""
			echo "      --- Use HswfQMC ---"
			echo "setdir - Make the current folder a working folder (with all necessary input files and folders)"
			echo "clean - Clean all old datas from previous simulations"
			echo "wash - Clean all old datas from previous simulations but the optimized wf and lattice positions"
			
			echo ""
			echo "      --- For developers only --- "
			echo "commit - commit and push on github"
			
			exit
			;;
		setdir)
			echo "Make the current folder a working folder (with all necessary input files and folders)"
			mkdir estimatori ottimizzazione ottimizzazione/gradiente posizioni reticolo
			cp ${pilot_PATH}/input_templates/* .
			PATHRANDOM="${pilot_PATH}/random_seed"
			sed -i.sedbak "s|RANDOM_SEED_FOLDER|${PATHRANDOM}|" dati_mc.d
			rm *.sedbak
			exit
			;;	
		build)
			cd ${pilot_PATH}
			echo "Build the executable file and necessary folders"
			cd source
			make
			mv HswfQMC* ../
			cd $CURRENT_PATH
			exit
			;;
		recompile)
			echo "Clean compiled files and compile again?"
			echo "Are you sure? [y/n] "
			read ANSW
			if [ "$ANSW" = "y" ]
			then
				cd $HswfQMC_PATH
				\rm HswfQMC*
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
		clean)
			echo "Clean all old datas from previous simulations!"
			echo "Are you sure? [y/n] "
			read ANSW
			if [ "$ANSW" = "y" ]
			then
				rm -v -f -r estimatori posizioni ottimizzazione
				rm -v -f output.d
				rm -v -f nohup.out
				rm -v -f reticolo/*.d
				mkdir -v estimatori
				mkdir -v posizioni
				mkdir -v estimatori/gradiente
				mkdir -v ottimizzazione
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
				rm -v -f -r estimatori posizioni
				rm -v -f ottimizzazione/*.dat
				rm -v -f output.d
				rm -v -f nohup.out
				rm -v -f reticolo/SR_Rp-*.d
				mkdir -v estimatori
				mkdir -v posizioni
				mkdir -v estimatori/gradiente
			else
				echo "Aborted"
			fi
			exit
			;;
		setpath)
			echo "Set the PATH variable in order to be able to use pilot-HswfQMC"
			if [ "${OS_NAME}" = "Darwin" ]
			then
				echo "-Recognized an OS X system. Modify .bash_profile instead of .bashrc"
				FILE_TO_SET="bash_profile"
			else
				echo "-Standard platform"
				FILE_TO_SET="bashrc"
			fi
			chmod u+x pilot-HswfQMC.sh
			COUNT_SETPATH=$(grep ~/.${FILE_TO_SET} -e "PATH=${CURRENT_PATH}:\$PATH" | wc -l)
			if [ $COUNT_SETPATH == 0 ]
			then
				PATH=${CURRENT_PATH}:\$PATH
				echo "" >> ~/.${FILE_TO_SET}
				echo "#add path for HswfQMC" >> ~/.${FILE_TO_SET} 
				echo "PATH=${CURRENT_PATH}:\$PATH" >> ~/.${FILE_TO_SET}
			else
				echo "-The variable PATH was already correctly set"
			fi
			exit
			;;
		install_lapack)
			echo "Download and compile the lapack library"
			cd ${pilot_PATH}
			echo "Which fortran compiler do you use? "
                        read FF
			echo "How many cores does your computer have? [if you are not sure type 1] "
                        read NUM_CPU
			svn co https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk
			mv trunk ${LAPACK_FOLDER}
			cd ${LAPACK_FOLDER}
			sed -i.bak "s/FORTRAN  = gfortran/FORTRAN = ${FF}/" make.inc.example
			sed -i.bak "s/OPTS     = -O2 -frecursiv/OPTS     = -O0 -march=native -frecursiv/" make.inc.example
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
		setmakefile)
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
					sed -i -e "/ifeq (\$(IDENTIFIER),${IDENTIFIER})/,+5d" source/makefile.users_settings   #might not work for Os X
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
				echo "Did you use pilotHswfQMC install_lapack? [y/n]"
				read ANSW
				if [ "${ANSW}" = "y" ]
				then
					echo "        LIBS=-llapack${HswfQMC_NAME} -lblas${HswfQMC_NAME}" >> source/makefile.users_settings
					echo "        LDFLAGS=-L${pilot_PATH}/${LAPACK_FOLDER}" >> source/makefile.users_settings
				else
					echo "Provide the followings:"
					echo "LIBS=[press enter for default: -llapack -lblas]"
					read ANSWLIBS
					if [ "${ANSWLIBS}" = "" ]
					then
						echo "        LIBS=-llapack -lblas" >> source/makefile.users_settings
					else
						echo "        LIBS=${ANSWLIBS}" >> source/makefile.users_settings
					fi
					echo echo "LDFLAGS=[press enter for default: -empty-]"
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
