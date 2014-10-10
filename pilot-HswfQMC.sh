#!/bin/bash

HswfQMC_PATH=$(which HswfQMC_exe  | sed -e "s/\/HswfQMC_exe//")
#echo "Path to the executable: "${HswfQMC_PATH}
CURRENT_PATH=$(pwd)
#echo "Current path: "${CURRENT_PATH}
pilot_PATH=$(which pilot-HswfQMC.sh | sed -e "s/\/pilot-HswfQMC.sh//")
#echo "Path to the pilot executable: "${pilot_PATH}

while :
do
	case $1 in
		help)
			echo "These are the options you have:
								"
			echo "install_lapack - Download and compile the lapack library (recommended) "
			echo "build - Compile the source code in order to have the HswfQMC_exe executable"
			echo "recompile - Recompile the code from scratch"
			
			echo "setdir - Make the current folder a working folder (with all necessary input files and folders)"
			echo "clean - Clean all old datas from previous simulations"
			echo "wash - Clean all old datas from previous simulations but the optimized wf and lattice positions"

			echo "commit - commit and push on github"
			
			exit
			;;
		setdir)
			if [ "$HswfQMC_PATH" == "" ]
			then
				echo "You first have to compile the code. Run \"pilot-HswfQMC.sh build\""
			else
				echo "Make the current folder a working folder (with all necessary input files and folders)"
				mkdir estimatori ottimizzazione ottimizzazione/gradiente posizioni reticolo
				cp ${HswfQMC_PATH}/input_templates/* .
				PATHRANDOM="${HswfQMC_PATH}/random_seed"
				sed -i.sedbak "s|RANDOM_SEED_FOLDER|${PATHRANDOM}|" dati_mc.d
				rm *.sedbak
			fi
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
			chmod u+x pilot-HswfQMC.sh
			echo "" >> ~/.bashrc
			echo "#add path for HswfQMC" >> ~/.bashrc 
			echo "PATH=${CURRENT_PATH}:\$PATH" >> ~/.bashrc
			exit
			;;
		install_lapack)
			echo "Download and compile the lapack library"
			LAPACK_FOLDER="lapack_lib"
			cd ${pilot_PATH}
			echo "Which fortran compiler do you use? "
                        read FF
			echo "How many cores does your computer have? [if you are not sure type 1] "
                        read NUM_CPU
			svn co https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk
			mv trunk ${LAPACK_FOLDER}
			cd ${LAPACK_FOLDER}
			sed -i.bak "s/FORTRAN  = gfortran/FORTRAN = ${FF}/" make.inc.example
			sed -i.bak "s/OPTS     = -O2 -frecursiv/OPTS     = -O3 -march=native -frecursiv/" make.inc.example
			\rm make.inc.example.bak
			mv make.inc.example make.inc
			make -j${NUM_CPU} blaslib
			mv librefblas.a librefmyblas.a
			make -j${NUM_CPU} lapacklib
			mv liblapack.a libmylapack.a
			cd $CURRENT_PATH
			echo ""
			echo "Lapack library compiled! In order to use it, insert in the Makefile:"
			echo "LIBS=-lmylapack -lmyblas"
			echo "LDFLAGS=-L${pilot_PATH}/${LAPACK_FOLDER}"
			exit
			;;
		*)
			echo "The command you typed is not supported. Run \"pilot-HswfQMC.sh help\" to read the supported ones"
			exit
			;;
	esac
done
