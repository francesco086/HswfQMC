cp scf.in OUT.save/K00001/
cp scf.out OUT.save/K00001/
cd OUT.save
for CARTELLA in K*/
do
        echo 'cartella ' $CARTELLA
        cd $CARTELLA
		rm -f *.xml
		iotk convert gkvectors.dat gkvectors.xml
		iotk convert evc.dat evc.xml
		rm -f *.dat
		cd ..
done
cd ..