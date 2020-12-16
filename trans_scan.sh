
rm -f -v ./data/R_L8.dat
#rm -f -v ./data/R_L10.dat
rm -f -v ./data/E*.dat

for h in $(seq 2.0 0.2 5.0)
do
	./trans.exe --nh 10000 --h $h --L 8
	#./trans.exe --nh 1000 --h $h --L 10
done

