NP=3
echo ${NP}
mpirun -np ${NP} src/mesher.py 300 400
mpirun -np ${NP} src/run.py  300 400
mpirun -np ${NP} src/plot.py 300 400
