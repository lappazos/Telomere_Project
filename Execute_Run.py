import os

for i in range(27):
    os.system("sbatch --mem=8g -c1 --time=15:0:0 --wrap 'python3 ../../../Telomere_Project/Telo_Detecter.py "+str(4*i)+" "+str(4*(i+1))+"'")