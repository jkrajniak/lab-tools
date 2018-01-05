mkdir analyze
cd analyze
for f in ../*.txt; do
    ln -s $f `basename $f`.xvg
done
ln -s ../../../../../../analysis_code/*.py .

sysid=$(echo $PWD | cut -f11 -d'/')
cl=$(echo $PWD | cut -f14 -d'/' | tr -d 'cl_')
axis=$(echo $PWD | cut -f15 -d'/' | tr -d 'pull_step_')

echo $sysid $cl $axis

python data_energy.py ${sysid}_${axis}_${cl}_energy
python data_strain.py ${sysid}_${axis}_${cl}
