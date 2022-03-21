cd /s/chris/grow
rm ./bin/conductor_mpi.exe
echo 'delete'
make
cd /s/chris/sim
qsub job.pbs && echo 'DONE' && qsub ppost.pbs && qsub mov.pbs
echo 'DONE'
wait
