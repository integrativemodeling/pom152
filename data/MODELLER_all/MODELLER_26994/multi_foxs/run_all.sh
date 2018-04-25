nohup ./lyre_run_multi_foxs1.sh > lyre_run_multi_foxs1.log &
sleep 3
nohup ./lyre_run_multi_foxs2.sh > lyre_run_multi_foxs2.log &
sleep 3
nohup ./lyre_run_multi_foxs3.sh > lyre_run_multi_foxs3.log &

#qsub ./chef_run_multi_foxs1.sh
#qsub ./chef_run_multi_foxs2.sh
#qsub ./chef_run_multi_foxs3.sh
