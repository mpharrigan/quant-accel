mkdir -p logs
accelerator interact --set_beta=0.0
for i in `seq 1 30`
do
    if [ $i -eq 15 ]
    then
        accelerator interact --set_beta=2.0
    fi

    accelerator TMat &>logs/round-$i-para-0.log &
    sleep 10

    wait

    accelerator model &>logs/model-$i.log
done
