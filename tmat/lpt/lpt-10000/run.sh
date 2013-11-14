mkdir -p logs
accelerator interact --set_beta=0.0
for i in `seq 1 30`
do
    if [ $i -eq 15 ]
    then
        accelerator interact --set_beta=2.0
    fi

    accelerator TMat &>logs/round-$i-para-0 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-1 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-2 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-3 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-4 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-5 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-6 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-7 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-8 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-9 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-10 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-11 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-12 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-13 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-14 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-15 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-16 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-17 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-18 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-19 &
    sleep 10

    wait

    wait

    accelerator model &>logs/model-$i.log
done
