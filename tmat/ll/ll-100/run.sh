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

    accelerator TMat &>logs/round-$i-para-20 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-21 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-22 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-23 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-24 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-25 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-26 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-27 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-28 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-29 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-30 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-31 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-32 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-33 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-34 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-35 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-36 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-37 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-38 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-39 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-40 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-41 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-42 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-43 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-44 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-45 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-46 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-47 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-48 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-49 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-50 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-51 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-52 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-53 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-54 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-55 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-56 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-57 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-58 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-59 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-60 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-61 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-62 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-63 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-64 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-65 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-66 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-67 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-68 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-69 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-70 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-71 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-72 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-73 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-74 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-75 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-76 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-77 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-78 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-79 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-80 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-81 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-82 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-83 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-84 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-85 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-86 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-87 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-88 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-89 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-90 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-91 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-92 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-93 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-94 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-95 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-96 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-97 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-98 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-99 &
    sleep 10

    wait

    wait

    accelerator model &>logs/model-$i.log
done
