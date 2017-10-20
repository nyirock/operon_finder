res1=$(date +%s.%N)
sleep 1
res2=$(date +%s.%N)
echo "Start time: $res1"
echo "Stop time:  $res2"
echo "Elapsed:    $(( res2 - res1 ))"
