file=2517287028.genes.faa
name=2517287028.genes
chunk_size="10 20 50 100 250 500 1000 1500 2000 2500 3000 3500 4000 5000 6000 8000 10000"

for i in $chunk_size
do
   #name=`echo "$file" | cut -d'.' -f1`
   #echo $file
   command_str=' ./operon_finder_v.0.3.3.py -i '$file' -o '$name'_'$i' --fragment_size '$i' -t all'
   #command_str='blastn -db ../nt/nt -query '$file' -out '$name'.xml -evalue 1e-10 -outfmt 5 -num_threads 8 -max_target_seqs 100'
   #echo $command_str
   res1=$(date +%s.%N)
   eval $command_str
   res2=$(date +%s.%N)
   elapsed=$( echo $res2 - $res1 | bc )
   echo -e $command_str'\t'$i'\t'$elapsed

##           fi
done
