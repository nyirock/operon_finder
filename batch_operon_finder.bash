
for file in *.faa
do
   #name=`echo "$file" | cut -d'.' -f1`
   echo $file
   command_str='./operon_finder_v.0.3.2.py -i '$file
   #command_str='blastn -db ../nt/nt -query '$file' -out '$name'.xml -evalue 1e-10 -outfmt 5 -num_threads 8 -max_target_seqs 100'
   echo $command_str
   eval $command_str

##           fi
done
