

# This script provide the summary table for unique UMIs, final UMIs, umi collision and insert error
# to run execute the command 
# navigate to the folder contraining  final_ccs_reads, umi_collision, and error_insert, and umi_stats folders
# project: path to a directorty contraining  final_ccs_reads, umi_collision, and error_insert folders


project=$1
gen_type=$2

echo 'sample name,unique umi count prior cutoff,unique umi count post cutoff,umi in final reads,umi collision,umi error,final post curation' > $project"Count_stat.csv"

ls $project | grep ".fasta" | rev | cut -d "." -f2- | rev | while read -r file 
do
  sample=$(echo $file)
  total_count=$(wc -l < $project"umi_stats/"$file"_counts_UMI.txt")
  total_count_post_cutoff=$(wc -l < $project"umi_stats/"$file"_umi_seq.txt")
  final_reads=$(grep ">" $project"final_ccs_reads/"$file"read.fasta" | wc -l)
  umi_file=$project"umi_collision/"$file"collision.txt"
  umi_collision=$([ -f $umi_file ] && wc -l < $umi_file || echo 0)
  error_file=$project"error_insert/"$file"error_insert.txt"
  error_insert=$([ -f $error_file ] && wc -l < $error_file || echo 0)
  final_post_cur=$(grep ">" $project"final_post_curation/"$gen_type"/"$file"_final.fasta" | wc -l)


  echo $sample","$total_count","$total_count_post_cutoff","$final_reads","$umi_collision","$error_insert","$final_post_cur >> $project"Count_stat.csv"
done

echo "Done!"
