# 1/5/19: Running this manually.

dir_in=/home/603/cxh603/PhD/PartB/Scripts/Model/exp03/

## Find jobs that didn't finish
grep -l 'Exit Status:        271' *.pbs.o* > ${dir_in}Incomplete_weeks_exp03_Try1.txt

## Modify their pbs scripts
python ${dir_in}Change_pbs_scripts_incomplete_weeks.py

## I'VE MANUALLY SPLIT THE pbs_to_submit.txt FILE IN TWO, SO I DON'T OVERLOAD THE QUEUE.
head -n 100 pbs_to_submit.txt > pbs_to_submit_part1.txt
tail -n +101 pbs_to_submit.txt > pbs_to_submit_part2.txt

## Resubmit them to the queue
while read line;do 
qsub ${dir_in}pbs_scripts/$line; 
echo $line 'submitted to queue';
done < pbs_to_submit_part1.txt
