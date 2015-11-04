# scripts_public - a set of usefule utilities


# Each script __should__ feature a summary line and some keyword 'tags'


# To find scripts with specific tags use:-
	./tagFinder.R --tags 'tag1,tag2'

To generate the list of script summaries in SCRIPT_SUMMARY.txt, do:-
	 grep -e "---SUMMARY---" *.R | sed s/---SUMMARY---// > SCRIPT_SUMMARY.txt



