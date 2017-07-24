grep -e "aliquot_id\|entity_submitter_id" $1 | sed 's/\(\s*\"entity_submitter_id\": \"\)\([0-9A-Za-z-]*\)",/\2/' | sed 's/\(.*aliquot_id\": \"\)\([A-Z0-9a-z-]*\)\",/\2\t@/' |tr "\n" "\t" |sed 's/@/\n/g' |sed 's/^\s*//'
#
#grep -e "\"submitter_id\"\|entity_submitter_id" $1 |grep -B 1 "entity" |grep submitter |sed 's/\(\s*\"submitter_id\": \"\)\([0-9A-Za-z-]*\)",/\2/' | sed 's/\(.*entity_submitter_id\": \"\)\([A-Z0-9-]*\)\",/\2\t@/' |tr "\n" "\t" |sed 's/@/\n/g' |sed 's/^\s*//' 

