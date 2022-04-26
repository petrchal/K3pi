echo $1
grep ENTRIES $1 | awk '{print $2}' | sed 's/=/ /g' | awk '{print $2}'| awk '{ sum += $1 } END { print sum }'
