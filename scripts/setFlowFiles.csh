cp ${1} FemtoList.txt
sed -i -- 's|'${2}'|'${3}'|g' FemtoList.txt
sed -i -- 's|/[^/]*$||' FemtoList.txt
cat FemtoList.txt | xargs -I % cp -r % .

set flowFileList=""
ls */*.ep.root >& /dev/null
if ( ${status} == 0 ) then
  foreach flowFile ( `pwd`/*/*.ep.root )
    set flowFileList="${flowFileList}${flowFile};"
  end
  set flowFileList=`echo ${flowFileList} | sed s'|.$||'` 
endif    

echo $flowFileList > flowFileList.txt