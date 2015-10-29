while read sampleid
do 
  [ ! -e ./$sampleid ] || rm -rf ./$sampleid
  mkdir $sampleid
  cd $sampleid
  cat ../$2 | grep $sampleid | cut -f 11 | tr ';' ' ' | xargs wget
  cd ..
  echo "\n"
done < $1
