#!/bin/sh

##This script makes html  files and underlined thos with tags <:t:>  and <:/t:>

dummy1=dummy1.txt
dummy2=dummy2.txt

dir=$1


for filename in `ls $dir`
do 



new_filename=$filename.html
html_header="<!DOCTYPE html> <html> <body>"
html_heading="<h1>$filename</h1>"
html_end_tag="</p> </body> </html>"

echo $html_header >> $new_filename 
echo $html_heading >> $new_filename 

cat $filename > $dummy1 

####color tags
sed 's/<:t:>/ <font color="red"><b> /g' $dummy1 > $dummy2
sed 's/<:\/t:>/ <\/font><\/b> /g' $dummy2 > $dummy1

###Now put them intp into new html file
cat $dummy1 >>  $new_filename 
echo $html_end_tag >> $new_filename 

rm $dummy1 
rm $dummy2

done   
