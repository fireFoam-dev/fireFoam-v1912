#!/bin/bash
# adds appropriate file and rule to wmake rules directory

destination=$WM_DIR/rules/General/
file=version2

if cmp $destination/$file $file >/dev/null 2>&1
    then
    cmp $WM_DIR/rules/General/version2 version2
    echo "./$file and $destination/$file are identical"
else
    cmp $WM_DIR/rules/General/version2 version2 
    if [ -e $file -a -d $destination ]
	then
	echo "installing $file in $destination"
	cp -rp $file $destination
    fi
fi

file=standard
if [ -e  $destination/$file ]
then
    if ! grep -q "fireFoam" "$destination/$file" ;then
	echo "adding rule to $destination/$file"
	echo "#fireFoam version stamping" >> $destination/$file
	echo "include \$(GENERAL_RULES)/version2" >> $destination/$file
    fi    
fi
