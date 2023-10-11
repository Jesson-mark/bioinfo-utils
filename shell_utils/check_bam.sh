bamfile=$1

right_chars="1f8b08040000000000ff0600424302001b0003000000000000000000"
last_28_chars=$(tail -c 28 $bamfile| xxd -p)
if [ $last_28_chars == $right_chars ] ; then
	echo "file is complete"
else
	echo "file isn't complete"
fi

