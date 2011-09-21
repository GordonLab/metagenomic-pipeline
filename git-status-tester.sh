#!/bin/sh
# Author: Nick Semenkovich <semenko@alum.mit.edu>

gitstatus=`git status -s | grep -v "^ M " | grep -v "^ ?? " `

if [ "$gitstatus" ]
then
    if tty -s
    then
	USER_TTY=`tty`
	echo "\a@@@ WARNING @@@:\nYou are using an outdated version of this software." > $USER_TTY
	echo "There is an updated version available in git.\n" > $USER_TTY
	echo "Please run:\n\t git pull" > $USER_TTY
	sleep 2
    fi
fi

#END
