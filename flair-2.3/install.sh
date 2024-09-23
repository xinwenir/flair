#!/bin/bash

checkout() {
	#svn co http://svn.cern.ch/guest/flair/trunk/flair .
	#svn co http://svn.cern.ch/guest/flair/trunk/geoviewer geoviewer
	svn co http://svn.cern.ch/guest/flair/tags/release/flair .
	svn co http://svn.cern.ch/guest/flair/tags/release/geoviewer geoviewer
	update
}

install() {
	make install-bin DESTDIR=.
	cd geoviewer
	make -j
	make install DESTDIR=..
	make clean
	cd ..
}

update() {
	svn up
	svn up geoviewer
}

help() {
	echo "flair installation script"
	echo
	echo "syntax: install.sh option [directory]"
	echo
	echo "Options:"
	echo "  -i <dir> | install <dir>	Download svn version and install to (dir)ectory"
	echo "  -u | up | update		Update svn version and recompile"
	echo
	exit
}

case $1 in
	-i|install)
		if [ .$2 = . ]; then
			echo
			echo "ERROR: Please specify directory where to install"
			echo
			help
			exit 1
		fi
		mkdir -p $2
		pushd $2
		checkout
		install
		popd
		;;

	-u|u|up|update)
		update
		install
		;;

	*)
		help
		;;
esac
