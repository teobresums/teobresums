#!/bin/sh

## function to print error message and exit
fail () {
echo "!!! $0: $1" >&2
exit 1
}

## remove M4 cache files
rm -rf autom4te.cache/
rm -f aclocal.m4
autoreconf --install --force

echo "
==================================================
bootstrap has been run successfully.
Now run './configure' with appropriate options
to configure TEOBResum.
==================================================
"
