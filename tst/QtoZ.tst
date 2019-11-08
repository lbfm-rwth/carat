#!/bin/sh
set -e

# Testing numbering of Z-classes
echo "### Test QtoZ"
../bin/QtoZ grp2.dat
echo "### QtoZ return code $?"
