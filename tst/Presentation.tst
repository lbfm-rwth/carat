#!/bin/sh
set -e

echo "### Test Presentation"
../bin/Presentation grp1.dat
echo "### Presentation return code $?"
