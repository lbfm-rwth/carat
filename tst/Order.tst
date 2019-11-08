#!/bin/sh
set -e

echo "### Test Order"
../bin/Order grp1.dat
echo "### Order return code $?"
