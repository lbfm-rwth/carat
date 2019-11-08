#!/bin/sh
set -e

echo "### Test Ex8-1"
../bin/Aut_grp Ex8_f
echo "### Ex8-1 return code $?"

echo "### Test Ex8-2"
../bin/Bravais_inclusions Ex8_gf
echo "### Ex8-2 return code $?"
