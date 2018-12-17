#!/usr/bin/env bash
# Inspired by `tst/testspecial/run_all.sh` from GAP
# (see https://github.com/gap-system/gap/)
#
# CARAT test files must have extension `.tst`.
# Comments and blank lines are ignored. For a test file named
# `testfile.tst`, the correct output of `bash testfile.tst`
# must be placed in the file `testfile.out`.
#
# This script iterates over all `*.tst` files, running each test
# and comparing its observed and expected output.
#
#set -e
retvalue=0
for testfile in *.tst; do
    echo "### Testing ${testfile}"
    name="$(basename -s .tst ${testfile})"
    #bash -e "./${name}.tst" > "${name}.tmp" 2>&1
    bash "./${name}.tst" > "${name}.tmp" 2>&1
    if ! diff -b "${name}.out" "${name}.tmp"; then
        echo "${testfile} failed, see observed/expected output above"
        retvalue=1
    else
        rm -f "${name}.tmp"
    fi
done
exit ${retvalue}

