#!/bin/bash

default=default

info=$($1 --version 2>/dev/null)
# += because certain compilers return non-0 status with --version
test $? -ne 0 && info+=$($1 -qversion 2>/dev/null) || true
test $? -ne 0 && { echo "unknown compiler"; exit 0; } || true

echo $info | grep -i "intel" > /dev/null && { echo "intel"; exit 0; }
echo $info | grep -i "ibm" > /dev/null && { echo "ibm"; exit 0; }
echo $info | grep -i "pgi" > /dev/null && { echo "pgi"; exit 0; }
echo $info | grep -i "clang" > /dev/null && { echo "clang"; exit 0; }
echo $info | grep -iE "gnu|gcc" > /dev/null && { echo "gnu"; exit 0; }
echo $info | grep -i "cray" > /dev/null && { echo "cray"; exit 0; }

echo $default
