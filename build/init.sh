#!/bin/bash
# Check for sudo
if [ "$EUID" -ne 0 ]
  then echo "Needs to be run as root user"
  exit
fi
cpupower frequency-set -g performance
