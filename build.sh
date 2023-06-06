#!/bin/bash

cr="ghcr.io/colossal-compsci"
software=$(awk '/<!-- SOFTWARE/{print $3}' README.md)
version=$(awk '/<!-- VERSION/{print $3}' README.md)
cmd="docker build $@ -t ${software}:${version} -t ${software}:latest -t ${cr}/${software}:${version} -t ${cr}/${software}:latest ."

echo $cmd
read -p "Build? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then eval $cmd; fi
