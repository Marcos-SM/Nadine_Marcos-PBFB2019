#! /bin/bash

find ./*/ -type f ! -name '*.ab1' -delete
rm Quality.txt bestqs.txt log.txt final.fasta
