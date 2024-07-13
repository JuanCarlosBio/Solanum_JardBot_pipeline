#!/usr/bin/env bash

paste <(ls -1 data/GVCF/* \
          | grep -v '\.idx$' \
          | sed 's/data\/GVCF\///; s/\.gvcf//') \
      <(ls -1 data/GVCF/* \
          | grep -v '\.idx$') 

