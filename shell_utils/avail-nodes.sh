#!/usr/bin/bash

pestat | awk '$9!=28 && $9!=80'
