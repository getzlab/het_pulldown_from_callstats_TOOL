#!/bin/bash

docker build -t getzlab/het_pulldown_from_callstats:v1 .

docker tag getzlab/het_pulldown_from_callstats:v1 gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:v1
docker tag getzlab/het_pulldown_from_callstats:v1 gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:latest

docker push gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:v1
docker push gcr.io/broad-getzlab-workflows/het_pulldown_from_callstats:latest
