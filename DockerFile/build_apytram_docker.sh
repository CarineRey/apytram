#!/bin/bash

set -eo

IMAGE_NAME=apytram
DOCKERFILE_DIR=.
TAG=v1.2-alpha
REPO=carinerey/$IMAGE_NAME:$TAG
docker build --no-cache -t $REPO $DOCKERFILE_DIR
docker push $REPO
