#!/bin/bash

docker build -t quay.io/nbarkas/abc-general-container .
echo You can now push the image with "docker push quay.io/nbarkas/abc-general-container:latest"