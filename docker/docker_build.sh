#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REPO_NAME="daisymut/sammyR"
IMAGE="ghcr.io/daisymut/sammyr"
TAG="${1}"

echo "Building ${IMAGE}:${TAG} da ${REPO_NAME}@${TAG}"

docker buildx build \
    --platform linux/amd64,linux/arm64 \
    -f "${SCRIPT_DIR}/sammyR.dockerfile" \
    --build-arg repo_name=${REPO_NAME} \
    --build-arg tag_name=${TAG} \
    -t ${IMAGE}:${TAG} \
    -t ${IMAGE}:latest \
    --push \
    "${SCRIPT_DIR}"
