#! /usr/bin/env bash

function message {
    echo -e "[$(date "+%Y-%m-%dT%H:%M:%S")] $1"
}

export -f message
