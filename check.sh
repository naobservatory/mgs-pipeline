#!/usr/bin/env bash

CHANGED=0
if [[ $# = 1 ]] && [[ "$1" = "--fix" ]]; then
    if ! black --check . &> /dev/null; then
        black . || exit 1
        CHANGED=1
    fi

    if ! isort --check . &> /dev/null; then
        isort . || exit 1
        CHANGED=1
    fi
elif [[ $# -gt 0 ]]; then
    echo "Usage: $0 [--fix]" > /dev/stderr
    exit 1
fi

# if this fails with "black: command not found" you need to install black
#    python3 -m pip install black
echo Running black to check formatting...
if ! black --check --diff --quiet .; then
    echo "FAIL: formatting"
    exit 1
fi


