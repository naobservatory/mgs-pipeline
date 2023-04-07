#!/usr/bin/env bash

ps aux | grep -v 'SC[R]EEN' | grep 'run[.]py' | sed 's/.*run[.]py//' | grep .
