#!/usr/bin/env bash

ps aux | grep 'run[.]py' | sed 's/.*run[.]py//'
