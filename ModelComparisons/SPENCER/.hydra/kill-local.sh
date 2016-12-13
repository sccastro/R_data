#!/bin/bash

ps=`ps a | grep $1 | grep -v grep`
kill ${ps:0:5}

