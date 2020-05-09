#!/usr/bin/bash
######################################################
## @file tests.sh
## @author Philipp Scholz <pscholz@outlook.de>
## @date 2020-05-09
## @brief Script to test Sapphire
######################################################

build/sapphire
build/sapphire reaction
build/sapphire decayer
build/sapphire random 
build/sapphire old
build/sapphire template
build/sapphire help
build/sapphire reaction examples/test_reaction.ini
build/sapphire decayer examples/test_decayer.ini
build/sapphire decayer examples/test_random.ini
