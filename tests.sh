#!/usr/bin/bash
######################################################
## @file tests.sh
## @author Philipp Scholz <pscholz@outlook.de>
## @date 2020-05-09
## @brief Script to test Sapphire
######################################################

build/sapphire || echo "Sapphire test failed!"
build/sapphire reaction || echo "Sapphire reaction test failed!"
build/sapphire decayer || echo "Sapphire decayer test failed!"
build/sapphire random || echo "Sapphire random test failed!"
build/sapphire old || echo "Sapphire old test failed!"
build/sapphire template || echo "Sapphire template test failed!"
build/sapphire help || echo "Sapphire help test failed!"
build/sapphire reaction examples/test_reaction.ini || echo "Sapphire reaction examples/test_reaction.ini failed!"
build/sapphire reaction examples/test_reaction2.ini || echo "Sapphire reaction examples/test_reaction2.ini failed!"
build/sapphire reaction examples/test_reaction3.ini || echo "Sapphire reaction examples/test_reaction3.ini failed!"
build/sapphire decayer examples/test_decayer.ini || echo "Sapphire decayer examples/test_decayer.ini failed!"
build/sapphire random examples/test_random.ini || echo "Sapphire random examples/test_random.ini failed!"