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
build/sapphire random examples/test_random.ini || echo "Sapphire random examples/test_random.ini failed!"
build/sapphire old 25Mg+p --average-rad-width || echo "Sapphire old 25Mg+p --average-rad-width failed!"
build/sapphire decayer examples/Decay_92Mo/decayer.ini || echo "Example Decayer_92Mo failed."
build/sapphire reaction examples/Reaction_N50/reaction.ini || echo "Example Reaction_N50 failed."
