#!/usr/bin/bash
######################################################
## @file tests.sh
## @author Philipp Scholz <pscholz@outlook.de>
## @date 2020-05-09
## @brief Script to test Sapphire
######################################################

#Basic functionality tests
build/sapphire || echo "Sapphire test failed!"
build/sapphire reaction || echo "Sapphire reaction test failed!"
build/sapphire decayer || echo "Sapphire decayer test failed!"
build/sapphire random || echo "Sapphire random test failed!"
build/sapphire gsf || echo "Sapphire gsf test failed!"
build/sapphire nld || echo "Sapphire nld test failed!"
build/sapphire old || echo "Sapphire old test failed!"
build/sapphire template || echo "Sapphire template test failed!"
build/sapphire help || echo "Sapphire help test failed!"

#Testing of gsf implementations
build/sapphire gsf 63Cu 0 0 0 || echo "Sapphire gsf_test_0 failed!"
build/sapphire gsf 63Cu 1 0 0 || echo "Sapphire gsf_test_1 failed!"
build/sapphire gsf 63Cu 2 2 2 || echo "Sapphire gsf_test_2 failed!"
build/sapphire gsf 63Cu 3 3 2 || echo "Sapphire gsf_test_3 failed!"

#Testing of nld implementations
build/sapphire nld 63Cu 0 || echo "Sapphire nld_test_0 failed!"
build/sapphire nld 64Cu 0 || echo "Sapphire nld_test_1 failed!"
build/sapphire nld 63Cu 1 || echo "Sapphire nld_test_2 failed!"
build/sapphire nld 64Cu 1 || echo "Sapphire nld_test_3 failed!"

#Decayer Module Test
build/sapphire decayer examples/Decay_92Mo/decayer.ini || echo "Sapphire decayer_test_0 failed!"
build/sapphire decayer examples/Decay_70Ge/decayer.ini || echo "Sapphire decayer_test_1 failed!"

#Reaction Module Test
build/sapphire reaction examples/Reaction_N50/reaction.ini || echo "Sapphire reaction_test_0 failed!"
