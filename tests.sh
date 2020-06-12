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
build/sapphire gsf 63Cu 0 0 0 || echo "Sapphire gsf_test_0 failed!"
build/sapphire gsf 63Cu 1 0 0 || echo "Sapphire gsf_test_0 failed!"
build/sapphire gsf 63Cu 2 2 2 || echo "Sapphire gsf_test_0 failed!"
build/sapphire gsf 63Cu 3 3 2 || echo "Sapphire gsf_test_0 failed!"
