#include "gtest/gtest.h"

extern void Initialize();

int main(int argc, char **argv){
    ::testing::InitGoogleTest(&argc, argv);
    Initialize();
    return RUN_ALL_TESTS();
}