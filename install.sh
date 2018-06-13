#!/bin/bash -e -x
CC=${PREFIX}/bin/gcc
cd mtsv/ext
#cargo test
cargo build --release
