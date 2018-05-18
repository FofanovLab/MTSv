#!/bin/bash -e -x
cd ext
cargo update
cargo test
cargo build --release
