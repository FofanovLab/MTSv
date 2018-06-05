extern crate cc;

fn main() {
cc::Build::new()
    .file("src/ssw.c")
    .compile("libssw.a");
 }