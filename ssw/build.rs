extern crate gcc;

fn main() {
    gcc::compile_library("libssw.a", &["src/ssw.c"]);
}
