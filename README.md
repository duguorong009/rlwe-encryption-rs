# RLWE Encryption Scheme  
This repo is for implementing this encryption scheme(https://eprint.iacr.org/2014/725.pdf) with Rustlang.  

Inspired by [this work](https://github.com/jnortiz/RLWE).  

**IMPORTANT**: For testing, we should use the `cargo r --release` command.  
The reason is that there is overflow of `i32` in `knuth-yao` sampling code.  
We need to disable overflow check for successful run of program.  
Also, it is to match the behavior of original cpp code.  

