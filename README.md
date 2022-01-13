# fxd-rs

![build](https://github.com/jiangzhe/fxd-rs/actions/workflows/build.yml/badge.svg)
![codecov](https://codecov.io/gh/jiangzhe/fxd-rs/branch/main/graph/badge.svg?token=8WGYSWBJ23)

Port MySQL decimal implementation in Rust.

The initial version uses 9 4-byte integers for units and 2 bytes for precision/scale.
Plan to reduce memory footprint using const generics.

### License

This project is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or
   https://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or
   https://opensource.org/licenses/MIT)

at your option.
