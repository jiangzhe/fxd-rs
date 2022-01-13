# fxd-rs

Port MySQL decimal implementation in Rust.

The initial version uses 9 4-byte integers for units and 2 bytes for precision/scale.
Plan to reduce memory footprint using const generics.
