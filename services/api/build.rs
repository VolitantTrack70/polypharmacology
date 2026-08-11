//! Compiles the chemworker protobuf definitions into Rust.
//!
//! Only the client is generated -- the server side of this contract is the
//! Python worker, not us. Requires `protoc` on PATH.

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proto = "../chemworker/proto/chemworker.proto";

    tonic_build::configure()
        .build_server(false)
        .build_client(true)
        .compile_protos(&[proto], &["../chemworker/proto"])?;

    println!("cargo:rerun-if-changed={proto}");
    Ok(())
}
