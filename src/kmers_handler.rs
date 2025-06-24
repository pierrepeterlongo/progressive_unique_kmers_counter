use bitnuc::{as_2bit, from_2bit, from_2bit_alloc};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Pack a sequence into a u64
    let packed = as_2bit(b"ACGT")?;
    assert_eq!(packed, 0b11100100);

    // Unpack back to a sequence using a reusable buffer
    let mut unpacked = Vec::new();
    from_2bit(packed, 4, &mut unpacked)?;
    assert_eq!(&unpacked, b"ACGT");
    unpacked.clear();

    // Unpack back to a sequence with a reallocation
    let unpacked = from_2bit_alloc(packed, 4)?;
    assert_eq!(&unpacked, b"ACGT");

    Ok(())
}