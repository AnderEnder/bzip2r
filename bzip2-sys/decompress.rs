use crate::private_ffi::DState;

#[no_mangle]
pub extern "C" fn makeMaps_d(s: &mut DState) {
    s.nInUse = 0;
    for i in 0..256 {
        if s.inUse[i] > 0 {
            s.seqToUnseq[s.nInUse as usize] = i as u8;
            s.nInUse += 1;
        }
    }
}
