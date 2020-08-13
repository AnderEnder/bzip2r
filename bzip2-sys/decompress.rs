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

fn GET_BITS(s: &mut DState, lll: i32, mut vvv: i32, nnn: i32) {
    if lll > 0 {
        s.state = lll;
        loop {
            if s.bsLive >= nnn {
                let v = (s.bsBuff >> (s.bsLive - nnn)) & ((1 << nnn) - 1);
                s.bsLive -= nnn;
                vvv = v as i32;
                break;
            }
            if (unsafe { *s.strm }.avail_in) == 0 {
                // RETURN(BZ_OK);
            }
            s.bsBuff = (s.bsBuff << 8) | (unsafe { *((*s.strm).next_in as *mut u8) } as u32);
            s.bsLive += 8;
            unsafe {
                *(*s.strm).next_in += 1;
                (*s.strm).avail_in -= 1;
                (*s.strm).total_in_lo32 += 1;
                if (*s.strm).total_in_lo32 == 0 {
                    (*s.strm).total_in_hi32 += 1;
                }
            }
        }
    }
}
