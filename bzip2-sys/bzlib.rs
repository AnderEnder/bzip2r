use std::ffi::CString;
use std::mem::size_of;
use std::os::raw::c_void;
use std::process::exit;

use crate::compress::BZ2_compressBlock;
use crate::crctable::BZ2_crc32Table;
use crate::private_ffi::{
    bz_stream, default_bzalloc, default_bzfree, EState, BZ_CONFIG_ERROR, BZ_MEM_ERROR,
    BZ_M_FINISHING, BZ_M_FLUSHING, BZ_M_RUNNING, BZ_N_OVERSHOOT, BZ_OK, BZ_PARAM_ERROR, BZ_S_INPUT,
    BZ_S_OUTPUT,
};
use std::slice::from_raw_parts_mut;

const BZ_VERSION: &str = "1.0.8, 13-Jul-2019";
const TRUE: u8 = 1;
const FALSE: u8 = 0;

#[no_mangle]
pub unsafe extern "C" fn BZ2_bzlibVersion() -> *const i8 {
    CString::new(BZ_VERSION).unwrap().as_ptr()
}

pub fn bz2_bzlib_version() -> &'static str {
    BZ_VERSION
}

#[no_mangle]
pub fn BZ2_bz__AssertH__fail(errcode: i32) {
    eprintln!(
        "

bzip2/libbzip2: internal error number {}.
This is a bug in bzip2/libbzip2, {}.
Please report it to: bzip2-devel@sourceware.org.  If this happened
when you were using some program which uses libbzip2 as a
component, you should also report this bug to the author(s)
of that program.  Please make an effort to report this bug;
timely and accurate bug reports eventually lead to higher
quality software.  Thanks.

",
        errcode,
        bz2_bzlib_version()
    );

    if errcode == 1007 {
        eprintln!(
            "
*** A special note about internal error number 1007 ***

Experience suggests that a common cause of i.e. 1007
is unreliable memory or other hardware.  The 1007 assertion
just happens to cross-check the results of huge numbers of
memory reads/writes, and so acts (unintendedly) as a stress
test of your memory system.

I suggest the following: try compressing the file again,
possibly monitoring progress in detail with the -vv flag.

* If the error cannot be reproduced, and/or happens at different
  points in compression, you may have a flaky memory system.
  Try a memory-test program.  I have used Memtest86
  (www.memtest86.com).  At the time of writing it is free (GPLd).
  Memtest86 tests memory much more thorougly than your BIOSs
  power-on test, and may find failures that the BIOS doesn't.

* If the error can be repeatably reproduced, this is a bug in
  bzip2, and I would very much like to hear about it.  Please
  let me know, and, ideally, save a copy of the file causing the
  problem -- without which I will be unable to investigate it.
"
        );
    }

    exit(3);
}

#[no_mangle]
pub extern "C" fn bz_config_ok() -> i32 {
    if size_of::<i32>() != 4 {
        return 0;
    }
    if size_of::<i16>() != 2 {
        return 0;
    }
    if size_of::<i8>() != 1 {
        return 0;
    }
    return 1;
}

#[no_mangle]
pub extern "C" fn prepare_new_block(s: &mut EState) {
    s.nblock = 0;
    s.numZ = 0;
    s.state_out_pos = 0;
    s.blockCRC = 0xffffffff;
    for i in 0..256 {
        s.inUse[i] = 0;
    }
    s.blockNo += 1;
}

#[no_mangle]
pub extern "C" fn init_RL(s: &mut EState) {
    s.state_in_ch = 256;
    s.state_in_len = 0;
}

#[no_mangle]
pub extern "C" fn isempty_RL(s: &EState) -> u8 {
    if s.state_in_ch < 256 && s.state_in_len > 0 {
        0
    } else {
        1
    }
}

fn bz_update_crc(crcVar: u32, cha: u8) -> u32 {
    (crcVar << 8) ^ BZ2_crc32Table[((crcVar >> 24) ^ cha as u32) as usize]
}

#[no_mangle]
pub extern "C" fn add_pair_to_block(s: &mut EState) {
    let block = unsafe { from_raw_parts_mut(s.block, s.nblockMAX as usize) };
    let ch = s.state_in_ch as u8;
    for _ in 0..s.state_in_len {
        s.blockCRC = bz_update_crc(s.blockCRC, ch);
    }
    s.inUse[s.state_in_ch as usize] = 1;
    match s.state_in_len {
        1 => {
            block[s.nblock as usize] = ch;
            s.nblock += 1;
        }
        2 => {
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
        }
        3 => {
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
        }
        _ => {
            s.inUse[(s.state_in_len - 4) as usize] = 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = ch;
            s.nblock += 1;
            block[s.nblock as usize] = (s.state_in_len - 4) as u8;
            s.nblock += 1;
        }
    }
}

#[no_mangle]
pub extern "C" fn flush_RL(s: &mut EState) {
    if s.state_in_ch < 256 {
        add_pair_to_block(s);
    }
    init_RL(s);
}

fn add_char_to_block(zs: &mut EState, zchh: u32) {
    let block = unsafe { from_raw_parts_mut(zs.block, zs.nblockMAX as usize) };
    // fast track the common case
    if zchh != zs.state_in_ch && zs.state_in_len == 1 {
        let ch = zs.state_in_ch as u8;
        zs.blockCRC = bz_update_crc(zs.blockCRC, ch);
        zs.inUse[zs.state_in_ch as usize] = 1;
        block[zs.nblock as usize] = ch;
        zs.nblock += 1;
        zs.state_in_ch = zchh;
    } else
    // general, uncommon cases
    if zchh != zs.state_in_ch || zs.state_in_len == 255 {
        if zs.state_in_ch < 256 {
            add_pair_to_block(zs);
        }
        zs.state_in_ch = zchh;
        zs.state_in_len = 1;
    } else {
        zs.state_in_len += 1;
    }
}

#[no_mangle]
pub extern "C" fn copy_input_until_stop(s: &mut EState) -> u8 {
    let mut progress_in = FALSE;
    let strm = unsafe { s.strm.as_mut() }.unwrap();
    // let next_in = unsafe { (strm.next_in as *mut u8).as_mut() }.unwrap();
    let mut next =
        unsafe { from_raw_parts_mut(strm.next_in, strm.avail_in as usize + 1) }.iter_mut();

    if s.mode == BZ_M_RUNNING as i32 {
        // fast track the common case
        loop {
            //block full?
            if s.nblock >= s.nblockMAX {
                break;
            }
            // no input
            if strm.avail_in == 0 {
                break;
            }
            progress_in = TRUE;
            let add = next.next().unwrap();
            add_char_to_block(s, *add as u32);
            strm.next_in = add;
            strm.avail_in -= 1;
            strm.total_in_lo32 += 1;
            if strm.total_in_lo32 == 0 {
                strm.total_in_hi32 += 1;
            }
        }
    } else {
        // general, uncommon case
        loop {
            // block full?
            if s.nblock >= s.nblockMAX {
                break;
            }
            // no input?
            if strm.avail_in == 0 {
                break;
            }
            // flush/finish end?
            if s.avail_in_expect == 0 {
                break;
            }
            progress_in = TRUE;
            let add = next.next().unwrap();
            add_char_to_block(s, *add as u32);
            // strm.next_in = add;
            strm.avail_in -= 1;
            strm.total_in_lo32 += 1;
            if strm.total_in_lo32 == 0 {
                strm.total_in_hi32 += 1;
            }
            s.avail_in_expect -= 1;
        }
    }
    return progress_in;
}

#[no_mangle]
pub extern "C" fn copy_output_until_stop(s: &mut EState) -> u8 {
    let mut progress_out = FALSE;
    let strm = unsafe { s.strm.as_mut() }.unwrap();
    let mut next =
        unsafe { from_raw_parts_mut(strm.next_out, strm.avail_in as usize + 1) }.iter_mut();
    let zbits = unsafe { from_raw_parts_mut(s.zbits, s.nblockMAX as usize) };

    loop {
        // no output space?
        if strm.avail_out == 0 {
            break;
        }

        // block done?
        if s.state_out_pos >= s.numZ {
            break;
        }

        if let Some(current) = next.next() {
            progress_out = TRUE;
            *current = zbits[s.state_out_pos as usize] as i8;
            s.state_out_pos += 1;
            strm.avail_out -= 1;
            // strm.next_out += 1;

            strm.total_out_lo32 += 1;
            if strm.total_out_lo32 == 0 {
                strm.total_out_hi32 += 1;
            }
        } else {
            break;
        }
    }

    return progress_out;
}

#[no_mangle]
pub extern "C" fn handle_compress(strm: &mut bz_stream) -> u8 {
    let mut progress_in = FALSE;
    let mut progress_out = FALSE;
    let s = unsafe { (strm.state as *mut EState).as_mut() }.unwrap();
    loop {
        if s.state == BZ_S_OUTPUT as i32 {
            progress_out |= copy_output_until_stop(s);
            if s.state_out_pos < s.numZ {
                break;
            }
            if s.mode == BZ_M_FINISHING as i32 && s.avail_in_expect == 0 && isempty_RL(s) > 0 {
                break;
            }
            prepare_new_block(s);
            s.state = BZ_S_INPUT as i32;
            if s.mode == BZ_M_FLUSHING as i32 && s.avail_in_expect == 0 && isempty_RL(s) > 0 {
                break;
            }
        }

        if s.state == BZ_S_INPUT as i32 {
            progress_in |= copy_input_until_stop(s);
            if s.mode != BZ_M_RUNNING as i32 && s.avail_in_expect == 0 {
                flush_RL(s);
                BZ2_compressBlock(s, (s.mode == BZ_M_FINISHING as i32) as u8);
                s.state = BZ_S_OUTPUT as i32;
            } else if s.nblock >= s.nblockMAX {
                BZ2_compressBlock(s, FALSE);
                s.state = BZ_S_OUTPUT as i32;
            } else if unsafe { *s.strm }.avail_in == 0 {
                break;
            }
        }
    }

    if progress_in > 0 || progress_out > 0 {
        1
    } else {
        0
    }
}

#[no_mangle]
pub extern "C" fn BZ2_indexIntoF(indx: i32, cftab: *mut i32) -> i32 {
    let cftab = unsafe { from_raw_parts_mut(cftab, 256) };
    let mut nb = 0;
    let mut na = 256;
    while na - nb != 1 {
        let mid = (nb + na) >> 1;
        if indx >= cftab[mid as usize] {
            nb = mid;
        } else {
            na = mid;
        }
    }
    nb
}

#[no_mangle]
pub extern "C" fn BZ2_bzCompressInit(
    strm: &mut bz_stream,
    blockSize100k: i32,
    verbosity: i32,
    mut workFactor: i32,
) -> i32 {
    if bz_config_ok() == 0 {
        return BZ_CONFIG_ERROR;
    }

    // check strm on null
    if blockSize100k < 1 || blockSize100k > 9 || workFactor < 0 || workFactor > 250 {
        return BZ_PARAM_ERROR;
    }

    if workFactor == 0 {
        workFactor = 30;
    }

    if strm.bzalloc.is_none() {
        strm.bzalloc = Some(default_bzalloc);
    }
    if strm.bzfree.is_none() {
        strm.bzfree = Some(default_bzfree);
    }

    let s_ptr = BZALLOC(strm, size_of::<EState>() as i32);
    if s_ptr.is_null() {
        return BZ_MEM_ERROR;
    }

    let s = unsafe { (s_ptr as *mut EState).as_mut() }.unwrap();

    s.strm = strm;

    s.arr1 = std::ptr::null_mut();
    s.arr2 = std::ptr::null_mut();
    s.ftab = std::ptr::null_mut();

    let n = 100000 * blockSize100k;
    s.arr1 = BZALLOC(strm, n * size_of::<u32>() as i32) as *mut u32;
    s.arr2 = BZALLOC(strm, (n + BZ_N_OVERSHOOT as i32) * size_of::<u32>() as i32) as *mut u32;
    s.ftab = BZALLOC(strm, 65537 * size_of::<u32>() as i32) as *mut u32;

    if s.arr1.is_null() || s.arr2.is_null() || s.ftab.is_null() {
        if !s.arr1.is_null() {
            BZFREE(strm, s.arr1 as *mut c_void);
        }
        if !s.arr2.is_null() {
            BZFREE(strm, s.arr2 as *mut c_void);
        }
        if !s.ftab.is_null() {
            BZFREE(strm, s.ftab as *mut c_void);
        }
        if !s_ptr.is_null() {
            BZFREE(strm, s_ptr as *mut c_void);
        }
        return BZ_MEM_ERROR;
    }

    s.blockNo = 0;
    s.state = BZ_S_INPUT as i32;
    s.mode = BZ_M_RUNNING as i32;
    s.combinedCRC = 0;
    s.blockSize100k = blockSize100k;
    s.nblockMAX = 100000 * blockSize100k - 19;
    s.verbosity = verbosity;
    s.workFactor = workFactor;

    s.block = s.arr2 as *mut u8;
    s.mtfv = s.arr1 as *mut u16;
    s.zbits = std::ptr::null_mut();
    s.ptr = s.arr1 as *mut u32;

    strm.state = s_ptr;
    strm.total_in_lo32 = 0;
    strm.total_in_hi32 = 0;
    strm.total_out_lo32 = 0;
    strm.total_out_hi32 = 0;
    init_RL(s);
    prepare_new_block(s);
    return BZ_OK as i32;
}

fn BZALLOC(strm: &mut bz_stream, nnn: i32) -> *mut c_void {
    let bzalloc = strm.bzalloc.unwrap();
    unsafe { bzalloc(strm.opaque, nnn, 1) }
}

fn BZFREE(strm: &mut bz_stream, ppp: *mut c_void) {
    let bzfree = strm.bzfree.unwrap();
    unsafe { bzfree(strm.opaque, ppp) }
}
