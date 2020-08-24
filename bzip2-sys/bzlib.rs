use std::ffi::CString;
use std::mem::size_of;
use std::os::raw::c_void;
use std::process::exit;

use crate::compress::{asserth, BZ2_compressBlock};
use crate::crctable::BZ2_crc32Table;
use crate::decompress::{BZ_GET_FAST, BZ_RAND_UPD_MASK};
use crate::private_ffi::{
    bz_stream, ferror, free, fwrite, malloc, unRLE_obuf_to_output_FAST, BZ2_bzWriteClose64,
    BZ2_decompress, DState, EState, __sFILE, BZFILE, BZ_CONFIG_ERROR, BZ_DATA_ERROR, BZ_FINISH,
    BZ_FINISH_OK, BZ_FLUSH, BZ_FLUSH_OK, BZ_IO_ERROR, BZ_MAX_UNUSED, BZ_MEM_ERROR, BZ_M_FINISHING,
    BZ_M_FLUSHING, BZ_M_IDLE, BZ_M_RUNNING, BZ_N_OVERSHOOT, BZ_OK, BZ_PARAM_ERROR, BZ_RUN,
    BZ_RUN_OK, BZ_SEQUENCE_ERROR, BZ_STREAM_END, BZ_S_INPUT, BZ_S_OUTPUT, BZ_X_BLKHDR_1, BZ_X_IDLE,
    BZ_X_MAGIC_1, BZ_X_OUTPUT, EOF,
};
use std::slice::from_raw_parts_mut;

const BZ_VERSION: &str = "1.0.8, 13-Jul-2019";
const TRUE: u8 = 1;
const FALSE: u8 = 0;

#[repr(C)]
pub struct bzFile {
    handle: *mut libc::FILE,
    buf: [i8; BZ_MAX_UNUSED as usize],
    bufN: i32,
    writing: u8,
    strm: bz_stream,
    lastErr: i32,
    initialisedOk: u8,
}

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

pub fn BZALLOC(strm: &mut bz_stream, nnn: i32) -> *mut c_void {
    let bzalloc = strm.bzalloc.unwrap();
    unsafe { bzalloc(strm.opaque, nnn, 1) }
}

fn BZFREE(strm: &mut bz_stream, ppp: *mut c_void) {
    let bzfree = strm.bzfree.unwrap();
    unsafe { bzfree(strm.opaque, ppp) }
}

#[no_mangle]
pub extern "C" fn default_bzalloc(_opaque: *mut c_void, items: i32, size: i32) -> *mut c_void {
    unsafe { libc::malloc((items * size) as usize) }
}

#[no_mangle]
pub extern "C" fn default_bzfree(_opaque: *mut c_void, addr: *mut c_void) {
    if !addr.is_null() {
        unsafe { libc::free(addr) };
    }
}

#[no_mangle]
pub extern "C" fn BZ2_bzCompress(strm: *mut bz_stream, action: i32) -> i32 {
    //    Bool progress;
    //    EState *s;
    if strm.is_null() {
        return BZ_PARAM_ERROR as i32;
    }
    let strm = unsafe { strm.as_mut() }.unwrap();
    // let s = unsafe { (s_ptr as *mut EState).as_mut() }.unwrap();

    let s_ptr = strm.state as *mut EState;
    if s_ptr.is_null() {
        return BZ_PARAM_ERROR as i32;
    }

    let s = unsafe { s_ptr.as_mut() }.unwrap();
    let s_strm = unsafe { s.strm.as_mut() }.unwrap();

    if s.strm != strm {
        return BZ_PARAM_ERROR as i32;
    }

    loop {
        match s.mode as u32 {
            BZ_M_IDLE => return BZ_SEQUENCE_ERROR as i32,

            BZ_M_RUNNING => {
                if action == BZ_RUN as i32 {
                    let progress = handle_compress(strm);
                    if progress > 0 {
                        return BZ_RUN_OK as i32;
                    } else {
                        return BZ_PARAM_ERROR as i32;
                    };
                } else if action == BZ_FLUSH as i32 {
                    s.avail_in_expect = strm.avail_in;
                    s.mode = BZ_M_FLUSHING as i32;
                    continue;
                } else if action == BZ_FINISH as i32 {
                    s.avail_in_expect = strm.avail_in;
                    s.mode = BZ_M_FINISHING as i32;
                    continue;
                } else {
                    return BZ_PARAM_ERROR;
                }
            }

            BZ_M_FLUSHING => {
                if action != BZ_FLUSH as i32 {
                    return BZ_SEQUENCE_ERROR;
                }
                if s.avail_in_expect != s_strm.avail_in as u32 {
                    return BZ_SEQUENCE_ERROR;
                }
                let _progress = handle_compress(strm);
                if s.avail_in_expect > 0 || isempty_RL(s) == 0 || s.state_out_pos < s.numZ {
                    return BZ_FLUSH_OK as i32;
                }
                s.mode = BZ_M_RUNNING as i32;
                return BZ_RUN_OK as i32;
            }

            BZ_M_FINISHING => {
                if action != BZ_FINISH as i32 {
                    return BZ_SEQUENCE_ERROR;
                }
                if s.avail_in_expect != s_strm.avail_in {
                    return BZ_SEQUENCE_ERROR;
                }
                let progress = handle_compress(strm);
                if progress == 0 {
                    return BZ_SEQUENCE_ERROR;
                }
                if s.avail_in_expect > 0 || isempty_RL(s) == 0 || s.state_out_pos < s.numZ {
                    return BZ_FINISH_OK as i32;
                }
                s.mode = BZ_M_IDLE as i32;
                return BZ_STREAM_END as i32;
            }
            _ => return BZ_OK as i32,
        }
    }
    // return BZ_OK as i32; // not reached
}

#[no_mangle]
pub extern "C" fn BZ2_bzCompressEnd(strm: *mut bz_stream) -> i32 {
    if strm.is_null() {
        return BZ_PARAM_ERROR;
    }

    let strm = unsafe { strm.as_mut() }.unwrap();

    let s = strm.state as *mut EState;
    if s.is_null() {
        return BZ_PARAM_ERROR;
    }

    let s = unsafe { s.as_mut() }.unwrap();

    if s.strm != strm {
        return BZ_PARAM_ERROR;
    }

    if !s.arr1.is_null() {
        BZFREE(strm, s.arr1 as *mut c_void);
    }

    if !s.arr2.is_null() {
        BZFREE(strm, s.arr2 as *mut c_void);
    }

    if s.ftab.is_null() {
        BZFREE(strm, s.ftab as *mut c_void);
    }

    BZFREE(strm, strm.state);

    strm.state = std::ptr::null_mut();

    return BZ_OK as i32;
}

#[no_mangle]
pub extern "C" fn BZ2_bzDecompressInit(strm: *mut bz_stream, verbosity: i32, small: i32) -> i32 {
    if !bz_config_ok() > 0 {
        return BZ_CONFIG_ERROR;
    };

    if strm.is_null() {
        return BZ_PARAM_ERROR;
    }
    if small != 0 && small != 1 {
        return BZ_PARAM_ERROR;
    }
    if verbosity < 0 || verbosity > 4 {
        return BZ_PARAM_ERROR;
    }

    let strm = unsafe { strm.as_mut() }.unwrap();

    if strm.bzalloc.is_none() {
        strm.bzalloc = Some(default_bzalloc);
    }
    if strm.bzfree.is_none() {
        strm.bzfree = Some(default_bzfree);
    }

    let s_raw = BZALLOC(strm, size_of::<DState>() as i32) as *mut DState;
    if s_raw.is_null() {
        return BZ_MEM_ERROR;
    }

    let s = unsafe { s_raw.as_mut() }.unwrap();

    s.strm = strm;
    strm.state = s_raw as *mut c_void;
    s.state = BZ_X_MAGIC_1 as i32;
    s.bsLive = 0;
    s.bsBuff = 0;
    s.calculatedCombinedCRC = 0;
    strm.total_in_lo32 = 0;
    strm.total_in_hi32 = 0;
    strm.total_out_lo32 = 0;
    strm.total_out_hi32 = 0;
    s.smallDecompress = small as u8;
    s.ll4 = std::ptr::null_mut();
    s.ll16 = std::ptr::null_mut();
    s.tt = std::ptr::null_mut();
    s.currBlockNo = 0;
    s.verbosity = verbosity;

    return BZ_OK as i32;
}

fn unRLE_obuf_to_output_FAST2(s: &mut DState) -> u8 {
    //    UChar k1;
    let strm = unsafe { s.strm.as_mut() }.unwrap();

    if s.blockRandomised > 0 {
        loop {
            // try to finish existing run
            loop {
                if strm.avail_out == 0 {
                    return FALSE;
                }
                if s.state_out_len == 0 {
                    break;
                }
                unsafe { *strm.next_out = s.state_out_ch as i8 };
                s.calculatedBlockCRC = bz_update_crc(s.calculatedBlockCRC, s.state_out_ch);

                s.state_out_len -= 1;
                unsafe { *strm.next_out += 1 };
                strm.avail_out -= 1;
                strm.total_out_lo32 += 1;
                if strm.total_out_lo32 == 0 {
                    strm.total_out_hi32 += 1;
                }
            }

            // can a new run be started?
            if s.nblock_used == s.save_nblock + 1 {
                return FALSE;
            }

            // Only caused by corrupt data stream?
            if s.nblock_used > s.save_nblock + 1 {
                return TRUE;
            }

            s.state_out_len = 1;
            s.state_out_ch = s.k0 as u8;
            let mut k1 = BZ_GET_FAST(s);

            BZ_RAND_UPD_MASK(s);

            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 != s.k0 {
                s.k0 = k1;
                continue;
            };

            s.state_out_len = 2;
            k1 = BZ_GET_FAST(s);
            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 != s.k0 {
                s.k0 = k1;
                continue;
            };

            s.state_out_len = 3;
            k1 = BZ_GET_FAST(s);
            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 != s.k0 {
                s.k0 = k1;
                continue;
            };

            k1 = BZ_GET_FAST(s);
            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            s.state_out_len = k1 as i32 + 4;
            s.k0 = BZ_GET_FAST(s);
            BZ_RAND_UPD_MASK(s);
            s.k0 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
        }
    } else {
        // restore
        let mut c_calculatedBlockCRC = s.calculatedBlockCRC as u32;
        let mut c_state_out_ch = s.state_out_ch as u8;
        let mut c_state_out_len = s.state_out_len as i32;
        let mut c_nblock_used = s.nblock_used as i32;
        let mut c_k0 = s.k0 as i32;
        let c_tt = s.tt as *mut u32;
        let c_tt_safe = unsafe { from_raw_parts_mut(c_tt, s.tPos as usize + 1) };
        let mut c_tPos = s.tPos as u32;
        let cs_next_out = strm.next_out as *mut i8;
        let mut cs_avail_out = strm.avail_out as u32;
        let ro_blockSize100k = s.blockSize100k as i32;
        // end restore

        let avail_out_INIT = cs_avail_out;
        let s_save_nblockPP = s.save_nblock + 1;
        //   unsigned int total_out_lo32_old;

        'return_notr: loop {
            loop {
                // try to finish existing run
                if c_state_out_len > 0 {
                    loop {
                        if cs_avail_out == 0 {
                            break 'return_notr;
                        }
                        if c_state_out_len == 1 {
                            break;
                        }
                        unsafe { *cs_next_out = c_state_out_ch as i8 };
                        c_calculatedBlockCRC = bz_update_crc(c_calculatedBlockCRC, c_state_out_ch);

                        c_state_out_len -= 1;
                        unsafe { *cs_next_out += 1 };
                        cs_avail_out -= 1;
                    }
                    //  s_state_out_len_eq_one:

                    if cs_avail_out == 0 {
                        c_state_out_len = 1;
                        break 'return_notr;
                    };
                    unsafe { *cs_next_out = c_state_out_ch as i8 };
                    c_calculatedBlockCRC = bz_update_crc(c_calculatedBlockCRC, c_state_out_ch);

                    unsafe { *cs_next_out += 1 };
                    cs_avail_out -= 1;
                }
                // Only caused by corrupt data stream?
                if c_nblock_used > s_save_nblockPP {
                    return TRUE;
                }

                // can a new run be started?
                if c_nblock_used == s_save_nblockPP {
                    c_state_out_len = 0;
                    break 'return_notr;
                };
                c_state_out_ch = c_k0 as u8;

                let (kt, tPos) = BZ_GET_FAST_C(c_tt_safe, ro_blockSize100k as u32, c_tPos);
                let mut k1 = kt;
                c_tPos = tPos;

                c_nblock_used += 1;
                if k1 != c_k0 {
                    c_k0 = k1;
                    // goto s_state_out_len_eq_one;
                };
                if c_nblock_used == s_save_nblockPP {
                    // goto s_state_out_len_eq_one;
                }

                c_state_out_len = 2;

                let (kt, tPos) = BZ_GET_FAST_C(c_tt_safe, ro_blockSize100k as u32, c_tPos);
                k1 = kt;
                c_tPos = tPos;

                c_nblock_used += 1;
                if c_nblock_used == s_save_nblockPP {
                    continue;
                }
                if k1 != c_k0 {
                    c_k0 = k1;
                    continue;
                };

                c_state_out_len = 3;

                let (kt, tPos) = BZ_GET_FAST_C(c_tt_safe, ro_blockSize100k as u32, c_tPos);
                k1 = kt;
                c_tPos = tPos;

                c_nblock_used += 1;
                if c_nblock_used == s_save_nblockPP {
                    continue;
                }
                if k1 as i32 != c_k0 {
                    c_k0 = k1 as i32;
                    continue;
                };

                let (kt, tPos) = BZ_GET_FAST_C(c_tt_safe, ro_blockSize100k as u32, c_tPos);
                k1 = kt;
                c_tPos = tPos;

                c_nblock_used += 1;
                c_state_out_len = k1 as i32 + 4;

                let (kt, tPos) = BZ_GET_FAST_C(c_tt_safe, ro_blockSize100k as u32, c_tPos);
                c_k0 = kt;
                c_tPos = tPos;
                c_nblock_used += 1;
            }
            break;
        }

        //    return_notr:
        let total_out_lo32_old = strm.total_out_lo32;
        strm.total_out_lo32 += avail_out_INIT - cs_avail_out;
        if strm.total_out_lo32 < total_out_lo32_old {
            strm.total_out_hi32 += 1;
        }

        // save
        s.calculatedBlockCRC = c_calculatedBlockCRC;
        s.state_out_ch = c_state_out_ch;
        s.state_out_len = c_state_out_len;
        s.nblock_used = c_nblock_used;
        s.k0 = c_k0;
        s.tt = c_tt;
        s.tPos = c_tPos;
        unsafe { *strm.next_out = cs_next_out as i8 };
        strm.avail_out = cs_avail_out;
        // end save
    }
    return FALSE;
}

fn BZ_GET_FAST_C(c_tt: &mut [u32], ro_blockSize100k: u32, c_tPos: u32) -> (i32, u32) {
    /* c_tPos is unsigned, hence test < 0 is pointless. */
    if c_tPos >= 100000 * ro_blockSize100k {
        // return TRUE;
    }
    let mut l_tPos = c_tt[c_tPos as usize];
    let cccc = (l_tPos & 0xff) as i32;
    l_tPos >>= 8;
    (cccc, l_tPos)
}

#[no_mangle]
pub extern "C" fn unRLE_obuf_to_output_SMALL(s: &mut DState) -> u8 {
    //    UChar k1;
    let strm = unsafe { s.strm.as_mut() }.unwrap();
    let mut k1 = 0;
    if s.blockRandomised > 0 {
        loop {
            // try to finish existing run
            loop {
                if strm.avail_out == 0 {
                    return FALSE;
                }
                if s.state_out_len == 0 {
                    break;
                }
                unsafe { *strm.next_out = s.state_out_ch as i8 };
                s.calculatedBlockCRC = bz_update_crc(s.calculatedBlockCRC, s.state_out_ch);
                s.state_out_len -= 1;
                unsafe { *strm.next_out += 1 };
                strm.avail_out -= 1;
                strm.total_out_lo32 += 1;
                if strm.total_out_lo32 == 0 {
                    strm.total_out_hi32 += 1;
                }
            }

            // can a new run be started?
            if s.nblock_used == s.save_nblock + 1 {
                return FALSE;
            }

            // Only caused by corrupt data stream?
            if s.nblock_used > s.save_nblock + 1 {
                return TRUE;
            }

            s.state_out_len = 1;
            s.state_out_ch = s.k0 as u8;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            s.state_out_len = 2;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            s.state_out_len = 3;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            BZ_RAND_UPD_MASK(s);
            k1 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
            s.state_out_len = k1 as i32 + 4;
            s.k0 = BZ_GET_SMALL(s, s.nblock_used as u32) as i32;

            BZ_RAND_UPD_MASK(s);
            s.k0 ^= if s.rNToGo == 1 { 1 } else { 0 };
            s.nblock_used += 1;
        }
    } else {
        loop {
            /* try to finish existing run */
            loop {
                if strm.avail_out == 0 {
                    return FALSE;
                }
                if s.state_out_len == 0 {
                    break;
                }
                // *((UChar *)(s.strm.next_out)) = s.state_out_ch;
                // BZ_UPDATE_CRC(s.calculatedBlockCRC, s.state_out_ch);
                s.state_out_len -= 1;
                unsafe { *strm.next_out += 1 };
                strm.avail_out -= 1;
                strm.total_out_lo32 += 1;
                if strm.total_out_lo32 == 0 {
                    strm.total_out_hi32 += 1;
                }
            }

            // can a new run be started?
            if s.nblock_used == s.save_nblock + 1 {
                return FALSE;
            }

            // Only caused by corrupt data stream?
            if s.nblock_used > s.save_nblock + 1 {
                return TRUE;
            }

            s.state_out_len = 1;
            s.state_out_ch = s.k0 as u8;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            s.state_out_len = 2;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            s.state_out_len = 3;
            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            s.nblock_used += 1;
            if s.nblock_used == s.save_nblock + 1 {
                continue;
            }
            if k1 as i32 != s.k0 {
                s.k0 = k1 as i32;
                continue;
            };

            k1 = BZ_GET_SMALL(s, s.nblock_used as u32);

            s.nblock_used += 1;
            s.state_out_len = k1 as i32 + 4;
            s.k0 = BZ_GET_SMALL(s, s.nblock_used as u32) as i32;
            s.nblock_used += 1;
        }
    }
}

fn BZ_GET_SMALL(s: &mut DState, nblock: u32) -> u32 {
    let ll4 = unsafe { from_raw_parts_mut(s.ll4, nblock as usize) };
    let ll16 = unsafe { from_raw_parts_mut(s.ll16, nblock as usize) };

    // c_tPos is unsigned, hence test < 0 is pointless.
    if s.tPos >= 100000 * s.blockSize100k as u32 {
        //return TRUE;
    }
    let cccc = BZ2_indexIntoF(s.tPos as i32, s.cftab.as_mut_ptr()) as u32;
    s.tPos = GET_LL(ll16, ll4, s.tPos as usize);
    cccc
}

fn GET_LL(ll16: &mut [u16], ll4: &mut [u8], i: usize) -> u32 {
    (ll16[i] as u32) | GET_LL4(ll4, i as u8) << 16
}

fn GET_LL4(ll4: &mut [u8], i: u8) -> u32 {
    ((ll4[(i >> 1) as usize] as u32) >> (i << 2) & 0x4) & 0xF
}

#[allow(unreachable_code)]
#[no_mangle]
pub extern "C" fn BZ2_bzDecompress(strm: *mut bz_stream) -> i32 {
    if strm.is_null() {
        return BZ_PARAM_ERROR;
    }
    let strm = unsafe { strm.as_mut() }.unwrap();

    if strm.state.is_null() {
        return BZ_PARAM_ERROR;
    }

    let s = unsafe { (strm.state as *mut DState).as_mut() }.unwrap();

    if s.strm != strm {
        return BZ_PARAM_ERROR;
    }

    loop {
        if s.state == BZ_X_IDLE as i32 {
            return BZ_SEQUENCE_ERROR;
        }
        if s.state == BZ_X_OUTPUT as i32 {
            let corrupt = if s.smallDecompress > 0 {
                unRLE_obuf_to_output_SMALL(s)
            } else {
                unsafe { unRLE_obuf_to_output_FAST(s) }
            };
            if corrupt > 0 {
                return BZ_DATA_ERROR;
            }
            if s.nblock_used == s.save_nblock + 1 && s.state_out_len == 0 {
                s.calculatedBlockCRC = BZ_FINALISE_CRC(s.calculatedBlockCRC);
                if s.verbosity >= 3 {
                    println!(
                        " {{0x{:08x}, 0x{:08x}}}",
                        s.storedBlockCRC, s.calculatedBlockCRC
                    );
                }
                if s.verbosity >= 2 {
                    println!("]");
                }
                if s.calculatedBlockCRC != s.storedBlockCRC {
                    return BZ_DATA_ERROR;
                }
                s.calculatedCombinedCRC =
                    (s.calculatedCombinedCRC << 1) | (s.calculatedCombinedCRC >> 31);
                s.calculatedCombinedCRC ^= s.calculatedBlockCRC;
                s.state = BZ_X_BLKHDR_1 as i32;
            } else {
                return BZ_OK as i32;
            }
        }
        if s.state >= BZ_X_MAGIC_1 as i32 {
            let r = unsafe { BZ2_decompress(s) };
            if r == BZ_STREAM_END as i32 {
                if s.verbosity >= 3 {
                    println!(
                        "\n    combined CRCs: stored = 0x{:08x}, computed = 0x{:08x}",
                        s.storedCombinedCRC, s.calculatedCombinedCRC,
                    );
                }
                if s.calculatedCombinedCRC != s.storedCombinedCRC {
                    return BZ_DATA_ERROR;
                }
                return r;
            }
            if s.state != BZ_X_OUTPUT as i32 {
                return r;
            }
        }
    }

    asserth(false, 6001);

    // NOTREACHED
    return 0;
}

fn BZ_FINALISE_CRC(crcVar: u32) -> u32 {
    !crcVar
}

#[no_mangle]
pub extern "C" fn BZ2_bzDecompressEnd(strm: *mut bz_stream) -> i32 {
    if strm.is_null() {
        return BZ_PARAM_ERROR;
    }
    let strm = unsafe { strm.as_mut() }.unwrap();

    if strm.state.is_null() {
        return BZ_PARAM_ERROR;
    }
    let s = unsafe { (strm.state as *mut DState).as_mut() }.unwrap();

    if s.strm != strm {
        return BZ_PARAM_ERROR;
    }

    if !s.tt.is_null() {
        BZFREE(strm, s.tt as *mut c_void);
    }
    if !s.ll16.is_null() {
        BZFREE(strm, s.ll16 as *mut c_void);
    }
    if !s.ll4.is_null() {
        BZFREE(strm, s.ll4 as *mut c_void);
    }

    BZFREE(strm, strm.state);
    strm.state = std::ptr::null_mut();

    return BZ_OK as i32;
}

#[no_mangle]
pub extern "C" fn myfeof(f: *mut libc::FILE) -> u8 {
    let c = unsafe { libc::fgetc(f) };
    if c == EOF {
        return TRUE;
    }
    unsafe { libc::ungetc(c, f) };
    return FALSE;
}

#[no_mangle]
pub extern "C" fn BZ2_bzWriteOpen(
    bzerror: *mut i32,
    f: *mut libc::FILE,
    blockSize100k: i32,
    verbosity: i32,
    mut workFactor: i32,
) -> *mut BZFILE {
    let bzf = std::ptr::null_mut() as *mut bzFile;

    let null = std::ptr::null_mut();

    BZ_SETERR(bzf, bzerror, BZ_OK as i32);

    if f.is_null()
        || (blockSize100k < 1 || blockSize100k > 9)
        || (workFactor < 0 || workFactor > 250)
        || (verbosity < 0 || verbosity > 4)
    {
        BZ_SETERR(bzf, bzerror, BZ_PARAM_ERROR);
        return null;
    };

    if unsafe { libc::ferror(f) > 0 } {
        BZ_SETERR(bzf, bzerror, BZ_IO_ERROR);
        return null;
    };

    let bzf_raw = unsafe { malloc(size_of::<bzFile>() as u64) as *mut bzFile };
    if bzf_raw.is_null() {
        BZ_SETERR(bzf, bzerror, BZ_MEM_ERROR);
        return null;
    };

    let bzf = unsafe { (bzf_raw as *mut bzFile).as_mut() }.unwrap();

    BZ_SETERR(bzf, bzerror, BZ_OK as i32);
    bzf.initialisedOk = FALSE;
    bzf.bufN = 0;
    bzf.handle = f;
    bzf.writing = TRUE;
    bzf.strm.bzalloc = None;
    bzf.strm.bzfree = None;
    bzf.strm.opaque = null;

    if workFactor == 0 {
        workFactor = 30;
    }
    let ret = BZ2_bzCompressInit(&mut bzf.strm, blockSize100k, verbosity, workFactor);
    if ret != BZ_OK as i32 {
        BZ_SETERR(bzf, bzerror, ret);
        unsafe {
            free(bzf_raw as *mut c_void);
        }
        return null;
    };

    bzf.strm.avail_in = 0;
    bzf.initialisedOk = TRUE;
    return bzf_raw as *mut BZFILE;
}

fn BZ_SETERR(bzf: *mut bzFile, bzerror: *mut i32, eee: i32) {
    if let Some(bzerror) = unsafe { bzerror.as_mut() } {
        *bzerror = eee;
    }
    if let Some(bzf) = unsafe { bzf.as_mut() } {
        bzf.lastErr = eee;
    }
}

#[no_mangle]
pub extern "C" fn BZ2_bzWrite(bzerror: *mut i32, b: *mut BZFILE, buf: *mut c_void, len: i32) {
    //    Int32 n, n2, ret;
    let bzf = b as *mut bzFile;

    BZ_SETERR(bzf, bzerror, BZ_OK as i32);
    if bzf.is_null() || buf.is_null() || len < 0 {
        BZ_SETERR(bzf, bzerror, BZ_PARAM_ERROR);
        return;
    };

    let bzf = unsafe { bzf.as_mut() }.unwrap();

    if !(bzf.writing > 0) {
        BZ_SETERR(bzf, bzerror, BZ_SEQUENCE_ERROR);
        return;
    };
    if unsafe { ferror(bzf.handle as *mut __sFILE) } > 0 {
        BZ_SETERR(bzf, bzerror, BZ_IO_ERROR);
        return;
    };

    if len == 0 {
        BZ_SETERR(bzf, bzerror, BZ_OK as i32);
        return;
    };

    bzf.strm.avail_in = len as u32;
    bzf.strm.next_in = buf as *mut i8;

    loop {
        bzf.strm.avail_out = BZ_MAX_UNUSED;
        bzf.strm.next_out = bzf.buf.as_mut_ptr();
        let ret = BZ2_bzCompress(&mut bzf.strm, BZ_RUN as i32);
        if ret != BZ_RUN_OK as i32 {
            BZ_SETERR(bzf, bzerror, ret);
            return;
        };

        if bzf.strm.avail_out < BZ_MAX_UNUSED {
            let n = BZ_MAX_UNUSED - bzf.strm.avail_out;
            let n2 = unsafe {
                fwrite(
                    bzf.buf.as_mut_ptr() as *mut c_void,
                    size_of::<u8>() as u64,
                    n as u64,
                    bzf.handle as *mut __sFILE,
                )
            };
            if n as u64 != n2 || unsafe { ferror(bzf.handle as *mut __sFILE) } > 0 {
                BZ_SETERR(bzf, bzerror, BZ_IO_ERROR);
                return;
            };
        }

        if bzf.strm.avail_in == 0 {
            BZ_SETERR(bzf, bzerror, BZ_OK as i32);
            return;
        };
    }
}

#[no_mangle]
pub extern "C" fn BZ2_bzWriteClose(
    bzerror: *mut i32,
    b: *mut BZFILE,
    abandon: i32,
    nbytes_in: *mut u32,
    nbytes_out: *mut u32,
) {
    let null = std::ptr::null_mut();
    unsafe { BZ2_bzWriteClose64(bzerror, b, abandon, nbytes_in, null, nbytes_out, null) };
}
