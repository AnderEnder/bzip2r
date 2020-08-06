use crate::private_ffi::{BZ_MAX_ALPHA_SIZE, BZ_MAX_CODE_LEN};
use std::slice::{from_raw_parts, from_raw_parts_mut};

fn asserth(cond: bool, errcode: i32) {
    if !cond {
        let version = "0.1";
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
            errcode, version
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
  problem -- without which I will be unable to investigate it."
            );
        }
    }
}

fn upheap(z: usize, heap: &mut [i32], weight: &[i32]) {
    let mut zz = z;
    let tmp = heap[zz] as usize;
    while weight[tmp] < weight[heap[zz >> 1] as usize] {
        heap[zz] = heap[zz >> 1];
        zz >>= 1;
    }
    heap[zz] = tmp as i32;
}

fn downheap(z: usize, nHeap: usize, heap: &mut [i32], weight: &[i32]) {
    let mut zz = z;
    let tmp = heap[zz] as usize;
    loop {
        let mut yy = zz << 1;
        if yy > nHeap {
            break;
        }

        if yy < nHeap && weight[heap[yy + 1] as usize] < weight[heap[yy] as usize] {
            yy += 1;
        }

        if weight[tmp] < weight[heap[yy] as usize] {
            break;
        }

        heap[zz] = heap[yy];
        zz = yy;
    }
    heap[zz] = tmp as i32;
}

#[allow(overflowing_literals)]
fn weightof(zz0: i32) -> i32 {
    zz0 & 0xffffff00
}

fn depthof(zz1: i32) -> i32 {
    zz1 & 0x000000ff
}

fn mymax(zz2: i32, zz3: i32) -> i32 {
    if zz2 > zz3 {
        zz2
    } else {
        zz3
    }
}

fn addweights(zw1: i32, zw2: i32) -> i32 {
    (weightof(zw1) + weightof(zw2)) | (1 + mymax(depthof(zw1), depthof(zw2)))
}

#[no_mangle]
pub extern "C" fn BZ2_hbMakeCodeLengths(len: *mut u8, freq: *mut i32, alphaSize: i32, maxLen: i32) {
    let len = unsafe { from_raw_parts_mut(len, alphaSize as usize) };
    let freq = unsafe { from_raw_parts_mut(freq, alphaSize as usize) };
    bz2_hb_make_code_lengths(len, freq, alphaSize, maxLen);
}

fn bz2_hb_make_code_lengths(len: &mut [u8], freq: &mut [i32], alphaSize: i32, maxLen: i32) {
    let alpha = alphaSize as usize;
    let max = maxLen as usize;
    let mut heap = vec![0_i32; BZ_MAX_ALPHA_SIZE as usize + 2];
    let mut weight = vec![0_i32; BZ_MAX_ALPHA_SIZE as usize + 2];
    let mut parent = vec![0_i32; BZ_MAX_ALPHA_SIZE as usize + 2];

    for i in 0..alpha {
        weight[i + 1] = if freq[i] == 0 { 1 } else { freq[i] } << 8;
    }

    // let weight = freq
    //     .iter()
    //     .map(|x| if *x == 0 { 1_i32 } else { *x })
    //     .collect::<Vec<_>>();

    loop {
        let mut nNodes = alpha;
        let mut nHeap = 0;

        heap[0] = 0;
        weight[0] = 0;
        parent[0] = -2;

        for i in 1..alpha + 1 {
            parent[i] = -1;
            nHeap += 1;
            heap[nHeap] = i as i32;
            upheap(nHeap, heap.as_mut_slice(), &weight);
        }

        asserth(nHeap < (BZ_MAX_ALPHA_SIZE as usize + 2), 2001);

        while nHeap > 1 {
            let n1 = heap[1] as usize;
            heap[1] = heap[nHeap];
            nHeap -= 1;
            downheap(1, nHeap, heap.as_mut_slice(), &weight);
            let n2 = heap[1] as usize;
            heap[1] = heap[nHeap];
            nHeap -= 1;
            downheap(1, nHeap, heap.as_mut_slice(), &weight);
            nNodes += 1;
            parent[n1] = nNodes as i32;
            parent[n2] = nNodes as i32;
            weight[nNodes] = addweights(weight[n1], weight[n2]) as i32;
            parent[nNodes] = -1;
            nHeap += 1;
            heap[nHeap] = nNodes as i32;
            upheap(nHeap, heap.as_mut_slice(), &weight);
        }

        asserth(nNodes < (BZ_MAX_ALPHA_SIZE as usize * 2), 2002);

        let mut tooLong = false;
        for i in 1..alpha + 1 {
            let mut j = 0_usize;
            let mut k = i;
            while parent[k] >= 0 {
                k = parent[k] as usize;
                j += 1;
            }
            len[i - 1] = j as u8;
            if j > max {
                tooLong = true;
            }
        }

        if !tooLong {
            break;
        }

        for i in 1..alpha + 1 {
            let mut j = weight[i] >> 8;
            j = 1 + (j / 2);
            weight[i] = j << 8;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn BZ2_hbAssignCodes(
    code: *mut i32,
    length: *mut u8,
    minLen: i32,
    maxLen: i32,
    alphaSize: i32,
) {
    let code = from_raw_parts_mut(code, alphaSize as usize);
    let length = from_raw_parts(length, alphaSize as usize);
    bz2_hb_assign_codes(code, length, minLen, maxLen);
}

fn bz2_hb_assign_codes(code: &mut [i32], length: &[u8], minLen: i32, maxLen: i32) {
    let mut vec: i32 = 0;
    for n in minLen..(maxLen + 1) {
        for (length, code) in length.iter().zip(code.iter_mut()) {
            if *length as i32 == n {
                *code = vec;
                vec += 1;
            }
        }
        vec <<= 1;
    }
}

fn bz2_hb_assign_codes2(code: &mut [i32], length: &[u8], minLen: i32, maxLen: i32, alphaSize: i32) {
    let mut n: i32 = minLen;
    let mut vec: i32 = 0;

    while n <= maxLen {
        let mut i = 0;
        while i < alphaSize {
            if length[i as usize] as i32 == n {
                code[i as usize] = vec;
                vec += 1;
            }
            i += 1
        }
        vec <<= 1;
        n += 1
    }
}

#[no_mangle]
pub unsafe extern "C" fn BZ2_hbCreateDecodeTables(
    limit: *mut i32,
    base: *mut i32,
    perm: *mut i32,
    length: *mut u8,
    minLen: i32,
    maxLen: i32,
    alphaSize: i32,
) {
    let limit = from_raw_parts_mut(limit, BZ_MAX_CODE_LEN as usize);
    let base = from_raw_parts_mut(base, BZ_MAX_CODE_LEN as usize);
    let perm = from_raw_parts_mut(perm, alphaSize as usize);
    let length = from_raw_parts_mut(length, alphaSize as usize);
    bz2_hb_create_decode_tables(limit, base, perm, length, minLen, maxLen, alphaSize);
}

pub fn bz2_hb_create_decode_tables(
    limit: &mut [i32],
    base: &mut [i32],
    perm: &mut [i32],
    length: &mut [u8],
    minLen: i32,
    maxLen: i32,
    alphaSize: i32,
) {
    let min = minLen as usize;
    let max = maxLen as usize;
    let alpha = alphaSize as usize;

    let mut pp = 0_usize;

    for i in min..max + 1 {
        for j in 0..alpha {
            if length[j] == (i as u8) {
                perm[pp] = j as i32;
                pp += 1;
            }
        }
    }

    // base.iter_mut()
    //     .take(BZ_MAX_CODE_LEN as usize)
    //     .for_each(|x| *x = 0);

    for i in 0..BZ_MAX_CODE_LEN as usize {
        base[i] = 0;
    }

    for i in 0..alpha {
        base[length[i] as usize + 1] += 1;
    }

    for i in 1..BZ_MAX_CODE_LEN as usize {
        base[i] += base[i - 1];
    }

    for i in 0..BZ_MAX_CODE_LEN as usize {
        limit[i] = 0;
    }

    let mut vec = 0;

    for i in min..(max + 1) {
        vec += base[i + 1] - base[i];
        limit[i] = vec - 1;
        vec <<= 1;
    }

    for i in (min + 1)..(max + 1) {
        base[i] = ((limit[i - 1] + 1) << 1) - base[i];
    }
}

pub fn bz2_hb_create_decode_tables2(
    limit: &mut [i32],
    base: &mut [i32],
    perm: &mut [i32],
    length: &mut [u8],
    minLen: i32,
    maxLen: i32,
    alphaSize: i32,
) {
    let mut pp: i32 = 0;
    let mut i: i32 = minLen;
    let mut vec: i32 = 0;

    while i <= maxLen {
        let mut j = 0;
        while j < alphaSize {
            if length[j as usize] as i32 == i {
                perm[pp as usize] = j;
                pp += 1
            }
            j += 1
        }
        i += 1
    }

    i = 0;
    while i < BZ_MAX_CODE_LEN as i32 {
        base[i as usize] = 0;
        i += 1
    }

    i = 0;
    while i < alphaSize {
        base[(length[i as usize] as usize + 1)] += 1;
        i += 1
    }

    i = 1;
    while i < BZ_MAX_CODE_LEN as i32 {
        base[i as usize] += base[(i - 1) as usize];
        i += 1
    }

    i = 0;
    while i < BZ_MAX_CODE_LEN as i32 {
        limit[i as usize] = 0;
        i += 1
    }

    i = minLen;
    while i <= maxLen {
        vec += base[(i + 1) as usize] - base[i as usize];
        limit[i as usize] = vec - 1;
        vec <<= 1 as i32;
        i += 1
    }

    i = minLen + 1;
    while i <= maxLen {
        base[i as usize] = ((limit[(i - 1) as usize] + 1) << 1) - base[i as usize];
        i += 1
    }
}
